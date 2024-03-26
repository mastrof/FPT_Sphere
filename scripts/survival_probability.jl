using Distributed
@everywhere using DrWatson
@everywhere begin
    @quickactivate :FPT_Sphere
    using Agents, MicrobeAgents
    using DataFrames
end
using CSV
using LsqFit

##
@everywhere function initialize(;
    N = 1000,
    # Φ = 1 - 1e-4, # free volume fraction -> (1-Φ)L^3 = V // for later
    L = 1000.0,
    R = 50.0,
    C₀ = 0.0,
    C₁ = 1.0,
    U = 20.0,
    Δt = 0.1,
)
    D = 3
    periodic = true
    # define L by fixing excluded volume // for later
    # V = 4π/3 * R^3 # μm^3
    # L = (V/(1-Φ))^(1/3) # μm
    space = ContinuousSpace(ntuple(_ -> L, D); periodic)
    properties = Dict(
        :chemoattractant => DiffusiveField(C₀, C₁, R, spacesize(space) ./ 2),
        :last_position => [zero(SVector{D,Float64}) for _ in 1:N]
    )
    model = StandardABM(BrownBerg{D}, space, Δt;
        properties, container=Dict,
        model_step! = encounter!,
    )
    motility = RunTumble(; speed=[U])
    for i in 1:N
        add_agent!(model; motility)
        # after adding agent to model, save its assigned position
        model.last_position[i] = position(model[i])
    end
    model
end

##
@everywhere begin
    mdata = [nagents]
    Δt = 0.1 # s
    T = 2*60*60 # s
    nsteps = round(Int, T/Δt)
    stop(model, s) = s >= nsteps || nagents(model) <= 1 
    parameters = Dict(
        :Δt => [Δt],
        #:Φ => [1-1e-4],
        :L => [1000.0],
        :R => [10.0, 25.0, 50.0],
        :C₁ => [0, 0.1, 0.5, 1.0],
        :U => [15.0, 30.0, 60.0]
    )
end
_, mdf = paramscan(parameters, initialize; mdata, n=stop, parallel=true)
CSV.write(datadir("sims", "survival.csv"), mdf)
