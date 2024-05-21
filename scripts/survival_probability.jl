using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :FPT_Sphere
@everywhere begin
    using Agents, MicrobeAgents
    using DataFrames
end
using CSV
using LsqFit

##
@everywhere function initialize(;
    N = 2500,
    # Φ = 1 - 1e-4, # free volume fraction -> (1-Φ)L^3 = V // for later
    L = 1000.0,
    R = 50.0,
    C0 = 0.0,
    C1 = 1.0,
    U = 20.0,
    Δt = 0.05,
    kwargs...
)
    D = 3
    periodic = true
    # define L by fixing excluded volume // for later
    # V = 4π/3 * R^3 # μm^3
    # L = (V/(1-Φ))^(1/3) # μm
    space = ContinuousSpace(ntuple(_ -> L, D); periodic)
    properties = Dict(
        :chemoattractant => DiffusiveField(C0, C1, R, spacesize(space) ./ 2),
        :last_position => [zero(SVector{D,Float64}) for _ in 1:N]
    )
    model = StandardABM(SonMenolascina{D,4}, space, Δt;
        properties, container=Dict,
        model_step! = encounter!,
    )
    motility = RunReverseFlick(1.0, [U], 1.0, [U])
    for i in 1:N
        add_agent!(model; motility)
        # after adding agent to model, save its assigned position
        model.last_position[i] = position(model[i])
    end
    model
end

@everywhere when_model(model, s) = (s % 1200 == 0)
@everywhere stop(model, s) = (s >= nsteps || nagents(model) <= 1)

@everywhere function production_run(config)
    model = initialize(; config...)
    dt = model.timestep
    T = 24*60*60 # s
    nsteps = round(Int, T/dt)
    mdata = [nagents]
    _, mdf = run!(model, stop; mdata, when_model)
    fout = Datadir("sims", savename("survival", config, "csv"))
    CSV.write(fout, mdf)
end

_, mdf = paramscan(parameters, initialize; mdata, n=stop, parallel=true)
CSV.write(Datadir("sims", "survival.csv"), mdf)

L = [500.0]
R = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
C1 = [0.0, 0.1, 0.25, 0.5, 1.0]
U = [15.0, 25.0, 40.0, 60.0]
num_replicates = 10
rep = collect(1:num_replicates)

allparams = @dict L R C1 U rep
dicts = dict_list(allparams)

pmap(production_run, dicts)
