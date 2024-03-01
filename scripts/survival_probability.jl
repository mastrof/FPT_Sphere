using Distributed
@everywhere using DrWatson
@everywhere begin
    @quickactivate :FPT_Sphere
    using Agents, MicrobeAgents
    using DataFrames
end
using CSV
using LsqFit
using GLMakie

##
@everywhere function initialize(;
    N = 2500,
    L = 1000.0,
    R = 50.0,
    C₀ = 0.0,
    C₁ = 1.0,
    U = 20.0,
    Δt = 0.1,
)
    D = 3
    periodic = true
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
        :Δt => Δt,
        :L => 1000.0,
        :R => 100.0,
        #:C₁ => [0.0, 0.5, 1.0, 2.0, 5.0, 10.0],
        :C₁ => [0.0, 1.0, 2.0],
        :U => [15.0, 30.0, 60.0]
    )
end
_, mdf = paramscan(parameters, initialize; mdata, n=stop, parallel=true)
CSV.write(datadir("sims"), "survival.csv")

##
fig = Figure(; size=(800,600), fontsize=32)
ax = Axis(fig[1,1];
    xgridvisible=false, ygridvisible=false,
    yscale=log10,
    xlabel = "time (s)",
    ylabel = "survival probability"
)
gdf_U = groupby(mdf, :U)
linestyles = [:solid, :dash, :dot] # changes with U
linecolors = cgrad(:roma, length(parameters[:C₁]); categorical=true) # changes with C
for (i,h) in enumerate(gdf_U)
    gdf_C = groupby(h, :C₁)
    for (j,g) in enumerate(gdf_C)
        t = g.step .* Δt
        φ = g.nagents ./ g.nagents[1]
        C = g.C₁[1]
        U = g.U[1]
        lab = "U = $U μm/s; C₁ = $(C) μM"
        lines!(ax, t, φ;
            linewidth=5, label=lab,
            linestyle=linestyles[i],
            color=linecolors[j],
        )
    end
end
fig

##
fig = Figure(; size=(800,600), fontsize=24)
makeax(i,j) = Axis(fig[i,j];
        xgridvisible=false, ygridvisible=false,
        yscale=log10,
        xlabel="time (s)", ylabel="survival probability",
    )
gdf_U = groupby(mdf, :U)
linecolors = cgrad(:viridis, length(parameters[:C₁]); categorical=true)
for (i,h) in enumerate(gdf_U)
    ax = makeax(1,i)
    U = h.U[1]
    ax.title = "U = $U μm/s"
    for (j,g) in enumerate(groupby(h, :C₁))
        t = g.step .* Δt
        φ = g.nagents ./ g.nagents[1]
        C = g.C₁[1]
        lines!(ax, t, φ; linewidth=1, color=linecolors[j])
        # fit
        _model(x,p) = -x/p[1]
        τ = curve_fit(_model, t, log10.(φ), [1.0]).param[1]
        lines!(ax, t[1:10:end], t -> exp10(-t/τ);
            linewidth=3, linestyle=:dash, color=linecolors[j],
            label="C₁ = $C μM; τ = $(round(Int, τ/60)) min"
        )
    end
    ylims!(ax, (1e-3, 1.01))
    axislegend(ax; position=:lb, patchsize=(30,20))
end
fig
