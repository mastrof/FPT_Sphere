using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :FPT_Sphere
@everywhere begin
     using CSV, DataFrames
     using StatsBase
     using LsqFit
end


function fit_expdecay(dataset::AbstractDataFrame)
    grouped_dataset = groupby(dataset, [:U, :R, :C1])
    df = @distributed append! for g in grouped_dataset
        replicates = unique(g.rep)
        τs = zeros(length(replicates))
        δτs = zeros(length(replicates))
        U = g.U[1]
        R = g.R[1]
        C1 = g.C1[1]
        Δt = g.Δt[1]
        for rep in replicates
            data = subset(g, :rep => r -> r .== rep; view=true)
            n = data.nagents ./ data.nagents[1]
            t = data.time .* Δt
            τs[rep], δts[rep] = fit_expdecat(n, t)
        end
        n = g.nagents ./ g.nagents[1]
        t = g.time .* Δt
        τ = mean(τs)
        δτ = mean(δτs)
        DataFrame(; U, R, C1, τ, δτ)
    end
    return df
end
@everywhere function fit_expdecay(S::AbstractVector, t::AbstractVector)
    fit = curve_fit(expdecay_log, t, log.(S), [1.0])
    τ = fit.param[1]
    δτ = stderror(fit)[1]
    return τ, δτ
end
@everywhere expdecay_log(t, p) = -t / p[1]

datafiles = filter(
    f -> contains(f, r"survival.*csv"),
    readdir(Datadir("sims"); join=true)
)
raw_df = @distributed append! for filename in datafiles
    CSV.read(filename, DataFrame)
end
df = fit_expdecay(raw_df)
CSV.write(Datadir("proc", "survival.csv"), DataFrame)
