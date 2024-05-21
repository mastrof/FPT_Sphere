using Distributed
@everywhere using DrWatson
@everywhere @quickactivate :FPT_Sphere
@everywhere begin
     using CSV, DataFrames
     using LsqFit
end


function fit_expdecay(dataset::AbstractDataFrame)
    grouped_dataset = groupby(dataset, [:U, :R, :C1])
    df = @distributed append! for g in grouped_dataset
        U = g.U[1]
        R = g.R[1]
        C1 = g.C1[1]
        Δt = g.Δt[1]
        n = g.nagents ./ g.nagents[1]
        t = g.time .* Δt
        τ, δτ = fit_expdecay(n, t)
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

dataset = CSV.read(Datadir("sims", "survival.csv"), DataFrame)
df = fit_expdecay(dataset)
CSV.write(Datadir("proc", "survival.csv"), df)
