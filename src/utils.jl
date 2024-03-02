export Datadir
# need custom datadir for cluster
# but cannot overwrite the one from DrWatson
function Datadir(args...)
    if haskey(ENV, "SCRATCH")
        return joinpath(ENV["SCRATCH"], projectname(), args...)
    else
        return datadir(args...)
    end
end
