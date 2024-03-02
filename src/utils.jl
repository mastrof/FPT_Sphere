# extension to use SCRATCH folder as datadir on cluster
if haskey(ENV, "SCRATCH")
    function DrWatson.datadir(args...)
        return joinpath(ENV["SCRATCH"], projectname(), args...)
    end
end
