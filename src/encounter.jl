function encounter!(microbe, model; callback=remove_agent!)
    F = field(model)
    R = source_size(F)
    origin = origin(F)
    a = radius(microbe)
    if euclidean_distance(position(microbe), origin, model) <= R + a
        callback(microbe, model)
    end
end
