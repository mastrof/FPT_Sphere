export encounter!, reinsert_agent!

function encounter!(model; callback=remove_agent!)
    for microbe in allagents(model)
        encounter!(microbe, model; callback)
        save_position(microbe, model)
    end
end

function encounter!(microbe, model; callback=remove_agent!)
    if is_encounter(microbe, model)
        callback(microbe, model)
    end
end

function is_encounter(microbe, model)
    sphere = source(model)
    P = source_origin(model)
    R = source_size(model)
    a = radius(microbe)
    S = R + a
    pos = position(microbe)
    d = euclidean_distance(pos, P, model)
    U = speed(microbe)
    dt = model.timestep
    if d <= S
        return true
    elseif d <= S + U*dt
        from = lastposition(model)[microbe.id]
        to = distancevector(from, pos, model)
        step_segment = Meshes.Ray(Meshes.Point(from), Meshes.Vec(to))
        sphere = Meshes.Ball(Meshes.Point(source_origin(model)), S)
        return Meshes.intersects(step_segment, sphere)
    else
        return false
    end
end

function reinsert_agent!(microbe, model)
    R = source_size(model)
    P = source_origin(model)
    a = radius(microbe)
    newpos = random_position(model)
    # new position must be outside the source volume
    while euclidean_distance(newpos, P, model) <= R + a
        newpos = random_position(model)
    end
    move_agent!(microbe, newpos, model)
    microbe.speed = random_speed(microbe, model)
    microbe.vel = random_velocity(model)
end


function save_position(microbe, model)
    lastposition(model)[microbe.id] = position(microbe)
end
lastposition(model) = model.last_position
