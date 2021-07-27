function draw_system(bd::CelestialObject, starttime, endtime=nothing; bodies=true, sois=true, objects=false)
    traces = AbstractTrace[]
    push!(traces, draw_central_body(bd)...)
    push!(traces, draw_central_body_orbit(bd, starttime))
    if bodies
        for body in bd.satellite_bodies
            push!(traces, draw_orbiting_body(body, starttime)...)
            push!(traces, draw_orbit(body, starttime, endtime))
        end
    end
    if sois
        for body in bd.satellite_bodies
            push!(traces, draw_soi(body, starttime)...)
        end
    end
    if objects
        for ob in bd.satellite_objects
            push!(traces, draw_orbit(ob, starttime, endtime))
        end
    end
    return traces
end