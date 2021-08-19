function draw_system(bd::CelestialObject, starttime, endtime=nothing; bodies=true, sois=true, objects=false, refline=false, kwargs...)
    traces = AbstractTrace[]
    push!(traces, draw_central_body(bd))
    if hasfield(typeof(bd), :orbit)
        push!(traces, draw_central_body_orbit(bd, starttime))
    end
    if bodies
        for body in bd.satellite_bodies
            push!(traces, draw_orbiting_body(body, starttime))
            append!(traces, draw_orbit(body, starttime, endtime; kwargs...))
        end
    end
    if sois
        for body in bd.satellite_bodies
            push!(traces, draw_soi(body, starttime))
        end
    end
    if objects
        for ob in bd.satellite_objects
            append!(traces, draw_orbit(ob, starttime, endtime; kwargs...))
        end
    end
    if refline
        push!(traces, draw_ref_line(orbit_plot_size(bd)))
    end
    return traces
end