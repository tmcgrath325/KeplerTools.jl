function draw_transfer(tfer::Transfer, t=tfer.departure_time; bodies=true, sois=true, arrows=true, angles=true, objects=false, refline=false, kwargs...)
    min_time = tfer.departure_time
    max_time = tfer.arrival_time
    traces = AbstractTrace[]
    append!(traces, draw_system(tfer.transfer_orbits[1].primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline, kwargs...))

    # starting orbit (or orbit of starting body)
    if tfer.startorbit.primary == tfer.transfer_orbits[1].primary
        append!(traces, draw_orbit(tfer.startorbit, t; kwargs...))
    end
    # ending orbit (or orbit of ending body)
    if tfer.endorbit.primary == tfer.transfer_orbits[1].primary
        append!(traces, draw_orbit(tfer.endorbit, t; kwargs...))
    end
    # transfer orbits
    for (i,torb) in enumerate(tfer.transfer_orbits)
        if i==length(tfer.transfer_orbits)
            etime = max_time
        else
            etime = tfer.transfer_orbits[i+1].epoch
        end
        append!(traces, draw_orbit(torb, torb.epoch, etime; kwargs...))
        if torb.epoch <= t <= etime
            push!(traces, draw_orbit_position(torb, t))
        end
    end

    if arrows
        systemburns = filter(x->x.orbit.primary == tfer.transfer_orbits[1].primary, tfer.burns)
        if !isempty(systemburns)
            biggestburn = maximum([norm(b.Δv̄) for b in systemburns])
            for brn in systemburns
                append!(traces, KeplerTools.draw_burn_arrow(brn; normalize=biggestburn))
            end
        end
    end

    if angles
        sorb = isempty(tfer.ejection_orbits)  ? tfer.startorbit : tfer.ejection_orbits[end].primary.orbit
        eorb = isempty(tfer.insertion_orbits) ? tfer.endorbit   : tfer.insertion_orbits[1].primary.orbit
        r̄ₛ = time_orbital_position(tfer.departure_time, sorb)
        r̄ₑ = time_orbital_position(min_time, eorb)
        append!(traces, KeplerTools.draw_angle_in_plane(r̄ₑ, r̄ₛ, sorb.basis;))

        r̄ₛ = time_orbital_position(tfer.arrival_time, tfer.transfer_orbits[end])
        r̄ₑ = time_orbital_position(tfer.departure_time, tfer.transfer_orbits[1])
        append!(traces, KeplerTools.draw_angle_in_plane(r̄ₑ, r̄ₛ, basis_rotation(r̄ₛ, r̄ₑ); wrap_min=0.))
    end
    return traces
end

function draw_ejection(dorb::Orbit, pkorb::Orbit, t=dorb.epoch, burns=Burn[]; drawpark=true, bodies=true, sois=true, objects=false, refline=false, arrows=false, angles=false, kwargs...)
    min_time = dorb.epoch
    max_time = ejection_time(dorb)
    traces = AbstractTrace[]
    append!(traces, draw_system(dorb.primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline, kwargs...))
    if drawpark
        append!(traces, draw_orbit(pkorb, min_time; kwargs...))
    end
    append!(traces, draw_orbit(dorb, min_time, max_time; kwargs...))
    if dorb.epoch <= t <= ejection_time(dorb) 
        push!(traces, KeplerTools.draw_orbit_position(dorb, t))
    end
    if arrows
        systemburns = filter(x->x.orbit.primary == pkorb.primary, burns)
        if !isempty(systemburns)
            biggestburn = maximum([norm(b.Δv̄) for b in systemburns])
            for brn in systemburns
                append!(traces, KeplerTools.draw_burn_arrow(brn; normalize=biggestburn))
            end
        end
    end
    if angles
        v̂ₚ = normalize(time_orbital_velocity(min_time, dorb.primary.orbit))
        r̄ₒ = time_orbital_position(min_time, dorb)
        append!(traces, KeplerTools.draw_angle_in_plane(r̄ₒ, v̂ₚ, pkorb.basis))
    end
    return traces
end

function draw_insertion(aorb::Orbit, pkorb::Orbit, t=aorb.epoch, burns=Burn[]; drawpark=true, bodies=true, sois=true, objects=false, refline=false, arrows=false, angles=false, kwargs...) # angle=false)
    min_time = insertion_time(aorb)
    max_time = aorb.epoch
    traces = AbstractTrace[]
    append!(traces, draw_system(aorb.primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline, kwargs...))
    if drawpark
        append!(traces, draw_orbit(pkorb, max_time; kwargs...))
    end
    append!(traces, draw_orbit(aorb, min_time, max_time; kwargs...))
    if insertion_time(aorb) <= t <= aorb.epoch
        push!(traces, KeplerTools.draw_orbit_position(aorb, t))
    end
    if arrows
        systemburns = filter(x->x.orbit.primary == pkorb.primary, burns)
        if !isempty(systemburns)
            biggestburn = maximum([norm(b.Δv̄) for b in systemburns])
            for brn in systemburns
                append!(traces, KeplerTools.draw_burn_arrow(brn; normalize=biggestburn))
            end
        end
    end
    # if angles
    #     v̂ₚ = normalize(time_orbital_velocity(min_time, aorb.primary.orbit))
    #     r̄ₒ = time_orbital_position(min_time, aorb)
    #     append!(traces, KeplerTools.draw_angle_in_plane(r̄ₒ, v̂ₚ, pkorb.basis))
    # end
    return traces
end