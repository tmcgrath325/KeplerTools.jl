function draw_transfer(tfer::Transfer, t=tfer.departure_time; bodies=true, sois=true, objects=false, refline=false, arrow=false, angle=false)
    min_time = tfer.departure_time
    max_time = tfer.arrival_time
    traces = AbstractTrace[]
    push!(traces, draw_system(tfer.transfer_orbits[1].primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline)...)
    # starting orbit (or orbit of starting body)
    if tfer.startorbit.primary == tfer.transfer_orbits[1].primary
        push!(traces, draw_orbit(tfer.startorbit, t))
        r̄ₛ = time_orbital_position(min_time, tfer.startorbit)
        sorb = tfer.startorbit
    else
        sorb = tfer.startorbit.primary.orbit
        while sorb.primary != tfer.transfer_orbits[1].primary
            sorb = eorb.primary.orbit
        end
        r̄ₛ = time_orbital_position(min_time, sorb)
    end
    # ending orbit (or orbit of ending body)
    if tfer.endorbit.primary == tfer.transfer_orbits[1].primary
        push!(traces, draw_orbit(tfer.endorbit, t))
        r̄ₑ = time_orbital_position(min_time, tfer.endorbit)
    else
        eorb = tfer.endorbit.primary.orbit
        while eorb.primary != tfer.transfer_orbits[1].primary
            eorb = eorb.primary.orbit
        end
        r̄ₑ = time_orbital_position(min_time, eorb)
    end
    # transfer orbits
    for (i,torb) in enumerate(tfer.transfer_orbits)
        if i==length(tfer.transfer_orbits)
            etime = max_time
        else
            etime = tfer.transfer_orbits[i+1].epoch
        end
        push!(traces, draw_orbit(torb, torb.epoch, etime))
        if torb.epoch <= t <= etime
            push!(traces, draw_orbit_position(torb, t))
        end
    end

    if arrow
        for brn in tfer.burns
            if brn.orbit ∈ tfer.transfer_orbits
                push!(traces, KeplerTools.draw_burn_arrow(brn; normalize=1)...)
                break
            end
        end
    end

    if angle
        push!(traces, KeplerTools.draw_angle_in_plane(r̄ₑ, r̄ₛ, sorb.basis)...)
        push!(traces, KeplerTools.draw_angle_in_plane(time_orbital_position(tfer.departure_time, tfer.transfer_orbits[1]), time_orbital_position(tfer.arrival_time, tfer.transfer_orbits[end]), tfer.transfer_orbits[1].basis)...)
    end
    return traces
end

function draw_ejection(eorb::Orbit, pkorb::Orbit, t=eorb.epoch, burns=Burn[]; bodies=true, sois=true, objects=false, refline=false, arrow=false, angle=false)
    min_time = eorb.epoch
    max_time = ejection_time(eorb)
    traces = AbstractTrace[]
    push!(traces, draw_system(eorb.primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline)...)
    push!(traces, draw_orbit(pkorb, min_time))
    push!(traces, draw_orbit(eorb, min_time, max_time))
    push!(traces, KeplerTools.draw_orbit_position(eorb, t))
    if arrow
        for brn in burns
            if brn.orbit == pkorb
                push!(traces, KeplerTools.draw_burn_arrow(brn)...)
                break
            end
        end
    end
    if angle
        v̂ₚ = normalize(time_orbital_velocity(min_time, eorb.primary.orbit))
        r̄ₒ = time_orbital_position(min_time, eorb)
        push!(traces, KeplerTools.draw_angle_in_plane(r̄ₒ, v̂ₚ, pkorb.basis)...)
    end
    return traces
end

function draw_insertion(aorb::Orbit, pkorb::Orbit, t=aorb.epoch, burns=Burn[]; bodies=true, sois=true, objects=false, refline=false, arrow=false,) # angle=false)
    min_time = insertion_time(aorb)
    max_time = aorb.epoch
    traces = AbstractTrace[]
    push!(traces, draw_system(aorb.primary, t; bodies=bodies, sois=sois, objects=objects, refline=refline)...)
    push!(traces, draw_orbit(pkorb, max_time))
    push!(traces, draw_orbit(aorb, min_time, max_time))
    push!(traces, KeplerTools.draw_orbit_position(aorb, t))
    if arrow
        for brn in burns
            if brn.orbit == pkorb
                push!(traces, KeplerTools.draw_burn_arrow(brn)...)
                break
            end
        end
    end
    # if angle
    #     v̂ₚ = normalize(time_orbital_velocity(min_time, aorb.primary.orbit))
    #     r̄ₒ = time_orbital_position(min_time, aorb)
    #     push!(traces, KeplerTools.draw_angle_in_plane(r̄ₒ, v̂ₚ, pkorb.basis)...)
    # end
    return traces
end