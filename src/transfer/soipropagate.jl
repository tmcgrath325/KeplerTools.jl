function soi_distance(t, orb::Orbit, bd::CelestialBody)
    if orb.primary != bd.orbit.primary
        throw(ArgumentError("orbit and body do not share a primary body"))
    end
    return norm(time_orbital_position(t, orb) - time_orbital_position(t, bd.orbit)) - bd.SoI
end

"""
    nextorb = propagate_across_soi(t, orb)

Searches for SoI changes starting at time `t`, and returns the orbit after the SoI
patch if one occurs.

For periodic (elliptical) orbits, the search begins at time `t` and ends at a full period
beyond time `t`. If no SoI change occurs, the original orbit is returned.

For hyperbolic orbits, the search begins at time `t` and ends at the time of SoI ejection.
If the orbit's primary body has no SoI (i.e. it is a Star), the search ends at the time when
the orbit is gauranteed to have passed beyond the outermost satellite of the star.
"""
function propagate_across_soi(t, orb::Orbit)
    if orb.e < 1
        endtime = orb.period
    else
        if typeof(orb.primary) == Star
            f_sat = max(orb.primary.satellite_bodies)
            endtime = true_to_time(orbital_angle(f_sat.a*(1+f_sat.e)+f_sat.SoI, orb), orb)
        else
            endtime = ejection_time(orb)
        end
    end

    intersects = Tuple{CelestialBody, Float64}[]
    for bd in orb.primary.satellite_bodies
        #TO DO: ensure that earliest intersection is found, if multiple exist
        obj_fun(t) = soi_distance(t, orb, bd)
        res = optimize(obj_fun, t, endtime, Brent())
        if Optim.minimum(res) < 0
            push!(intersects, (bd, Optim.minimizer(res)))
            endtime = Optim.minimizer(res)
        end
    end

    if !isempty(intersects)
        intsct_bd, mindist_time = intersects[end]
        obj_fun2(t) = abs(soi_distance(t, orb, intsct_bd))
        res2 = optimize(obj_fun2, t, mindist_time, Brent())
        intsct_time = Optim.minimizer(res2)
        r̄ = time_orbital_position(intsct_time, orb) - time_orbital_position(intsct_time, intsct_bd.orbit)
        v̄ = time_orbital_velocity(intsct_time, orb) - time_orbital_velocity(intsct_time, intsct_bd.orbit)
        @show r̄, v̄, intsct_time 
    elseif orb.e > 1
        intsct_bd = orb.primary.orbit.primary
        intsct_time = endtime
        r̄ = time_orbital_position(intsct_time, orb) + time_orbital_position(intsct_time, orb.primary.orbit)
        v̄ = time_orbital_velocity(intsct_time, orb) + time_orbital_velocity(intsct_time, orb.primary.orbit)
    else
        return orb
    end
    return Orbit(intsct_time, r̄, v̄, intsct_bd)
end