function angle_in_plane(vec::SVector{3}, mat::Union{MRP,SMatrix{3,3}})
    planevec = mat*vec
    return atan(planevec[2],planevec[1])
end

function angle_in_plane(vec::SVector{3}, orb::Orbit)
    planevec = inert_to_perif_bases(vec, orb)
    return atan(planevec[2],planevec[1])
end

angle_in_plane(vec::SVector{3}, orb::Orbit, t) = bound_angle(angle_in_plane(vec, orb) - angle_in_plane(state_vector(t, orb)[1], orb))
angle_in_plane(vec::SVector{3}, vecref::SVector{3}, orb::Orbit) = bound_angle(angle_in_plane(vec, orb) - angle_in_plane(vecref, orb))
angle_in_plane(orb::Orbit, orbref::Orbit, t) = bound_angle(angle_in_plane(state_vector(t, orb)[1], orb) - angle_in_plane(state_vector(t, orbref)[1], orbref))
