function angle_in_plane(vec::AbstractVector{<:Real}, mat::Union{MRP,AbstractMatrix}=MRP(1I))
    planevec = mat*vec
    return atan(planevec[2],planevec[1])
end

angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, mat::Union{MRP,AbstractMatrix}=MRP(1I)
    ) = bound_angle(angle_in_plane(vec,mat) - angle_in_plane(vecref,mat))
angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, xdir::AbstractVector{<:Real}, ydir::AbstractVector{<:Real}
    ) = angle_in_plane(vec, vecref, basis_MRP(xdir,ydir))
angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}
    ) = angle_in_plane(vec, vecref, basis_MRP(vec,vecref))

# angle between the vector `vec` projected to the orbital plane and the periapsis of `orb`
angle_in_plane(vec::AbstractVector{<:Real}, orb::Orbit
    ) = angle_in_plane(vec, basis_MRP!(orb))

# angle between the vectors `vec` and `vecref` after projecting both to the orbital plane of `orb`
angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, orb::Orbit
    ) = angle_in_plane(vec, vecref, basis_MRP!(orb))

# angle between the vector `vec` projected to the orbital plane and the orbit's position at time `t`
angle_in_plane(vec::AbstractVector{<:Real}, orb::Orbit, t::Real
    ) = bound_angle(angle_in_plane(vec, state_vector(t, orb)[1], orb))

# angle between the orbital positions of `orb` from `orbref` at time `t`, after projecting the position `orb` to the orbital plane of `orbref`
angle_in_plane(orb::Orbit, orbref::Orbit, t::Real
    ) = bound_angle(angle_in_plane(state_vector(t, orb)[1], orbref, t))
