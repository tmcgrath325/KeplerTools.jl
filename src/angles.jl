"""
    Δθ = angle_in_plane(vec, mat)
    Δθ = angle_in_plane(vec1, vec2, mat)
    Δθ = angle_in_plane(vec, orbit)
    Δθ = angle_in_plane(vec, vec, orbit)
    Δθ = angle_in_plane(vec, orbit, t)
    Δθ = angle_in_plane(orbit1, orbit2, t)

Computes the angle between two vectors after projecting them to the orbital plane defined by an orbit's
basis matrix `mat`. 

If only one vector is provided, the angle between the vector `vec` and the orbit's periapsis
is found. 

If a vector, orbit, and time `t` are provided, the in-plane angle between the vector and the orbit's 
position at `t` is found.

If two orbits and a time `t` are provide, the in-plane angle between the two orbits' positions at time
`t` is found.
"""
function angle_in_plane(vec::AbstractVector{<:Real}, mat::Union{Rotation{3,Float64},AbstractMatrix}=MRP(1I))
    planevec = mat\vec
    return atan(planevec[2],planevec[1])
end

angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, mat::Union{Rotation{3,Float64},AbstractMatrix}=MRP(1I)
    ) = wrap_angle(angle_in_plane(vec,mat) - angle_in_plane(vecref,mat))
angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, xdir::AbstractVector{<:Real}, ydir::AbstractVector{<:Real}
    ) = angle_in_plane(vec, vecref, basis_rotation(xdir,ydir))

# angle between the vector `vec` projected to the orbital plane and the periapsis of `orb`
angle_in_plane(vec::AbstractVector{<:Real}, orb::Orbit
    ) = angle_in_plane(vec, orb.basis)

# angle between the vectors `vec` and `vecref` after projecting both to the orbital plane of `orb`
angle_in_plane(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real}, orb::Orbit
    ) = angle_in_plane(vec, vecref, orb.basis)

# angle between the vector `vec` projected to the orbital plane and the orbit's position at time `t`
angle_in_plane(vec::AbstractVector{<:Real}, orb::Orbit, t::Real
    ) = wrap_angle(angle_in_plane(vec, time_state_vector(t, orb)[1], orb))

# angle between the orbital positions of `orb` from `orbref` at time `t`, after projecting the position `orb` to the orbital plane of `orbref`
angle_in_plane(orb::Orbit, orbref::Orbit, t::Real
    ) = wrap_angle(angle_in_plane(time_state_vector(t, orb)[1], orbref, t))
