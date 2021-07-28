struct Burn
    orbit::Orbit
    Δv̄::SVector{3,Float64}
    time::Float64
end

function Burn(time, prograde, normal, radial, orb::Orbit)
    rot = orientation_rotation(time_to_true(time, orb), orb)
    Δv̄ = rot*@SVector([prograde, normal, -radial])
    return Burn(orb, Δv̄, time)
end

# methods for instantaneous orientation (prograde, radial, normal) 

"""
    rotmat = orientation_rotation(pos, vel)
    rotmat = orientation_rotation(θ, orb)

Computes the rotation matrix which defines, the prograde, radial, and normal directions
for an orbital state vector `pos`, `vel`.

If a true anomaly `θ` and Orbit `orb` are provided, they are used to compute the state vector,
which will be used to compute the rotation matrix.
"""
orientation_rotation(pos, vel, h=cross(pos, vel)
    ) = basis_rotation(normalize(vel), normalize(h))

orientation_rotation(θ, orb::Orbit) = orientation_rotation(state_vector(θ, orb)...)

"""
    prograde, normal, radial = orientation_components(brn)

Computes the prograde, normal, and radial components of the Burn `brn`.
"""
function orientation_components(brn::Burn) 
    vec = time_orientation_rotation(brn.time, brn.orbit) \ brn.Δv̄
    return vec[1], vec[2], -vec[3]
end


apply_burn(brn::Burn
    ) = Orbit(brn.time, time_orbital_position(brn.time, brn.orbit), time_orbital_velocity(brn.time, brn.orbit)+brn.Δv̄, brn.orbit.primary)


### descriptive display ###

function Base.show(io::IO, brn::Burn)
    o_Δv̄ = orientation_components(brn)
    println(io,
        "Instantaneous burn in orbit above $(brn.orbit.primary.name):\n",
        "  UT:\t\t$(brn.time) s\n",
        "  total Δv:\t$(norm(brn.Δv̄)) m/s\n",
        "  prograde:\t$(o_Δv̄[1]) m/s\n",
        "  normal:\t$(o_Δv̄[2]) m/s\n",
        "  radial:\t$(-o_Δv̄[3]) m/s\n",
    )
    return
end