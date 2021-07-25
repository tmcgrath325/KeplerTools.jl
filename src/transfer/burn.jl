struct Burn
    orbit::Orbit
    Δv̄::SVector{3,Float64}
    time::Float64
end

orientation_components(brn::Burn) = time_orientation_rotation(brn.time, brn.orbit) \ brn.Δv̄


### descriptive display ###

function Base.show(io::IO, brn::Burn)
    o_Δv̄ = orientation_components(brn)
    println(io,
        "Instantaneous burn in orbit above $(brn.orbit.primary.name):\n",
        "  UT:\t\t$(brn.time) s\n",
        "  total Δv:\t$(norm(brn.Δv̄)) m/s\n",
        "  prograde:\t$(o_Δv̄[1]) m/s\n",
        "  normal:\t$(o_Δv̄[3]) m/s\n",
        "  radial:\t$(o_Δv̄[2]) m/s\n",
    )
    return
end