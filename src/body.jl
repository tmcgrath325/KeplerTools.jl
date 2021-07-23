# type definitions

struct SpaceObject
    name::AbstractString
    orbit::Orbit
    primary::CelestialObject
end

struct CelestialBody <: CelestialObject
    name::AbstractString
    eqradius::Float64
    μ::Float64
    SoI::Float64
    rotperiod::Float64
    rotinitial::Float64
    orbit::Orbit
    primary::CelestialObject
    satellite_bodies::Vector{CelestialBody}
    satellite_objects::Vector{SpaceObject}
end

struct Star <: CelestialObject
    name::AbstractString
    eqradius::Float64
    μ::Float64
    satellite_bodies::Vector{CelestialBody}
    satellite_objects::Vector{SpaceObject}
end

# alternate constructors

function SpaceObject(name::AbstractString, orbit::Orbit, primary::CelestialObject)
    ob = SpaceObject(name::AbstractString, orbit::Orbit, primary::CelestialObject)
    add_to_primary!(ob)
    return bd
end

set_soi(μ, μprim, a) = a *(μ/μprim)^(2/5)

function CelestialBody(name::AbstractString, eqradius, μ, SoI, rotperiod, rotinitial, orbit::Orbit)
    bd = CelestialBody(name, eqradius, μ, SoI, rotperiod, rotinitial, orbit, orbit.primary, CelestialBody[], SpaceObject[])
    add_to_primary!(bd)
    return bd
end

function CelestialBody(name::AbstractString, eqradius, μ, rotperiod, rotinitial, orbit::Orbit)
    bd = CelestialBody(name, eqradius, μ, set_soi(μ, orbit.primary.μ, orbit.a), rotperiod, rotinitial, orbit, orbit.primary, CelestialBody[], SpaceObject[])
    add_to_primary!(bd)
    return bd
end

Star(name::AbstractString, eqradius, μ) = Star(name, eqradius, μ, CelestialBody[], SpaceObject[])

# comparison and sorting methods

function Base.isless(orb1::Orbit, orb2::Orbit)
    fnames = fieldnames(Orbit)
    for fn in fnames
        if orb1.fn != orb2.fn
            return isless(orb1.fn, orb2.fn)
        end
    end
    return false
end

function Base.isless(bd1::CelestialObject, bd2::CelestialObject)
    if bd1.orbit.a == bd2.orbit.a
        return isless(bd1.name, bd2.name)
    else
        return isless(bd1.orbit.a, bd2.orbit.a)
    end
end

function add_to_primary!(bd::CelestialBody)
    push!(bd.primary.satellite_bodies, bd)
    sort!(bd.primary.satellite_bodies)
end

function add_to_primary!(ob::SpaceObject)
    push!(ob.primary.satellite_objects, bd)
    sort!(ob.primary.satellite_objects)
end

# descriptive display

Base.show(io::IO, st::Star) = println(io,
    """Star "$(st.name)":\n""",
    "  $(length(st.satellite_bodies)) orbiting celestial bodies\n",
    "  $(length(st.satellite_objects)) orbiting space objects"
)

Base.show(io::IO, bd::CelestialBody) = println(io,
    """CelestialBody "$(bd.name)" orbiting $(bd.primary.name):\n""",
    "  $(length(bd.satellite_bodies)) orbiting celestial bodies\n",
    "  $(length(bd.satellite_objects)) orbiting space objects"
)

Base.show(io::IO, ob::SpaceObject) = println(io,
    """SpaceObject "$(ob.name)" orbiting $(ob.primary.name):\n""",
    "  semi-major axis:\t\t$(ob.orb.a)\n",
    "  eccentricity:\t\t\t$(ob.orb.e)\n",
    "  inclination:\t\t\t$(rad2deg(ob.orb.i))°\n",
    "  RA of ascending node:\t\t$(rad2deg(ob.orb.Ω))°\n",
    "  Argument of periapsis:\t$(rad2deg(ob.orb.ω))°\n",
    "  mean anomaly at epoch:\t$(rad2deg(orb.Mo))°\n",
    "  epoch:\t\t\t$(ob.orb.epoch) s"
)
