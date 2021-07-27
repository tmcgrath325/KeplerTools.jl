### type definitions ###

"""
A space object considered to have negligible mass compared to the body around which it orbits
(e.g. a spacecraft or an asteroid).
"""
struct SpaceObject{P}
    name::AbstractString
    orbit::Orbit{P}
    primary::P
end

"""
A celestial body (e.g. a planet or moon) that is in orbit around a larger body.
"""
struct CelestialBody{P} <: CelestialObject
    name::AbstractString
    eqradius::Float64
    μ::Float64
    SoI::Float64
    rotperiod::Float64
    rotinitial::Float64
    orbit::Orbit{P}
    primary::CelestialObject
    color::Tuple{Int, Int, Int}
    satellite_bodies::Vector{<:CelestialBody}
    satellite_objects::Vector{<:SpaceObject}
end

"""
A star that is not in orbit around a larger body. Serves as the root of a tree defining a solar system.
"""
struct Star <: CelestialObject
    name::AbstractString
    eqradius::Float64
    μ::Float64
    color::Tuple{Int, Int, Int}
    satellite_bodies::Vector{CelestialBody{Star}}
    satellite_objects::Vector{SpaceObject{Star}}
end

### alternate constructors ###

# function SpaceObject(name::AbstractString, orbit::Orbit, primary::CelestialObject)
#     ob = SpaceObject(name::AbstractString, orbit::Orbit, primary::CelestialObject)
#     add_to_primary!(ob)
#     return bd
# end

"""
    soi = (μ, μprim, a)

Computes the sphere of influence of a body with semimajor axis `a` and gravity parameter `μ` in orbit around a larger body
with gravity parameter `μprim`.
"""
set_soi(μ, μprim, a) = a *(μ/μprim)^(2/5)

"""
Adds a `CelestialBody` or `SpaceObject` to the lists of satellites around its primary body.
"""
function add_to_primary!(bd::CelestialBody)
    push!(bd.primary.satellite_bodies, bd)
    sort!(bd.primary.satellite_bodies)
end

function add_to_primary!(ob::SpaceObject)
    push!(ob.primary.satellite_objects, bd)
    sort!(ob.primary.satellite_objects)
end

"""
Creates a `CelestialBody` and adds it to the lists of satellites around its primary body.
"""
function CelestialBody(name::AbstractString, eqradius, μ, SoI, rotperiod, rotinitial, orbit::Orbit, color::Tuple{Int,Int,Int}=(255,255,255))
    bd = CelestialBody(name, eqradius, μ, SoI, rotperiod, rotinitial, orbit, orbit.primary, color, CelestialBody[], SpaceObject[])
    add_to_primary!(bd)
    return bd
end

function CelestialBody(name::AbstractString, eqradius, μ, rotperiod, rotinitial, orbit::Orbit, color::Tuple{Int,Int,Int}=(255,255,255))
    bd = CelestialBody(name, eqradius, μ, set_soi(μ, orbit.primary.μ, orbit.a), rotperiod, rotinitial, orbit, orbit.primary, color, CelestialBody{<:CelestialBody}[], SpaceObject{<:CelestialBody}[])
    add_to_primary!(bd)
    return bd
end

Star(name::AbstractString, eqradius, μ, color::Tuple{Int,Int,Int}=(255,255,255)) = Star(name, eqradius, μ, color, CelestialBody{Star}[], SpaceObject{Star}[])


### helper methods ###

"""
    sats = collect_satellite_bodies(bd)

Collects all planets or moons that are contained within the sphere of influence of `bd`.
"""
function collect_satellite_bodies(bd::CelestialObject)
    if isempty(bd.satellite_bodies)
        return CelestialBody[]
    else
        gchildren = CelestialBody[]
        for st in bd.satellite_bodies
            push!(gchildren, collect_satellite_bodies(st)...)
        end
        return CelestialBody[bd.satellite_bodies..., gchildren...]
    end
end

"""
    sats = collect_satellite_bodies(bd)

Collects all space objects that are contained within the sphere of influence of `bd`.
"""
function collect_satellite_objects(bd::CelestialObject)
    children = collect_satellite_bodies(bd)
    sat_objects = copy(bd.satellite_objects)
    for child in children
        push!(sat_objects, child.satellite_objects...)
    end
    return sat_objects
end

"""
    parent = closest_common_parent(bd1, bd2)

Returns the `CelestialObject` corresponding to the smallest system to which both `bd1` and `bd2` belong.
"""
function closest_common_parent(bd1::Union{CelestialObject,SpaceObject}, bd2::Union{CelestialObject,SpaceObject})
    if bd1 == bd2 || bd2 ∈ collect_satellite_bodies(bd1)
        return bd1
    end
    bd = bd1
    while hasfield(typeof(bd), :orbit)
        bd = bd.orbit.primary
        if bd2 ∈ collect_satellite_bodies(bd) || bd==bd2
            return bd
        end
    end
    return nothing
end

closest_common_parent(orb1::Orbit, orb2::Orbit) = closest_common_parent(orb1.primary, orb2.primary)

"""
    path = path_to_parent(bd, parent)

Returns a vector of `CelestialObject`s traversed in the solar system tree to arrive at `parent`
starting from `bd`. `parent` is assumed to be at a higher level than `bd`
"""
function path_to_parent(bd::Union{CelestialObject,SpaceObject}, parent::CelestialObject)
    if bd == parent return [bd] end
    bodies = Union{CelestialObject,SpaceObject}[bd]
    pbd = bd
    while hasfield(typeof(pbd), :orbit)
        pbd = pbd.orbit.primary
        push!(bodies, pbd)
        if pbd == parent
            return bodies
        end
    end
    return
end

"""
    path = path_to_parent(startbd, endbd)

Returns a vector of `CelestialObject`s traversed in the solar system tree to arrive at `endbd`
starting from `startbd`.
"""
function path_to_body(startbd::Union{CelestialObject,SpaceObject}, endbd::Union{CelestialObject,SpaceObject})
    commonparent = closest_common_parent(startbd, endbd)
    pathup = path_to_parent(startbd, commonparent)
    pathdown = reverse(path_to_parent(endbd, commonparent))
    return [pathup..., pathdown[2:end]...]
end


### comparison and sorting methods ###

function Base.isless(bd1::CelestialObject, bd2::CelestialObject)
    if bd1.orbit.a == bd2.orbit.a
        return isless(bd1.name, bd2.name)
    else
        return isless(bd1.orbit.a, bd2.orbit.a)
    end
end


### descriptive display ###

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
