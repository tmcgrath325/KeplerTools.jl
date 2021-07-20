# type definitions

abstract type CelestialObject end

struct Orbit
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    Mo::Float64
    epoch::Float64
    primary::CelestialObject
    attrs::Dict{Symbol, Any}
end

struct StateVector
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
end

struct OrbitalState
    time::Float64
    statevector::StateVector
    orbit::Orbit
end

# alternate constructors

flight_path_angle(θ, orb::Orbit) = atan(orb.e*sin(θ)/(1+orb.e*cos(θ)))

function state_vector(t, orb::Orbit)
    θ = time_to_true(t, orb)
    r = orb.a*(1-orb.e^2)/(1+orb.e*cos(θ))              # height
    ō = r * @SVector([cos(θ),sin(θ),0])                 # position in orbital plane
    ϕ = flight_path_angle(θ, orb)                       # flight path angle
    v = √(orb.primary.μ*(2/r - 1/orb.a))               # speed
    dōdt = v * @SVector([cos(θ+π/2-ϕ),sin(θ+π/2-ϕ),0])  # velocity in orbital plane
    # rotate to reference frame
    r̄ = perif_to_inert_bases(ō, orb)
    v̄ = perif_to_inert_bases(dōdt, orb)
    return r̄, v̄
end

Orbit(a,body::CelestialObject) = Orbit(a,0,0,0,0,0,0,body,Dict{Symbol,Any}())
Orbit(a,e,i,Ω,ω,Mo,body::CelestialObject) = Orbit(a,e,i,Ω,ω,Mo,0,body,Dict{Symbol,Any}())
Orbit(a,e,i,Ω,ω,Mo,epoch,body::CelestialObject) = Orbit(a,e,i,Ω,ω,Mo,epoch,body,Dict{Symbol,Any}())

StateVector(t, orb::Orbit) = StateVector(state_vector(t,orb)...)

function Orbit(t::Real, r̄::SVector{3,<:Real}, v̄::SVector{3,<:Real}, prim::CelestialObject, epoch::Real=t)
    r = norm(r̄) 
    v = norm(v̄)
    μ = prim.μ
    if r == 0
        throw(ArgumentError("Invalid position"))
    end
    a = 1 / (2/r - (v^2)/μ)
    h̄ = cross(r̄,v̄)
    i = acos(h̄[3]/norm(h̄))
    ē = cross(v̄,h̄)/μ - r̄/r
    e = norm(ē)
    if isapprox(e, 0, atol=1e-15)       # if circular
        e = 0.
        ê = @SVector([1.,0.,0.])        # set periapsis at reference x
    else
        ê = ē/e
    end
    n̄ = cross(@SVector([0.,0.,1.]), h̄)
    n = norm(n̄)
    if isapprox(n, 0, atol=1e-15)       # if uninclined
        n = 0.
        n̂ = ê                           # set Ω at periapsis
    else
        n̂ = n̄/n
    end
    Ω = try bound_angle(copysign(1,n̂[2]) * acos(n̂[1]))
    catch error
        if isa(error, DomainError)
            Ω = bound_angle(copysign(1,n̂[2]) * acos(copysign(1,n̂[1])))
        end
    end
    ω = try bound_angle(copysign(1,ē[3]) * acos(dot(n̂,ê)))
    catch error
        if isa(error, DomainError)
            Ω = bound_angle(copysign(1,ē[3]) * acos(copysign(1,dot(n̂,ê))))
        end
    end
    θ = angle_in_plane(r̄, basis_matrix(Ω, ω, i))
    M = true_to_mean(θ, e)
    Δt = t - epoch 
    Mo = bound_angle(M - time_to_mean(Δt, period(a,μ)))
    return Orbit(a,e,i,Ω,ω,Mo,epoch,prim)
end

Orbit(t::Float64, statevec::StateVector, prim::CelestialObject, epoch::Real=t)= Orbit(t, statevec.position, statevec.velocity, prim, epoch)

OrbitalState(t::Real, orb::Orbit) = OrbitalState(t, StateVector(t,orb), orb)
OrbitalState(t::Real, statevec::StateVector, prim::CelestialObject, epoch::Real=t) = OrbitalState(t, StateVector, Orbit(t, statevec, prim, epoch))

# comparison and sorting methods

function Base.isapprox(orb1::Orbit, orb2::Orbit; atol=1e-6)
    a = isapprox(orb1.a, orb2.a, atol=orb2.a*atol)
    e = isapprox(orb1.e, orb2.e, atol=orb2.e*atol)
    i = isapproxangle(orb1.i, orb2.i, atol=atol)
    θ = isapproxangle(angle_in_plane(orb1, orb2, 0), 0, atol=atol)
    # Ω = isapproxangle(orb1.Ω, orb2.Ω, atol=atol)
    # ω = isapproxangle(orb1.ω, orb2.ω, atol=atol)
    return a && e && i && θ # && Ω && ω
end

# descriptive display


Base.show(io::IO, orb::Orbit) = println(io,
    "Orbit above $(orb.primary.name):\n",
    "  semi-major axis:\t\t$(orb.a) m\n",
    "  eccentricity:\t\t\t$(orb.e)\n",
    "  inclination:\t\t\t$(deg2rad(orb.i))°\n",
    "  RA of ascending node:\t\t$(deg2rad(orb.Ω))°\n",
    "  Argument of periapsis:\t$(deg2rad(orb.ω))°\n",
    "  epoch:\t\t\t$(orb.epoch) s"
)
