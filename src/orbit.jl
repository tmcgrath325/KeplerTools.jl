period(a::Real, μ::Real) = 2π*√(abs(a)^3 / μ)

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
    period::Float64
    basis::MRP{Float64}
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

flight_path_angle(θ::Real, orb::Orbit) = atan(orb.e*sin(θ)/(1+orb.e*cos(θ)))

function state_vector(t::Real, orb::Orbit)
    θ = time_to_true(t, orb)
    r = orb.a*(1-orb.e^2)/(1+orb.e*cos(θ))              # height
    ō = r * @SVector([cos(θ),sin(θ),0])                 # position in orbital plane
    ϕ = flight_path_angle(θ, orb)                       # flight path angle
    v = √(orb.primary.μ*(2/r - 1/orb.a))                # speed
    dōdt = v * @SVector([cos(θ+π/2-ϕ),sin(θ+π/2-ϕ),0])  # velocity in orbital plane
    # rotate to reference frame
    r̄ = perifocal_to_inertial_bases(ō, orb)
    v̄ = perifocal_to_inertial_bases(dōdt, orb)
    return r̄, v̄
end

Orbit(a,body::CelestialObject
    ) = Orbit(a,0.,0.,0.,0.,0.,0.,body,period(a,body.μ),basis_MRP(0.,0.,0.))
Orbit(a,e,i,Ω,ω,Mo,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,0,body,period(a,body.μ),basis_MRP(Ω, ω, i))
Orbit(a,e,i,Ω,ω,Mo,epoch,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,epoch,body,period(a,body.μ),basis_MRP(Ω, ω, i))

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
    θ = angle_in_plane(r̄, basis_MRP(Ω, ω, i))
    M = true_to_mean(θ, e)
    Δt = t - epoch 
    if e < 1
        Mo = bound_angle(M - time_to_mean(Δt, period(a,μ)))
    else
        Mo = M - time_to_mean(Δt, period(a,μ))
    end

    return Orbit(a,e,i,Ω,ω,Mo,epoch,prim)
end

Orbit(t::Real, statevec::StateVector, prim::CelestialObject, epoch::Real=t
    ) = Orbit(t, statevec.position, statevec.velocity, prim, epoch)

OrbitalState(t::Real, orb::Orbit
    ) = OrbitalState(t, StateVector(t,orb), orb)
OrbitalState(t::Real, statevec::StateVector, prim::CelestialObject, epoch::Real=t
    ) = OrbitalState(t, StateVector, Orbit(t, statevec, prim, epoch))

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

# TO DO: account for possibility of negative angles (i.e. inclination)
function Base.isapprox(orb1::Orbit, orb2::Orbit; atol=1e-6, rtol=1e-6)
    a = isapprox(orb1.a, orb2.a; atol=atol, rtol=rtol)
    e = isapprox(orb1.e, orb2.e; atol=atol, rtol=rtol)
    i = isapproxangle(orb1.i, orb2.i; atol=atol, rtol=rtol)
    M = isapproxangle(time_to_mean(0, orb1), time_to_mean(0, orb2), atol=atol, rtol=rtol)
    pos = isapprox(StateVector(0,orb1), StateVector(0,orb2); atol=atol, rtol=rtol)
    return a && e && i && M && pos
end

function Base.isapprox(stvec1::StateVector, stvec2::StateVector; atol=1e-6, rtol=1e-6)
    r = isapprox(stvec1.position, stvec2.position; atol=atol, rtol=rtol)
    v = isapprox(stvec1.velocity, stvec2.velocity; atol=atol, rtol=rtol)
    return r && v
end

# descriptive display

Base.show(io::IO, orb::Orbit) = println(io,
    "Orbit above $(orb.primary.name):\n",
    "  semi-major axis:\t\t$(orb.a) m\n",
    "  eccentricity:\t\t\t$(orb.e)\n",
    "  inclination:\t\t\t$(rad2deg(orb.i))°\n",
    "  RA of ascending node:\t\t$(rad2deg(orb.Ω))°\n",
    "  argument of periapsis:\t$(rad2deg(orb.ω))°\n",
    "  mean anomaly at epoch:\t$(rad2deg(orb.Mo))°\n",
    "  epoch:\t\t\t$(orb.epoch) s"
)
