### type definitions ###

abstract type CelestialObject end

struct Orbit{P}
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    Mo::Float64
    epoch::Float64
    primary::P
    period::Float64
    basis::RotZXZ{Float64}
end

struct StateVector
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
end

struct OrbitalState
    time::Real
    statevector::StateVector
    orbit::Orbit
end



### helper methods ###

period(a, μ) = 2π*√(abs(a)^3 / μ)


## angle-based methods ##

flight_path_angle(θ, e) = atan(e*sin(θ)/(1+e*cos(θ)))
flight_path_angle(θ, orb::Orbit) = flight_path_angle(θ, orb.e)

orbital_distance(θ, a, e) = a*(1-e^2)/(1+e*cos(θ)) 
orbital_distance(θ, orb::Orbit) = orbital_distance(θ, orb.a, orb.e)

# true anomaly of the orbit at the given distance (rising side)
orbital_angle(d, a, e) = wrap_acos((a*(1-e^2)/d - 1)/e)
orbital_angle(d, orb::Orbit) = orbital_angle(d, orb.a, orb.e)

# direction of the orbital velocity at the given true anomaly
orbital_direction(θ, e) = θ + π/2 - flight_path_angle(θ, e)
orbital_direction(d, a, e) = orbital_direction(orbital_angle(d,a,e), e)

# orbital position at the specified true anomaly
function orbital_position(θ, a, e, basis::Rotation{3,Float64})
    r = orbital_distance(θ, a, e) 
    ō = r * @SVector([cos(θ),sin(θ),0])
    return perifocal_to_inertial_bases(ō, basis)
end
orbital_position(θ, a, e, i, Ω, ω) = orbital_position(θ, a, e, basis_rotation(Ω, i, ω))
orbital_position(θ, orb::Orbit) = orbital_position(θ, orb.a, orb.e, orb.basis)

# orbital velocity at the specified true anomaly
function orbital_velocity(θ, a, e, μ, basis::Rotation{3,Float64})
    r = orbital_distance(θ, a, e) 
    ϕ = orbital_direction(θ, e)
    v = √(μ*(2/r - 1/a))
    dōdt = v * @SVector([cos(ϕ),sin(ϕ),0])
    return perifocal_to_inertial_bases(dōdt, basis)
end
orbital_velocity(θ, a, e, μ, i, Ω, ω) = orbital_velocity(θ, a, e, μ, basis_rotation(Ω, i, ω))
orbital_velocity(θ, orb::Orbit) = orbital_velocity(θ, orb.a, orb.e, orb.primary.μ, orb.basis)

state_vector(θ, a, e, μ, basis) = (orbital_position(θ, a, e, basis), orbital_velocity(θ, a, e, μ, basis))
state_vector(θ, a, e, μ, i, Ω, ω) = state_vector(θ, a, e, μ, basis_rotation(Ω, i, ω))
state_vector(θ, orb::Orbit) = state_vector(θ, orb.a, orb.e, orb.primary.μ, orb.basis)


## time-based methods ##

time_orbital_position(t, a, e, i, Ω, ω, Mo, epoch, μ, T=period(a,μ)
    ) = orbital_position(time_to_true(t, e, T, Mo, epoch), a, e, i, Ω, ω)
time_orbital_position(t, orb::Orbit
    ) = orbital_position(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.basis)

time_orbital_velocity(t, a, e, i, Ω, ω, Mo, epoch, μ, T=period(a,μ)
    ) = orbital_velocity(time_to_true(t, e, T, Mo, epoch), a, e, μ, i, Ω, ω)
time_orbital_velocity(t, orb::Orbit
    ) = orbital_velocity(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.primary.μ, orb.basis)

time_state_vector(t, a, e, i, Ω, ω, Mo, epoch, μ, T=period(a,μ)
    ) = state_vector(time_to_true(t, e, T, Mo, epoch), a, e, μ, i, Ω, ω)
time_state_vector(t, orb::Orbit
    ) = state_vector(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.primary.μ, orb.basis)



### Alternate Constructors ###

Orbit(a,body::CelestialObject
    ) = Orbit(a,0.,0.,0.,0.,0.,0.,body,period(a,body.μ),basis_rotation(0.,0.,0.))
Orbit(a,e,i,Ω,ω,Mo,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,0.,body,period(a,body.μ),basis_rotation(Ω, i, ω))
Orbit(a,e,i,Ω,ω,Mo,epoch,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,epoch,body,period(a,body.μ),basis_rotation(Ω, i, ω))

StateVector(args...) = StateVector(state_vector(args...)...)
StateVector(θ, orb::Orbit) = StateVector(state_vector(θ,orb)...)

time_StateVector(args...) = StateVector(time_state_vector(args...)...)
time_StateVector(t, orb::Orbit) = StateVector(time_state_vector(t,orb)...)

function Orbit(t, r̄::AbstractVector{<:Real}, v̄::AbstractVector{<:Real}, prim::CelestialObject, epoch=t)
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

    Ω = wrap_angle(copysign(1,n̂[2]) * wrap_acos(n̂[1]))
    ω = wrap_angle(copysign(1,ē[3]) * wrap_acos(dot(n̂,ê)))

    θ = angle_in_plane(r̄, basis_rotation(Ω, i, ω))
    M = true_to_mean(θ, e)
    Δt = t - epoch 
    if e < 1
        Mo = wrap_angle(M - time_to_mean(Δt, period(a,μ)))
    else
        Mo = M - time_to_mean(Δt, period(a,μ))
    end

    return Orbit(a,e,i,Ω,ω,Mo,epoch,prim)
end

Orbit(t, statevec::StateVector, prim::CelestialObject, epoch=t
    ) = Orbit(t, statevec.position, statevec.velocity, prim, epoch)

OrbitalState(t, orb::Orbit
    ) = OrbitalState(t, StateVector(t,orb), orb)
OrbitalState(t, statevec::StateVector, prim::CelestialObject, epoch=t
    ) = OrbitalState(t, StateVector, Orbit(t, statevec, prim, epoch))

### comparison and sorting methods ###

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
    b = isapprox(orb1.basis, orb2.basis;  atol=atol, rtol=rtol)
    pos = isapprox(time_StateVector(0,orb1), time_StateVector(0,orb2); atol=atol, rtol=rtol)
    if !a
        @show orb1.a, orb2.a
    end
    if !e
        @show orb1.e, orb2.e
    end
    if !b
        @show orb1.basis, orb2.basis
    end
    if !pos
        @show time_StateVector(0,orb1), time_StateVector(0,orb2)
    end
    return a && e && b && pos
end

function Base.isapprox(stvec1::StateVector, stvec2::StateVector; atol=1e-6, rtol=1e-6)
    r = isapprox(stvec1.position, stvec2.position; atol=atol, rtol=rtol)
    v = isapprox(stvec1.velocity, stvec2.velocity; atol=atol, rtol=rtol)
    return r && v
end

### descriptive display ###

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
