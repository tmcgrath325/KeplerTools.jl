### type definitions ###

"""
A celestial object, such as a planet, star, or moon.
"""
abstract type CelestialObject end

"""
A Keplerian orbit, defined, by semimajor axis (`.a`), eccentricity (`.e`), inclination(`.i`),
right ascension of the ascending node (`.Ω`), argument of the periapsis (`.ω`), 
mean anomaly at epoch(`.Mo`), and an epoch (`.epoch`), as well as the celestial body around
which the orbit takes place (.`primary`).
"""
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

### helper methods ###
"""
    Δt = period(a, μ)

Compute the period of an orbit with semimajor axis `a` around a body with
gravity parameter `μ`.
"""
period(a, μ) = 2π*√(abs(a)^3 / μ)


## angle-based methods ##
"""
    ϕ = flight_path_angle(θ, e)
    ϕ = flight_path_angle(θ, orbit)

Compute the flight path angle of an orbit with eccentricity `e` at true anomaly `θ`.
"""
flight_path_angle(θ, e) = atan(e*sin(θ)/(1+e*cos(θ)))
flight_path_angle(θ, orb::Orbit) = flight_path_angle(θ, orb.e)

"""
    r = orbital_distance(θ, a, e)
    r = orbital_distance(θ, orbit)

Compute the orbital distance of an orbit with semimajor axis `a` and eccentricity `e` 
at true anomaly `θ`, via the [conic-section orbit equation](https://en.wikipedia.org/wiki/Kepler_orbit#Development_of_the_laws).
"""
orbital_distance(θ, a, e) = a*(1-e^2)/(1+e*cos(θ)) 
orbital_distance(θ, orb::Orbit) = orbital_distance(θ, orb.a, orb.e)

"""
    θ = orbital_angle(r, a, e)
    θ = orbital_angle(r, orbit)

Compute the true anomaly of an orbit with semimajor axis `a` and eccentricity `e` 
at orbital distance `r`, via the [conic-section orbit equation](https://en.wikipedia.org/wiki/Kepler_orbit#Development_of_the_laws).
"""
orbital_angle(r, a, e) = wrap_acos((a*(1-e^2)/r - 1)/e)
orbital_angle(r, orb::Orbit) = orbital_angle(r, orb.a, orb.e)

"""
    ψ = orbital_direction(θ, e)
    ψ = orbital_direction(r, a, e)

Compute the direction of velocity of an orbit with semimajor axis `a` at true anomaly `θ`.

Alternatively to true anomaly, orbital distance `r` and eccentricity `e` can be provided.
"""
orbital_direction(θ, e) = θ + π/2 - flight_path_angle(θ, e)
orbital_direction(r, a, e) = orbital_direction(orbital_angle(r,a,e), e)

"""
    r̄ = orbital_position(θ, a, e, basis)
    r̄ = orbital_position(θ, a, e, i, Ω, ω)
    r̄ = orbital_position(θ, orbit)

Compute the position vector at true anomaly `θ` of an orbit with with semimajor axis `a`, 
eccentricity `e`, and a rotation `basis` defining the orbit's perifocal plane.

The rotation can be computed from inclinaion `i`, RAAN `Ω`, and the argument of the periapsis `ω`.
"""
function orbital_position(θ, a, e, basis::Rotation{3,Float64})
    r = orbital_distance(θ, a, e) 
    ō = r * @SVector([cos(θ),sin(θ),0])
    return perifocal_to_inertial_bases(ō, basis)
end
orbital_position(θ, a, e, i, Ω, ω) = orbital_position(θ, a, e, basis_rotation(Ω, i, ω))
orbital_position(θ, orb::Orbit) = orbital_position(θ, orb.a, orb.e, orb.basis)

"""
    v̄ = orbital_velocity(θ, a, e, μ, basis)
    v̄ = orbital_velocity(θ, a, e, μ, i, Ω, ω)
    v̄ = orbital_velocity(θ, orbit)

Compute the velocity vector at true anomaly `θ` of an orbit with with semimajor axis `a`, 
eccentricity `e`, gravity paremter `μ`, and a rotation `basis` defining the orbit's perifocal plane.

The rotation can be computed from inclinaion `i`, RAAN `Ω`, and the argument of the periapsis `ω`.
"""
function orbital_velocity(θ, a, e, μ, basis::Rotation{3,Float64})
    r = orbital_distance(θ, a, e) 
    ϕ = orbital_direction(θ, e)
    v = √(μ*(2/r - 1/a))
    dōdt = v * @SVector([cos(ϕ),sin(ϕ),0])
    return perifocal_to_inertial_bases(dōdt, basis)
end
orbital_velocity(θ, a, e, μ, i, Ω, ω) = orbital_velocity(θ, a, e, μ, basis_rotation(Ω, i, ω))
orbital_velocity(θ, orb::Orbit) = orbital_velocity(θ, orb.a, orb.e, orb.primary.μ, orb.basis)

"""
    r̄, v̄ = state_vector(θ, a, e, μ, basis)
    r̄, v̄ = state_vector(θ, a, e, μ, i, Ω, ω)
    r̄, v̄ = state_vector(θ, orbit)

Compute the state vector at true anomaly `θ` of an orbit with with semimajor axis `a`, 
eccentricity `e`, gravity paremter `μ`, and a rotation `basis` defining the orbit's perifocal plane.

The rotation can be computed from inclinaion `i`, RAAN `Ω`, and the argument of the periapsis `ω`.
"""
state_vector(θ, a, e, μ, basis) = (orbital_position(θ, a, e, basis), orbital_velocity(θ, a, e, μ, basis))
state_vector(θ, a, e, μ, i, Ω, ω) = state_vector(θ, a, e, μ, basis_rotation(Ω, i, ω))
state_vector(θ, orb::Orbit) = state_vector(θ, orb.a, orb.e, orb.primary.μ, orb.basis)


## time-based methods ##

"""
    r̄ = time_orbital_position(t, a, e, i, Ω, ω, Mo, epoch, μ)
    r̄ = time_orbital_position(t, orbit)

Compute the position vector at time `t` from a full set of orbital paremeters, as well as a gravity parameter `μ`.
"""
time_orbital_position(t, a, e, i, Ω, ω, Mo, epoch, μ
    ) = orbital_position(time_to_true(t, e, period(a,μ), Mo, epoch), a, e, i, Ω, ω)
time_orbital_position(t, orb::Orbit
    ) = orbital_position(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.basis)

"""
    r̄ = time_orbital_velocity(t, a, e, i, Ω, ω, Mo, epoch, μ)
    r̄ = time_orbital_velocity(t, orbit)

Compute the velocity vector at time `t` from a full set of orbital paremeters, as well as a gravity parameter `μ`.
"""
time_orbital_velocity(t, a, e, i, Ω, ω, Mo, epoch, μ
    ) = orbital_velocity(time_to_true(t, e, period(a,μ), Mo, epoch), a, e, μ, i, Ω, ω)
time_orbital_velocity(t, orb::Orbit
    ) = orbital_velocity(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.primary.μ, orb.basis)

"""
    r̄ = time_orbital_velocity(t, a, e, i, Ω, ω, Mo, epoch, μ)
    r̄ = time_orbital_velocity(t, orbit)

Compute the state vector at time `t` from a full set of orbital paremeters, as well as a gravity parameter `μ`.
"""
time_state_vector(t, a, e, i, Ω, ω, Mo, epoch, μ
    ) = state_vector(time_to_true(t, e, period(a,μ), Mo, epoch), a, e, μ, i, Ω, ω)
time_state_vector(t, orb::Orbit
    ) = state_vector(time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch), orb.a, orb.e, orb.primary.μ, orb.basis)


### Alternate Constructors ###

Orbit(a,body::CelestialObject
    ) = Orbit(a,0.,0.,0.,0.,0.,0.,body,period(a,body.μ),basis_rotation(0.,0.,0.))
Orbit(a,e,i,Ω,ω,Mo,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,0.,body,period(a,body.μ),basis_rotation(Ω, i, ω))
Orbit(a,e,i,Ω,ω,Mo,epoch,body::CelestialObject
    ) = Orbit(a,e,i,Ω,ω,Mo,epoch,body,period(a,body.μ),basis_rotation(Ω, i, ω))

"""
    orb = Orbit(t, r̄, v̄, prim, epoch)

Compute an orbit from a time `t`, state vector (`r̄, v̄`), primary body `prim`, and epoch.
"""
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
    pos = isapprox(time_orbital_position(0,orb1), time_orbital_position(0,orb2); atol=atol, rtol=rtol)
    return a && e && b && pos
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
