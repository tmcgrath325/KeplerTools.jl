function departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel::AbstractVector{<:Real}, c)
    # for a given parking position, and velocity at SoI, and eccentricity, define a departure/arrival orbit
    rₒ = norm(r̄ₒ)
    θₒ = c * orbital_angle(rₒ, a, e)
    θₛₚ = c * orbital_angle(rₛₚ, a, e)
    ω = wrap_angle(angle_in_plane(r̄ₒ, basis_rotation(Ω, i, 0.))-θₒ)

    v̄ₛₚ = orbital_velocity(θₛₚ,a,e,μ,i,Ω,ω)

    return 1-dot(v̄ₛₚ./norm(v̄ₛₚ), v̄rel./norm(v̄rel)), ω
end

function departarrive_true_anomaly(r̄ₒ::AbstractVector{<:Real}, prim::Union{CelestialBody,Star}, v̄rel::AbstractVector{<:Real}, tₛₚ, c)
    # optimize the eccentricity of the departure/arrival orbit against the departure/arrival vector
    μ = prim.μ
    rₒ = norm(r̄ₒ)
    rₛₚ = prim.SoI
    ĥ = normalize(cross(r̄ₒ,v̄rel))
    n̂ = normalize(cross(@SVector([0.,0.,1.]), ĥ))

    a =  1 / (2/rₛₚ - sum(abs2,v̄rel)/μ)
    i = wrap_acos(ĥ[3])
    Ω = wrap_angle(copysign(1,n̂[2]) * acos(n̂[1]))

    if a > 0
        lb = 0.
        ub = 1. - eps(Float64)
    else
        lb = 1. + eps(Float64)
        ub = max(1-rₒ/a, rₒ/a-1)
        ub = ub == 1 ? 1. + 2*eps(Float64) : ub
    end

    objectivefun(e) = departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, c)[1]
    res = optimize(objectivefun, lb, ub, Brent())
    e = Optim.minimizer(res)

    eccobj, ω = departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, c)

    θₒ = c * orbital_angle(rₒ, a, e)
    θₛₚ = c * orbital_angle(rₛₚ, a, e)
    T = period(a,μ)
    Mₒ = true_to_mean(θₒ, e)
    Δt = true_to_time(θₛₚ, e, T) - true_to_time(θₒ, e, T)
    tₒ = tₛₚ - Δt

    return Orbit(a,e,i,Ω,ω,Mₒ,tₒ,prim), eccobj, 0.
end

function departarrive_true_anomaly(θ, pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ, c)
    r̄ₒ = orbital_position(θ, pkorb)
    daorb, eccobj = departarrive_true_anomaly(r̄ₒ, pkorb.primary, v̄rel::AbstractVector{<:Real}, tₛₚ, c)[1:2]

    a, e, i, Ω, ω, Mₒ = daorb.a, daorb.e, daorb.i, daorb.Ω, daorb.ω, daorb.Mo

    tₒ = true_to_time(θ, pkorb, daorb.epoch - pkorb.period/2)

    v̄ₒ = time_orbital_velocity(daorb.epoch, daorb)
    v̄ₚₖ = orbital_velocity(θ, pkorb)
    Δv̄ = v̄ₒ .- v̄ₚₖ

    return Orbit(a,e,i,Ω,ω,Mₒ,tₒ,pkorb.primary), eccobj, Δv̄
end

function departarrive_objectivefun(θ, pkorb, v̄rel, tₛₚ, c)
    err, Δv̄ = departarrive_true_anomaly(θ, pkorb, v̄rel, tₛₚ, c)[2:3]
    return norm(Δv̄) * exp(1000*err)        # penalize mismatched direction
end

function departarrive_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; out=true)
    c = out ? 1 : -1
    # optimize dv against true anomaly, with each true anomaly being optimized for eccentricity
    if pkorb.e < 1           # elliptical case: any true anomaly
        lb = -2π
        ub = 2π
    else                    # hyperbolic case: true anomalies occuring within the SoI
        ub = orbital_angle(pkorb.primary.SoI, pkorb.a, pkorb.e) - eps(Float64)
        lb = -ub
    end

    objectivefun(θ) = departarrive_objectivefun(θ, pkorb, v̄rel, tₛₚ, c)

    res = optimize(objectivefun, lb, ub, Brent())
    θₚₖ = Optim.minimizer(res)

    daorb, eccobj, Δv̄ = departarrive_true_anomaly(θₚₖ, pkorb, v̄rel, tₛₚ, c)
    if eccobj > 1e-3
        @warn "departure/arrival velocity mismatch from target direction: $(pkorb.primary.name), $eccobj"
    end
    return daorb, eccobj, Δv̄
end

function departarrive_orbit(r̄ₒ::AbstractVector{<:Real}, tₒ, prim::Union{CelestialBody,Star}, v̄rel::AbstractVector{<:Real}, tₛₚ; out=true)
    # return daorb, 0, Δv̄
end


function fast_departarrive_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; out=true, tol=0.1, maxit=50)
    c = out ? 1 : -1
    μ = pkorb.primary.μ
    rₛₚ = pkorb.primary.SoI
    vₛₚ = norm(v̄rel)
    err = tol+1
    it = 0
    rₒnext = pkorb.a
    a = 1 / (2/rₛₚ - vₛₚ^2/μ)

    e, θₛₚ = 0., 0.
    r̄ₒ = @SVector zeros(Float64, 3)
    v̄ₒ = @SVector zeros(Float64, 3)
    while abs(err) > tol && it < maxit
        it = it+1
        rₒ = rₒnext
        vₒ = √(vₛₚ^2 + 2(μ/rₒ - μ/rₛₚ))
        e = √(1 + 2((vₒ^2)/2 - μ/rₒ) * rₒ^2 * vₒ^2 / μ^2)
        θₛₚ = c * orbital_angle(rₛₚ, a, e)
        ϕₛₚ = orbital_direction(θₛₚ, e)
        v̄ₛₚ = vₛₚ * @SVector([cos(ϕₛₚ),
                            sin(ϕₛₚ),
                            0])

        r̄ₒ = @SVector([rₒ, 0, 0])
        v̄ₒ = @SVector([0, vₒ, 0])
        v̄relplane = inertial_to_perifocal_bases(v̄rel, pkorb)
        # R = align_vectors(v̄ₛₚ, v̄relplane)               # align velocity at SoI with input direction,
        #                                                # valid only for undefined parking orbit
        ψ₁ = atan(v̄relplane[3], √abs(vₛₚ^2 - v̄ₛₚ[1]^2 - v̄relplane[3]^2))
        # R₁ = MRP([1        0        0;                   # rotate around x-axis to match i
        #           0  cos(ψ₁) -sin(ψ₁);
        #           0  sin(ψ₁)  cos(ψ₁)])
        R₁ = MRP((1.,0.,0.,0.,cos(ψ₁),sin(ψ₁),0.,-sin(ψ₁),cos(ψ₁)))
        R1v̄ₛₚ = R₁*v̄ₛₚ
        ψ₂ = atan(v̄relplane[2], v̄relplane[1]) - atan(R1v̄ₛₚ[2], R1v̄ₛₚ[1])
        # R₂ = MRP([cos(ψ₂) -sin(ψ₂)  0;                   # rotate around z-axis to match ω
        #           sin(ψ₂)  cos(ψ₂)  0;
        #           0        0        1]))
        R₂ = MRP((cos(ψ₂),sin(ψ₂),0.,-sin(ψ₂),cos(ψ₂),0.,0.,0.,1.))

        v̄ₒ = perifocal_to_inertial_bases(R₂*R₁*v̄ₒ, pkorb)
        r̄ₒ = perifocal_to_inertial_bases(R₂*R₁*r̄ₒ, pkorb)

        δθ = angle_in_plane(r̄ₒ, pkorb)
        r̄ₒpkorb = state_vector(δθ, pkorb)[1]

        errprev = err
        err = norm(r̄ₒ) - norm(r̄ₒpkorb)
        if abs(err/errprev) > 0.9
            rₒnext = (norm(r̄ₒpkorb) + rₒ)/2
        else
            rₒnext = norm(r̄ₒpkorb)
        end
    end
    Δt = true_to_time(θₛₚ, e, period(a, μ))
    θₚₖ = angle_in_plane(r̄ₒ, pkorb)
    tₚₖ = true_to_time(θₚₖ, pkorb, tₛₚ-Δt-pkorb.period/2)
    v̄ₚₖ = orbital_velocity(θₚₖ, pkorb)
    Δv̄ = c*(v̄ₒ .- v̄ₚₖ)
    return Orbit(tₚₖ, r̄ₒ, v̄ₒ, pkorb.primary), Δv̄
end


"""
    orb = departure_orbit(pkorb, v̄rel, tₛₚ)

Computes an optimized departure orbit from the parking orbit `pkorb`, the relative velocity
at the patch `v̄rel`, the time at the patch `tₛₚ`. 
"""
departure_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; kwargs...) = departarrive_orbit(pkorb, v̄rel, tₛₚ; out=true,  kwargs...)

"""
    orb = arrival_orbit(pkorb, v̄rel, tₛₚ)

Computes an optimized arrival orbit from the parking orbit `pkorb`, the relative velocity
at the patch `v̄rel`, the time at the patch `tₛₚ`.
"""
arrival_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; kwargs...)   = departarrive_orbit(pkorb, v̄rel, tₛₚ; out=false, kwargs...)


"""
    orb = fast_departure_orbit(pkorb, v̄rel, tₛₚ)

Computes a departure orbit from the parking orbit `pkorb`, the relative velocity
at the patch `v̄rel`, the time at the patch `tₛₚ`. 

This function is faster than `departure_orbit`, but may return sub-optimal trajectories
for highly elliptical parking orbits. 
"""
fast_departure_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; kwargs...) = fast_departarrive_orbit(pkorb, v̄rel, tₛₚ; out=true, kwargs...)

"""
    orb = fast_arrival_orbit(pkorb, v̄rel, tₛₚ)

Computes a arrival orbit from the parking orbit `pkorb`, the relative velocity
at the patch `v̄rel`, the time at the patch `tₛₚ`. 

This function is faster than `arrival_orbit`, but may return sub-optimal trajectories
for highly elliptical parking orbits. 
"""
fast_arrival_orbit(pkorb::Orbit{<:CelestialBody}, v̄rel::AbstractVector{<:Real}, tₛₚ; kwargs...)   = fast_departarrive_orbit(pkorb, v̄rel, tₛₚ; out=false, kwargs...)


# true anomaly/time at the SoI
patch_angle(daorb::Orbit, c) = c*orbital_angle(daorb.primary.SoI, daorb)
patch_time(daorb::Orbit, c) = true_to_time(patch_angle(daorb, c), daorb)

"""
    θ = ejection_angle(dorb)

Computes the true anomaly at SoI ejection for a departure orbit `dorb`.
"""
ejection_angle(dorb::Orbit)  = patch_angle(dorb,  1)

"""
    t = ejection_time(dorb)

Computes the time at SoI ejection for a departure orbit `dorb`.
"""
ejection_time(dorb::Orbit)   = patch_time(dorb,  1)

"""
    θ = insertion_angle(aorb)

Computes the true anomaly at SoI insertion for a arrival orbit `aorb`.
"""
insertion_angle(aorb::Orbit) = patch_angle(aorb, -1)

"""
    t = insertion_angle(aorb)

Computes the time at SoI insertion for a arrival orbit `aorb`.
"""
insertion_time(aorb::Orbit)  = patch_time(aorb, -1)

patch_position(daorb::Orbit, c) = orbital_position(patch_angle(daorb, c), daorb)
patch_velocity(daorb::Orbit, c) = orbital_velocity(patch_angle(daorb, c), daorb)
patch_state_vector(daorb::Orbit, c) = state_vector(patch_angle(daorb, c), daorb)

"""
    r̄ = ejection_position(dorb)

Computes the position vector at SoI ejection for a departure orbit `dorb`.
"""
ejection_position(dorb::Orbit) = orbital_position(ejection_angle(dorb), dorb)

"""
    v̄ = ejection_velocity(dorb)

Computes the velocity vector at SoI ejection for a departure orbit `dorb`.
"""
ejection_velocity(dorb::Orbit) = orbital_velocity(ejection_angle(dorb), dorb)

"""
    r̄, v̄ = ejection_state_vector(dorb)

Computes the state vector at SoI ejection for a departure orbit `dorb`.
"""
ejection_state_vector(dorb::Orbit) = state_vector(ejection_angle(dorb), dorb)


"""
    r̄ = insertion_position(aorb)

Computes the position vector at SoI insertion for a arrival orbit `aorb`.
"""
insertion_position(aorb::Orbit) = orbital_position(insertion_angle(aorb), aorb)


"""
    v̄= insertion_velocity(aorb)

Computes the velocity vector at SoI insertion for a arrival orbit `aorb`.
"""
insertion_velocity(aorb::Orbit) = orbital_velocity(insertion_angle(aorb), aorb)


"""
    r̄, v̄ = insertion_state_vector(aorb)

Computes the state vector at SoI insertion for a arrival orbit `aorb`.
"""
insertion_state_vector(aorb::Orbit) = state_vector(insertion_angle(aorb), aorb)

    