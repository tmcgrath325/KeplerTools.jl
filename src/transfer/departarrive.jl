function departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel::AbstractVector{<:Real}, c)
    # for a given parking position, and velocity at SoI, and eccentricity, define a departure/arrival orbit
    rₒ = norm(r̄ₒ)
    θₒ = c * orbital_angle(rₒ, a, e)
    θₛₚ = c * orbital_angle(rₛₚ, a, e)
    ω = wrap_angle(angle_in_plane(r̄ₒ, basis_MRP(Ω, 0., i))-θₒ)

    v̄ₛₚ = orbital_velocity(θₛₚ,a,e,μ,i,Ω,ω)

    return (a,e,i,Ω,ω,), 1-dot(v̄ₛₚ./norm(v̄ₛₚ), v̄rel./norm(v̄rel))
end

function departarrive_true_anomaly_objective(θ, pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ, c)
    # optimize the eccentricity of the departure/arrival orbit against the departure/arrival vector
    μ = pkorb.primary.μ

    r̄ₒ = state_vector(θ, pkorb)[1]
    rₒ = norm(r̄ₒ)
    rₛₚ = pkorb.primary.SoI
    h̄ = cross(r̄ₒ,v̄rel)
    n̄ = cross(@SVector([0.,0.,1.]), h̄)
    n̂ = n̄/norm(n̄)

    a = 1 / (2/rₛₚ - sum(abs2,v̄rel)/μ)
    i = acos(h̄[3]/norm(h̄))
    Ω = wrap_angle(copysign(1,n̂[2]) * acos(n̂[1]))

    if a > 0
        lb = 0.
        ub = 1. - eps(Float64)
    else
        lb = 1. + eps(Float64)
        ub = max(1-rₒ/a, rₒ/a-1)
        ub = ub == 1 ? 1. + 2*eps(Float64) : ub
    end

    objectivefun(e) = departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, c)[2]
    res = optimize(objectivefun, lb, ub, Brent())
    e = Optim.minimizer(res)

    (a,e,i,Ω,ω,), eccobj = departarrive_eccentricity_objective(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, c)
    θₒ = c * orbital_angle(rₒ, a, e)
    θₛₚ = c * orbital_angle(rₛₚ, a, e)
    Mₒ = true_to_mean(θₒ, e)
    T = period(a,μ)
    Δt = true_to_time(θₛₚ, e, T) - true_to_time(θₒ, e, T)
    tₒ = true_to_time(θ, pkorb, tₛₚ - Δt - pkorb.period/2)  # ensure `pkorb` and `daorb` are continuous
    # tₒ = tₛₚ - Δt                                           # ensure `tₛₚ` is matched exactly
    
    daorb = Orbit(a,e,i,Ω,ω,Mₒ,tₒ,pkorb.primary)

    # get the Δv at the burn for this eccentricity
    v̄ₒ = orbital_velocity(θₒ, daorb)
    v̄ₚₖ = orbital_velocity(θ, pkorb)
    Δv̄ = v̄ₒ .- v̄ₚₖ
    return daorb, eccobj, Δv̄
end

function departarrive_orbit(pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ::Real; out=true)
    c = out ? 1 : -1
    # optimize dv against true anomaly, with each true anomaly being optimized for eccentricity
    if pkorb.e < 1           # elliptical case: any true anomaly
        lb = -2π
        ub = 2π
    else                    # hyperbolic case: true anomalies occuring within the SoI
        ub = -lb = orbital_angle(pkorb.primary.SoI, pkorb.a, pkorb.e) - eps(Float64)
    end

    # objectivefun(θ) = norm(departarrive_true_anomaly_objective(θ, pkorb, v̄rel, tₛₚ, c)[3])
    function objectivefun(θ)
        err, Δv̄ = departarrive_true_anomaly_objective(θ, pkorb, v̄rel, tₛₚ, c)[2:3]
        return norm(Δv̄) * exp(1000*err)        # penalize mismatched direction
    end

    res = optimize(objectivefun, lb, ub, Brent())
    θₚₖ = Optim.minimizer(res)

    daorb, eccobj, Δv̄ = departarrive_true_anomaly_objective(θₚₖ, pkorb, v̄rel, tₛₚ, c)
    if eccobj > 1e-3
        @warn "departure/arrival velocity mismatch from target direction: $(pkorb.primary.name), $eccobj"
    end
    return daorb, eccobj, Δv̄
end

function quick_departarrive_orbit(pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ; out=true, tol=0.1, maxit=50)
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
        ψ₁ = atan(v̄relplane[3], √(sum(abs2, v̄ₛₚ) - v̄ₛₚ[1]^2 - v̄relplane[3]^2))
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
    v̄ₚₖ = orbital_velocity(θₚₖ, pkorb)
    Δv̄ = v̄ₒ - v̄ₚₖ
    return Orbit(tₛₚ - Δt, r̄ₒ, v̄ₒ, pkorb.primary), Δv̄
end

departure_orbit(pkorb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...) = departarrive_orbit(pkorb, v̄rel, tₛₚ; out=true,  kwargs...)
arrival_orbit(pkorb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...)   = departarrive_orbit(pkorb, v̄rel, tₛₚ; out=false, kwargs...)
