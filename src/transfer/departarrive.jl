function departarrive_direction(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel::AbstractVector{<:Real}, tₛₚ)
    # for a given parking position, and velocity at SoI, and eccentricity, define a departure/arrival orbit
    rₒ = norm(r̄ₒ)
    θₒ = orbital_angle(rₒ, a, e)
    θₛₚ = orbital_angle(rₛₚ, a, e)
    ω = bound_angle(angle_in_plane(r̄ₒ, basis_MRP(Ω, 0., i))-θₒ)
    Mₛₚ = true_to_mean(θₛₚ, e)

    Mₒ = true_to_mean(θₒ, e)
    Δt = true_to_time(θₛₚ, e, period(a, μ)) - true_to_time(θₒ, e, period(a, μ))
    tₒ = tₛₚ - Δt

    stvecₛₚ = StateVector(tₛₚ,a,e,i,Ω,ω,Mₛₚ,tₛₚ,μ)
    v̄ₛₚ = stvecₛₚ.position 

    return (a,e,i,Ω,ω,Mₒ,tₒ,), 1-dot(v̄ₛₚ./norm(v̄ₛₚ), v̄rel./norm(v̄rel))
end

function departarrive_eccentricity(θ, pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ)
    # optimize the eccentricity of the departure/arrival orbit against the departure/arrival vector
    μ = pkorb.primary.μ
    t = true_to_time(θ, pkorb)

    r̄ₒ = state_vector(t, pkorb)[1]
    rₒ = norm(r̄ₒ)
    rₛₚ = pkorb.primary.SoI
    h̄ = cross(r̄ₒ,v̄rel)
    n̄ = cross(@SVector([0.,0.,1.]), h̄)
    n̂ = n̄/norm(n̄)

    a = 1 / (2/rₛₚ - sum(abs2,v̄rel)/μ)
    i = acos(h̄[3]/norm(h̄))
    Ω = bound_angle(copysign(1,n̂[2]) * acos(n̂[1]))

    lb = a > 0 ? 0. : 1.
    ub = max(1-rₒ/a, rₒ/a-1)
    ub = ub == 1 ? 1. + eps(Float64) : ub

    # objectivefun(e) = departarrive_direction(e[1], a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, tₛₚ)[2]
    # res = optimize(objectivefun, [lb], [ub], [(lb+ub)/2], Fminbox(LBFGS()))

    objectivefun(e) = departarrive_direction(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, tₛₚ)[2]
    res = optimize(objectivefun, lb, ub, Brent())
    e = Optim.minimizer(res)

    (a,e,i,Ω,ω,Mₒ,tₒ), eccobj = departarrive_direction(e, a, Ω, i, r̄ₒ, rₛₚ, μ, v̄rel, tₛₚ)
    daorb = Orbit(a,e,i,Ω,ω,Mₒ,tₒ,pkorb.primary)
    stvecₒ = StateVector(tₒ, daorb)
    stvecₚₖ = StateVector(tₒ, pkorb)

    # get the Δv at the burn for this eccentricity
    Δv̄ = stvecₚₖ.velocity .- stvecₒ.velocity

    return daorb, eccobj, Δv̄
end

function departarrive_orbit(pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ::Real)
    # optimize dv against true anomaly, with each true anomaly being optimized for eccentricity
    if pkorb.e < 1           # elliptical case: any true anomaly
        lb = -1*π
        ub = π
    else                    # hyperbolic case: true anomalies occuring within the SoI
        ub = -lb = orbital_angle(pkorb.primary.SoI, pkorb.a, pkorb.e)
    end
    
    @show lb, ub

    # function objectivefun(θ)
    #     eccobj, Δv̄ = departarrive_eccentricity(θ[1], pkorb, v̄rel, tₛₚ)[2:3]
    #     return eccobj <= atol ? norm(Δv̄) : Inf
    # end

    # res = optimize(objectivefun, [lb], [ub], [0.], Fminbox(LBFGS()))
    # θₚₖ = Optim.minimizer(res)[1]

    function objectivefun(θ)
        eccobj, Δv̄ = departarrive_eccentricity(θ, pkorb, v̄rel, tₛₚ)[2:3]
        return exp(eccobj*1e6) * norm(Δv̄)
    end

    res = optimize(objectivefun, lb, ub, Brent())
    θₚₖ = Optim.minimizer(res)

    daorb, eccobj, Δv̄ = departarrive_eccentricity(θₚₖ, pkorb, v̄rel, tₛₚ)
    return daorb, eccobj, Δv̄
end

function depart_arrive_orbit()

end

function quick_departarrive_orbit(pkorb::Orbit, v̄rel::AbstractVector{<:Real}, tₛₚ; trajectory=:out, tol=0.1, maxit=50)
    if trajectory == :out
        c = 1
    elseif trajectory == :in
        c = -1
    else
        throw(ArgumentError("`trajectory` of "*string(trajectory)*" is invalid (must be either `:out` or `:in`)."))
    end
    μ = pkorb.primary.μ
    rₛₚ = pkorb.primary.SoI
    vₛₚ = norm(v̄rel)
    err = tol+1
    it = 0
    rₒnext = pkorb.a
    a = 1 / (2/rₛₚ - vₛₚ^2/μ)
    @show a
    e, θₛₚ = 0., 0.
    r̄ₒ = @SVector zeros(Float64, 3)
    v̄ₒ = @SVector zeros(Float64, 3)
    while abs(err) > tol && it < maxit
        it = it+1
        @show it
        rₒ = rₒnext
        vₒ = √(vₛₚ^2 + 2(μ/rₒ - μ/rₛₚ))
        e = √(1 + 2((vₒ^2)/2 - μ/rₒ) * rₒ^2 * vₒ^2 / μ^2)
        θₛₚ = c * orbital_angle(rₛₚ, a, e)
        @show rₛₚ, a, e
        ϕₛₚ = flight_path_angle(θₛₚ, e)
        v̄ₛₚ = vₛₚ * @SVector([cos(θₛₚ + π/2 - ϕₛₚ),
                            sin(θₛₚ + π/2 - ϕₛₚ),
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
        r̄ₒpkorb = state_vector(true_to_time(δθ, pkorb), pkorb)[1]
        errprev = err
        err = norm(r̄ₒ) - norm(r̄ₒpkorb)
        if abs(err/errprev) > 0.9
            rₒnext = (norm(r̄ₒpkorb) + rₒ)/2
        else
            rₒnext = norm(r̄ₒpkorb)
        end
    end
    @show θₛₚ, e, period(a, μ)
    Δt = true_to_time(θₛₚ, e, period(a, μ))
    @show Δt
    return Orbit(tₛₚ - Δt, r̄ₒ, v̄ₒ, pkorb.primary)
end

departure_orbit(pkorb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...) = departarrive_orbit(pkorb, v̄rel, tₛₚ; trajectory=:out, kwargs...)
arrive_orbit(pkorb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...) =    departarrive_orbit(pkorb, v̄rel, tₛₚ; trajectory=:in, kwargs...)
