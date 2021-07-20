function departarrive_orbit(porb::Orbit, v̄rel::SVector{3}, tₛₚ; trajectory=:out, tol=0.1, maxit=50)
    if trajectory == :out
        c = 1
    elseif trajectory == :in
        c = -1
    else
        throw(ArgumentError("`trajectory` of "*string(trajectory)*" is invalid (must be either `:out` or `:in`)."))
    end
    μ = porb.primary.μ
    rₛₚ = porb.primary.SoI
    err = tol+1
    it = 1
    rₒnext = porb.a
    while abs(err) > tol && it < maxit
        it = it+1
        rₒ = rₒnext
        vₒ = √(sum(abs2(v̄rel)) +2(μ/rₒ - μ/rₛₚ))
        e = √(1 + 2((vₒ^2)/2 - μ/rₒ) * rₒ^2 * vₒ^2 / μ^2)
        a = 1 / (2/rₒ - vₒ^2/μ)
        θₛₚ = c * acos(1/e * a*(1-e^2)/rₛₚ)
        ϕₛₚ = atan(e*sin(θₛₚ) / (1 + e*cos(θₛₚ)))
        v̄ₛₚ = √(μ * (2/rₛₚ - 1/a)) * @SVector([cos(θₛₚ + π/2 - ϕₛₚ),
                                             sin(θₛₚ + π/2 - ϕₛₚ),
                                             0])
        r̄ₒ = @SVector([rₒ, 0, 0])
        v̄ₒ = @SVector([0, vₒ, 0])
        v̄relplane = ref_to_plane_bases(porb, v̄rel)
        ψ₁ = atan(v̄relplane[3], √(sum(abs2(v̄ₛₚ)) - v̄ₛₚ[1]^2 - v̄relplane[3]^2))
        R₁ = MRP(@SMatrix([1        0        0;                   # rotate around x-axis to match i
                           0  cos(ψ₁) -sin(ψ₁);
                           0  sin(ψ₁)  cos(ψ₁)]))
        v̄ₒ = R₁*v̄ₒ
        ψ₂ = atan(v̄relplane[2], v̄relplane[1]) - atan(v̄ₒ[2], v̄ₒ[1])
        R₂ = MRP(@SMatrix([cos(ψ₂) -sin(ψ₂)  0;                   # rotate around z-axis to match ω
                           sin(ψ₂)  cos(ψ₂)  0;
                           0        0        1]))
        v̄ₒ = plane_to_ref_bases(porb, R₂*v̄ₒ)
        r̄ₒ = plane_to_ref_bases(porb, R₂*R₁*r̄ₒ)

        δθ = angle_in_plane(porb, r̄ₒ)
        r̄ₒporb = state_vector(porb, true_to_time(porb, δθ))
        errprev = err
        err = norm(r̄ₒ) - norm(r̄ₒporb)
        if abs(err/errprev) > 0.9
            rₒnext = (norm(r̄ₒporb) + rₒ)/2
        else
            rₒnext = norm(r̄ₒporb)
        end
    end
    # k̂ = plane_to_ref_bases(porb, @SVector([0 0 1]))
    # v̄ₚ = v(μ/rₒ) cross(k̂,r̄ₒ)/rₒ
    Δt = true_to_time(θₛₚ, e, period(a, μ))
    return Orbit(r̄ₒ, v̄ₒ, tₛₚ - Δt, primary(porb))
end

departure_orbit(porb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...) = departarrive_orbit(porb, v̄rel, tₛₚ; trajectory=:out, kwargs...)
arrive_orbit(porb::Orbit, v̄rel::SVector{3}, tₛₚ; kwargs...) =    departarrive_orbit(porb, v̄rel, tₛₚ; trajectory=:in, kwargs...)
