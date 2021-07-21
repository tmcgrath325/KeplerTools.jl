function p_lambert(r̄ₛ::SVector{3,<:Real}, r̄ₑ::SVector{3,<:Real}, Δt::Real, μ::Real, tol=1e-6::Real, maxit=200::Int)
    rₛ, rₑ = norm(r̄ₛ), norm(r̄ₑ)

    # true anomaly change for the transfer
    Δθ = angle_in_plane(r̄ₑ, r̄ₛ, MRP(1I))

    # p iteration constants
    k = rₛ * rₑ * (1-cos(Δθ))
    L = rₛ + rₑ
    m = rₛ * rₑ * (1+cos(Δθ))

    # bounds
    pj = k / (L + √(2m))
    pjj = k / (L - √(2m))
    if Δθ > π
        pmin = 0
        pmax = pjj
    else
        pmin = pj
        pmax = Inf
    end

    # Newton-p-iteration
    it = 0
    err = tol + 1
    p = (pj+pjj)/2
    pnext = p
    f = g = 0.
    while (err > tol) && (it < maxit)
        it = it + 1
        p = pnext
        a = m*k*p / ((2m-L^2)*(p^2) + 2k*L*p - k^2)
        f = 1 - rₑ/p * (1 - cos(Δθ))
        g = rₛ * rₑ * sin(Δθ) / √(μ*p)
        df = √(μ/p)*tan(Δθ/2)*((1-cos(Δθ))/p - 1/rₛ - 1/rₑ);
        if a > 0    # Elliptical case
            sinΔE = -rₛ*rₑ*df/√(μ*a)
            cosΔE = 1 - rₛ/a * (1-f)
            ΔE = atan(sinΔE, cosΔE)
            t = g + √(a^3/μ)*(ΔE-sinΔE)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) + √(a^3/μ)*(2k*sinΔE)/(p*(k-L*p))
        else        # Hyperbolic case
            ΔF = acosh(1 - rₛ/a * (1-f))
            t = g + √(-a^3/μ)*(sinh(ΔF)-ΔF)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) - √(-a^3/μ)*(2k*sinh(ΔF))/(p*(k-L*p))
        end
        err = abs(Δt - t)/Δt
        @show err
        pnext = p + (Δt - t)/dtdp
        # If the next guess is outside of allowed bounds, use bisection
        if pnext < pmin
            pnext = (p + pmin)/2
        elseif pnext > pmax
            pnext = (p + pmax)/2
        end
    end
    if it == maxit
        error("Maximum number of iterations exceeded: failed to converge")
    end
    v̄ = (r̄ₑ - f*r̄ₛ)/g
    return v̄
end

function p_lambert(sorb::Orbit, eorb::Orbit, stime, etime; kwargs...) 
    r̄ₛ, r̄ₑ = state_vector(stime,sorb)[1], state_vector(etime,eorb)[1]
    v̄ = p_lambert(r̄ₛ, r̄ₑ, etime-stime, sorb.primary.μ; kwargs...)
    return Orbit(stime, r̄ₛ, v̄, sorb.primary)
end

