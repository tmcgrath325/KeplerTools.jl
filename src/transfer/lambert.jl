function p_lambert(r̄ₛ::AbstractVector{<:Real}, r̄ₑ::AbstractVector{<:Real}, Δt, μ; dir=1, tol=1e-12, maxit=200)
    rₛ, rₑ = norm(r̄ₛ), norm(r̄ₑ)

    # true anomaly change for the transfer
    Δθ = atan(norm(cross(r̄ₛ, r̄ₑ)), dot(r̄ₛ, r̄ₑ))             # angle between start and end positions
    if bound_angle(dir*angle_in_plane(r̄ₑ, r̄ₛ, MRP(1I))) > π
        Δθ = 2π - Δθ
    end
    # p iteration constants
    k = rₛ * rₑ * (1-cos(Δθ))
    L = rₛ + rₑ
    m = rₛ * rₑ * (1+cos(Δθ))

    # bounds
    pj = k / (L + √(2m))
    pjj = k / (L - √(2m))
    if Δθ > π
        pmin = 0.
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
            ΔE = bound_angle(atan(sinΔE, cosΔE))
            t = g + √(a^3/μ)*(ΔE-sinΔE)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) + √(a^3/μ)*(2k*sinΔE)/(p*(k-L*p))
        else        # Hyperbolic case
            ΔF = acosh(1 - rₛ/a * (1-f))
            t = g + √(-a^3/μ)*(sinh(ΔF)-ΔF)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) - √(-a^3/μ)*(2k*sinh(ΔF))/(p*(k-L*p))
        end
        err = abs(Δt - t)/Δt
        pnext = p + (Δt - t)/dtdp
        # If the next guess is outside of allowed bounds, use bisection
        if pnext < pmin
            pnext = (p + pmin)/2
        elseif pnext > pmax
            pnext = (p + pmax)/2
        end
    end
    # if it == maxit
    #     @show err
    #     println("Warning! Maximum number of iterations exceeded: Lambert solver failed to converge")
    # end
    v̄ₛ = (r̄ₑ - f*r̄ₛ)/g
    return v̄ₛ
end

function p_lambert(sorb::Orbit, eorb::Orbit, stime, etime; kwargs...) 
    r̄ₛ, r̄ₑ = state_vector(stime,sorb)[1], state_vector(etime,eorb)[1]
    d = sorb.i<π/2 ? 1 : -1     # choose short/long way based on clockwise sense of `sorb`
    v̄ₛ = p_lambert(r̄ₛ, r̄ₑ, etime-stime, sorb.primary.μ; dir=d, kwargs...)
    return Orbit(stime, r̄ₛ, v̄ₛ, sorb.primary)
end

