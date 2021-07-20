function p_lambert(r̄ₛ::SVector{Real,3}, r̄ₑ::SVector{Real,3}, t::Real, μ::Real, tol=1e-6::Real, maxit=200::Int)
    rₛ, rₑ = norm(r̄ₛ), norm(rₑ)
    # true anomaly change for the transfer
    Δθ = atan2(cross(r̄ₛ, r̄ₑ),dot(r̄ₛ,r̄ₑ))
    phase = angle_in_plane(sorb, r̄ₛ, r̄ₑ)
    if phase > π
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
            Δt = g + √(a^3/μ)*(ΔE-sinΔE)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) + √(a^3/μ)*(2k*sinΔE)/(p*(k-L*p))
        else        # Hyperbolic case
            ΔF = acosh(1 - rₛ/a * (1-f))
            Δt = g + √(-a^3/μ)*(sinh(ΔF)-ΔF)
            dtdp = -g/(2p) - 1.5a*(t-g)*(k^2 + (2m-L^2)*p^2)/(m*k*p^2) - √(-a^3/μ)*(2k*sinh(ΔF))/(p*(k-L*p))
        end
        err = abs(Δt - t)
        pnext = p + t/dtdp
    end
    if it == maxit
        error("Maximum number of iterations reached")
    end
    v̄ = (r̄ₑ - f*r̄ₛ)/g
    return v̄
end

function p_lambert(sorb::Orbit, eorb::Orbit, stime, etime, tol=1e-6, maxit=200) 
    r̄ₛ, r̄ₑ = state_vector(sorb,stime), state_vector(eorb,etime)
    v̄ = p_lambert(r̄ₛ, r̄ₑ, etime-stime, primary(sorb).μ, tol, maxit)
    return Orbit(stime, r̄ₛ, v̄, primary(sorb))
end
