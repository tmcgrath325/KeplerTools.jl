"""
    v̄ = p_lambert_velocity(r̄ₛ, r̄ₑ, Δt, μ; dir=1)

Computes the single-revolution orbit around a body with gravity parameter `μ`, passing through positions `r̄ₛ` and `r̄ₑ`, 
and separated by time interval `Δt`.

Because there are two solutions, `dir` can be used to specifiy a clockwise or anti-clockwise orbit
relative to the coordinate system.
"""
function p_lambert_velocity(r̄ₛ::AbstractVector{<:Real}, r̄ₑ::AbstractVector{<:Real}, Δt, μ; dir=1, tol=1e-12, maxit=200)
    rₛ, rₑ = norm(r̄ₛ), norm(r̄ₑ)

    # true anomaly change for the transfer
    Δθ = atan(norm(cross(r̄ₛ, r̄ₑ)), dot(r̄ₛ, r̄ₑ))             # angle between start and end positions
    if wrap_angle(dir*angle_in_plane(r̄ₑ, r̄ₛ, MRP(1I))) > π
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
            ΔE = wrap_angle(atan(sinΔE, cosΔE))
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
    #     @warn "Lambert solver failed to converge: err = $err"
    # end
    v̄ₛ = (r̄ₑ - f*r̄ₛ)/g
    return v̄ₛ
end

"""
    torb = p_lambert_velocity(sorbₛ, eorb, stime, etime, μ, s_patch_position=@SVector(zeros(3)), e_patch_position=@SVector(zeros(3)); kwargs...)

Computes the single-revolution orbit around a body with gravity parameter `μ`, passing through the positions of orbits 
`sorb` and `eorb` at times `stime` and `etime`, respectively.

Positions `s_patch_position` and `e_patch_position` can be provided to modify the start and end positions, respectively.
This can be used to minimize patching error across spheres of influence.
"""
function p_lambert(sorb::Orbit, eorb::Orbit, stime, etime, s_patch_position=@SVector(zeros(3)), e_patch_position=@SVector(zeros(3)); kwargs...) 
    r̄ₛ, r̄ₑ = time_orbital_position(stime,sorb) + s_patch_position, time_orbital_position(etime,eorb) + e_patch_position
    d = sorb.i<π/2 ? 1 : -1     # choose short/long way based on clockwise sense of `sorb`
    v̄ₛ = p_lambert_velocity(r̄ₛ, r̄ₑ, etime-stime, sorb.primary.μ; dir=d, kwargs...)
    return Orbit(stime, r̄ₛ, v̄ₛ, sorb.primary)
end

function plane_change_p_lambert(sorb::Orbit, eorb::Orbit, stime, etime, s_patch_position=@SVector(zeros(3)), e_patch_position=@SVector(zeros(3)); kwargs...)
    r̄ₛ, r̄ₑ = time_orbital_position(stime,sorb) + s_patch_position, time_orbital_position(etime,eorb) + e_patch_position
    
    # change r̄ₑ so that it is in-plane with the starting orbit, at the same angle and distance.
    r̄ₑplane = inertial_to_perifocal_bases(r̄ₑ, sorb)
    r̄ₑplane = perifocal_to_inertial_bases(norm(r̄ₑ)*normalize(@SVector([r̄ₑplane[1], r̄ₑplane[2], 0])), sorb)
    
    d = sorb.i<π/2 ? 1 : -1     # choose short/long way based on clockwise sense of `sorb`
    v̄ₛ = p_lambert_velocity(r̄ₛ, r̄ₑplane, etime-stime, sorb.primary.μ; dir=d, kwargs...)
    torb1 = Orbit(stime, r̄ₛ, v̄ₛ, sorb.primary)
    
    # identify the true anomaly at plane change
    θs = time_to_true(stime,torb1)
    pcangle = min(wrap_angle(time_to_true(etime,torb1) - θs - π/2, 0.))
    θpc = θs + pcangle
    pctime = true_to_time(θpc, torb1, stime)

    # rotate velocity vector prior to plane change burn to get post-burn velocity
    r̄pc, v̄pc1 = state_vector(θpc, torb1)
    n̂1 = normalize(cross(r̄pc, v̄pc1))
    n̂2 = normalize(cross(r̄pc, r̄ₑ))
    Ri = align_vectors(n̂1, n̂2)
    v̄pc2 = Ri * v̄pc1
    torb2 = Orbit(pctime, r̄pc, v̄pc2, torb1.primary)

    return torb1, torb2
end