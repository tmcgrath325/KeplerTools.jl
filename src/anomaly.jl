"""
    M = time_to_mean(t, T, Mo=0., epoch=0.)

Computes the mean anomaly of an orbit with period `T` from a time `t`, the orbit's mean anomaly
at epoch `Mo`, and the epoch.
"""
function time_to_mean(t, T, Mo=0., epoch=0.) 
    return Mo + 2π*(t - epoch)/T
end

function time_to_mean(t, orb::Orbit)
    M = time_to_mean(t, orb.period, orb.Mo, orb.epoch)
    if orb.e < 1
        return wrap_angle(M)
    else
        return M
    end
end

"""
    t = mean_to_time(M, T, Mo=0, epoch=0, tmin=nothing)

Computes the time at which an orbit with period `T`, mean anomaly
at epoch `Mo`, and epoch `epoch` has the specified mean anomaly `M`.

Because this time is not uniquely defined for elliptical orbits, `tmin` can be specified
to ensure that the returned time falls between `tmin` and `tmin + T`.
"""
function mean_to_time(M, T, Mo=0, epoch=0, tmin=nothing) 
    t = epoch + (M-Mo)*T/(2π)
    if tmin === nothing
        return t
    else
        return t + T*ceil((tmin-t)/T)
    end
end

mean_to_time(M, orb::Orbit, tmin=nothing) = mean_to_time(M, orb.period, orb.Mo, orb.epoch, orb.e<1 ? tmin : nothing)

"""
    θ = mean_to_true(M, e)

Computes the true anomaly of an orbit with eccentricity `e` at the specified mean anomaly `M`.
"""
function mean_to_true(M, e)
    if e == 1       # Parabolic case
        throw(ArgumentError("Parabolic case (e=1) not implemented"))
    elseif e < 1    # Elliptical case
        E = kepler(M, e)
        return 2*atan(√(1+e)*sin(E/2),  √(1-e)*cos(E/2))
    else            # Hyperbolic case
        H = kepler(M, e)
        return 2*atan(√(e+1)*sinh(H/2), √(e-1)*cosh(H/2))
    end
end

"""
    θ = time_to_true(t, e, T, Mo=0., epoch=0.)

Computes the mean anomaly of an orbit with eccentricity `e` and period `T` from a time `t`, the orbit's mean anomaly
at epoch `Mo`, and the epoch.
"""
function time_to_true(t, e, T, Mo=0, epoch=0)
    M = time_to_mean(t, T, Mo, epoch)
    return mean_to_true(M, e)
end

time_to_true(t, orb::Orbit) = time_to_true(t, orb.e, orb.period, orb.Mo, orb.epoch)

"""
    M = true_to_mean(θ, e)

Computes the mean anomaly of an orbit with eccentricity `e` at the specified true anomaly `θ`.
"""
function true_to_mean(θ, e)
    if e == 1       # Parabolic case
        throw(ArgumentError("Parabolic case (e=1) not implemented"))
    elseif e < 1    # Elliptical case
        E = 2*atan(sin(θ/2)*√(1-e), cos(θ/2)*√(1+e))
        return wrap_angle(E - e*sin(E))
    else            # Hyperbolic case
        H = 2*atanh(tan(θ/2)*√((e-1)/(e+1)))
        # H = copysign(acosh((e+cos(θ))/(1+e*cos(θ))), θ)
        return e*sinh(H) - H
    end
end

true_to_mean(θ, orb::Orbit) = true_to_mean(θ, orb.e)

"""
    t = true_to_time(θ, e, T, Mo=0, epoch=0, tmin=nothing)

Computes the time at which an orbit with eccentricity `e`, period `T`, mean anomaly
at epoch `Mo`, and epoch `epoch` has the specified true anomaly `θ`.

Because this time is not uniquely defined for elliptical orbits, `tmin` can be specified
to ensure that the returned time falls between `tmin` and `tmin + T`.
"""
function true_to_time(θ, e, T, Mo=0, epoch=0, tmin=nothing)
    M = true_to_mean(θ, e)
    return mean_to_time(M, T, Mo, epoch, e<1 ? tmin : nothing)
end

# true_to_time(θ, e, a, μ, Mo=0, epoch=0, tmin=nothing) = true_to_time(θ, e, period(a, μ), Mo, epoch, tmin)
true_to_time(θ, orb::Orbit, tmin=nothing) = true_to_time(θ, orb.e, orb.period, orb.Mo, orb.epoch, tmin)
