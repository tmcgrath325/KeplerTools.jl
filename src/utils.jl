function hash_struct(struc, h::UInt, fields=fieldnames(typeof(struc)))
    for f in fields
        h = Base.hash(getfield(struc,f), h)
    end
    return h
end

const EarthTime =  Dict("year" =>   365,    # days per year
                        "day"  =>   24 ,    # hours per day
                        "hour" =>   60 ,    # minutes per hour
                        "minute" => 60 ,    # seconds per minute
)

const KerbalTime = Dict("year" =>   424,    # days per year
                        "day"  =>   6,      # hours per day
                        "hour" =>   60,     # minutes per hour
                        "minute" => 60,     # seconds per minute
)

function seconds_to_datetime(nseconds, timedefs=KerbalTime)
    year, day, hour, minute = timedefs["year"], timedefs["day"], timedefs["hour"], timedefs["minute"]

    y = floor(nseconds/(minute*hour*day*year))+1       # starting from year 1
    nseconds -= minute*hour*day*year*(y-1)
    d = floor(nseconds/(minute*hour*day))+1            # starting from day 1
    nseconds -= minute*hour*day*(d-1)
    h = floor(nseconds/(minute*hour))                  # starting from hour 0
    nseconds -= minute*hour*(h)
    m = floor(nseconds/minute)                         # starting from minute 0
    nseconds -= minute*(m)
    s = round(nseconds)
    return NTuple{5,Int}((y,d,h,m,s))
end

function seconds_to_datetime_string(nseconds, timedefs=KerbalTime)
    ny, nd, nh, nm, ns = seconds_to_datetime(nseconds, timedefs)
    return "Year $ny, Day $nd, $(lpad(nh, 2, '0')):$(lpad(nm, 2, '0')):$(lpad(ns, 2, '0'))"
end

function wrap_angle(θ::Real, min::Real = 0)
    if min <= θ < (min + 2π)
        return θ
    else
        return θ - 2π*floor((θ-min)/2π)
    end
end

function angledist(θ₁, θ₂)
    dist = wrap_angle(θ₁ - θ₂)
    return dist > π ? dist - 2π : dist
end

function isapproxangle(θ₁,θ₂; atol=1e-12, rtol=1e-12, kwargs...)
    return isapprox(angledist(θ₁, θ₂),0;atol=atol,rtol=rtol,kwargs...) || isapprox(rem,1;atol=atol,rtol=rtol,kwargs...)
end

function wrap_acos(x::Real)
    if abs(x) > 1
        x = copysign(1,x)
    end
    return acos(x)
end

function wrap_asin(x::Real)
    if abs(x) > 1
        x = copysign(1,x)
    end
    return asin(x)
end

