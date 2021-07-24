const SolarTime =  Dict("year" =>   365,    # days per year
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

function wrap_angle(θ::Real, min::Real = 0)
    if min <= θ < (min + 2π)
        return θ
    else
        return θ - 2π*floor((θ-min)/2π)
    end
end

function isapproxangle(θ₁,θ₂; atol=1e-12, rtol=1e-12, kwargs...)
    diff = θ₁ - θ₂
    rem = (diff%2π)/(2π)
    return isapprox(rem,0;atol=atol,rtol=rtol,kwargs...) || isapprox(rem,1;atol=atol,rtol=rtol,kwargs...)
end

function wrap_acos(x::Real)
    if abs(x) > 1
        x = copysign(1,x)
    end
    return acos(x)
end
