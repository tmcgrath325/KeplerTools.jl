struct Porkchop
    startorbit::Orbit
    endorbit::Orbit
    starttimes::Vector{Float64}
    flighttimes::Vector{Float64}
    # transfers::Matrix{Transfer}
    Δv::Matrix{Float64}
end

function Porkchop(startorb::Orbit, endorb::Orbit, starttimes::AbstractVector{<:Real}, flighttimes::AbstractVector{<:Real})
    Δv = fill(NaN, length(starttimes), length(flighttimes))
    for (i,st) in enumerate(starttimes)
        for (j,ft) in enumerate(flighttimes)
            Δv[i,j] = Transfer(startorb, endorb, st, st+ft).Δv
        end
    end
    return Porkchop(startorb, endorb, starttimes, flighttimes, Δv)
end

Porkchop(startorb::Orbit, endorb::Orbit, stime1, stime2, ftime1, ftime2; npts=100
    ) = Porkchop(startorb, endorb, collect(range(stime1, stop=stime2, length=npts)), collect(range(ftime1, stop=ftime2, length=npts)))