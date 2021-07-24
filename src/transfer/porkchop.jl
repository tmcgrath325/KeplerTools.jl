struct Porkchop
    startorbit::Orbit
    endorbit::Orbit
    starttimes::Vector{Float64}
    flighttimes::Vector{Float64}
    transferorbits::Matrix{Orbit}
end

function Porkchop(startorbit::Orbit, endorbit::Orbit, starttimes, flighttimes)
    transferorbits = fill(Orbit)
end
