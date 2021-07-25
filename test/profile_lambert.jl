using KeplerTools
using ProfileView

include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
rand_orbits = Orbit[]
rand_start_times = Float64[]
rand_end_times = Float64[]
for i=1:1000000
    orb = Orbit(moho.orbit.a + 3*duna.orbit.a*rand(),
                rand(),
                π*rand(),
                2π*rand(),
                2π*rand(),
                2π*rand(),
                1e6*rand(),
                sun
    )
    starttime = orb.epoch + 1e6*rand()
    endtime = starttime + orb.period*rand()
    push!(rand_orbits, orb)
    push!(rand_start_times, starttime)
    push!(rand_end_times, endtime)
end

function prof_p_lambert(rand_orbits::Vector{Orbit}, rand_start_times, rand_end_times, N=1000000)
    for i=1:N
        torb = p_lambert(rand_orbits[i], rand_orbits[i], rand_start_times[i], rand_end_times[i])
    end
end

ProfileView.@profview prof_p_lambert(rand_orbits, rand_start_times, rand_end_times)
