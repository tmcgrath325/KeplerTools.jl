using KeplerTools
using LinearAlgebra
using Test
using Statistics

const KT = KeplerTools

@testset "Orbits and Bodies" begin
    sun = Star("Kerbol", 261600000., 1.1723328e18);
    @test length(sun.satellite_bodies) == 0
    kerbin =    CelestialBody("Kerbin", 600000.,  3.5316000e12,       21549.425,      0.,     Orbit(13599840256., 0., 0., 0., 0., 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 1
    moho =      CelestialBody("Moho",   250000.,  168.609378655e9,    1210000.,       190.,   Orbit(5263138304.,0.200000003, deg2rad(7.), deg2rad(70.), deg2rad(15.), 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 2
    eve =       CelestialBody("Eve",    700000.,  8171.730229211e9,   80500.,         0.,     Orbit(9832684544.,0.01, deg2rad(2.099999905), deg2rad(15.000000000), deg2rad(0.), 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 3
    duna =      CelestialBody("Duna",   320000.,  301.363211975e9,    65517.859375,   90.,    Orbit(20726155264., 0.051, deg2rad(0.06), deg2rad(135.5), 0., 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 4
end

@testset "Kepler's Equation Solver" begin
    # elliptical case
    es = collect(range(0,stop=0.999,length=1000))
    Ms = collect(range(-π,stop=π,length=1000))
    errs = fill(NaN, length(es), length(Ms))
    for (i,e) in enumerate(es) for (j,M) in enumerate(Ms)
        θ = KT.mean_to_true(M, e)
        errs[i,j] = abs(KT.angledist(KT.true_to_mean(θ, e), M))
    end end
    @test maximum(errs) < 1e-9
    @test mean(errs) < 1e-12
    # hyperbolic case
    es = collect(range(1.01,stop=100,length=1000))
    Ms = collect(range(-100*π,stop=100*π,length=1000))
    errs = fill(NaN, length(es), length(Ms))
    for (i,e) in enumerate(es) for (j,M) in enumerate(Ms)
        θ = KT.mean_to_true(M, e)
        errs[i,j] = abs(KT.angledist(KT.true_to_mean(θ, e), M))
    end end
    @test maximum(errs) < 1e-9
    @test mean(errs) < 1e-12
end

@testset "Anomalies" begin
    # elliptical orbits from Kerbol system
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # mean anomaly
    for bd in sun.satellite_bodies
        t = rand()*10*bd.orbit.period  - 5*bd.orbit.period + bd.orbit.epoch
        M = time_to_mean(t, bd.orbit)
        @test isapprox(t, mean_to_time(M, bd.orbit, t-1), rtol=1e-12)
    end
    # true anomaly
    for bd in sun.satellite_bodies
        t = rand()*10*bd.orbit.period  - 5*bd.orbit.period + bd.orbit.epoch
        θ = time_to_true(t, bd.orbit)
        @test isapprox(t, true_to_time(θ, bd.orbit, t-1), rtol=1e-12)
    end

    # random elliptical orbits
    for i in 1:50
        orb = Orbit(moho.orbit.a + 3*duna.orbit.a*rand(),
                    rand(),
                    π*rand(),
                    2π*rand(),
                    2π*rand(),
                    2π*rand(),
                    1e6*rand(),
                    sun
                    )
        t = rand()*10*orb.period  - 5*orb.period + orb.epoch
        M = time_to_mean(t, orb)
        @test isapprox(t, mean_to_time(M, orb, t-1), rtol=1e-12)
        θ = time_to_true(t, orb)
        @test isapprox(t, true_to_time(θ, orb, t-1), rtol=1e-12)
    end

    # random hyperbolic orbits
    for i in 1:50
        e = 1+4*rand()
        a = -(moho.orbit.a + 3*duna.orbit.a*rand())/(e-1)
        orb = Orbit(a,
                    e,
                    π*rand(),
                    2π*rand(),
                    2π*rand(),
                    2π*rand(),
                    1e6*rand(),
                    sun
                    )
        t = rand()*10*orb.period  - 5*orb.period + orb.epoch
        M = time_to_mean(t, orb)
        @test isapprox(t, mean_to_time(M, orb, t-1), rtol=1e-12)
        θ = time_to_true(t, orb)
        @test isapprox(t, true_to_time(θ, orb, t-1), rtol=1e-12)
    end

end

@testset "Bases and Rotations" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end

@testset "State vectors" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # recover orbits from state vectors
    for bd in sun.satellite_bodies
        t = (rand()-2)*10*bd.orbit.period + bd.orbit.epoch
        pos, vel = time_state_vector(t, bd.orbit)
        orb = Orbit(t, pos, vel, bd.primary, bd.orbit.epoch)
        @test orb ≈ bd.orbit
        pos2, vel2 = time_state_vector(t, orb)
        @test isapprox(pos, pos2, atol=norm(pos)*1e-6)
        @test isapprox(vel, vel2, atol=norm(vel)*1e-6)
        orb2 = Orbit(t, pos2, vel2, bd.primary)
        @test orb2 ≈ bd.orbit
    end
    # generate random elliptical orbits for testing, and try to recover them from state vectors
    for i in 1:50
        orb = Orbit(moho.orbit.a + 3*duna.orbit.a*rand(),
                    rand(),
                    π*rand(),
                    2π*rand(),
                    2π*rand(),
                    2π*rand(),
                    1e6*rand(),
                    sun
                    )
        time = orb.epoch + 1e6*rand()
        pos, vel = time_state_vector(time, orb)
        @test Orbit(time, pos, vel, sun) ≈ orb
    end

    # generate random hyperbolic orbits for testing, and try to recover them from state vectors
    for i in 1:50
        e = 1+4*rand()
        a = -(moho.orbit.a + 3*duna.orbit.a*rand())/(e-1)
        orb = Orbit(a,
                    e,
                    π*rand(),
                    2π*rand(),
                    2π*rand(),
                    2π*rand(),
                    1e6*rand(),
                    sun
                    )
        time = orb.epoch + 1e6*rand()
        pos, vel = time_state_vector(time, orb)
        @test Orbit(time, pos, vel, sun) ≈ orb
    end
end

@testset "Plane angles" begin
    @test angle_in_plane([1,0,0], [1,0,0]) == 0.
    @test angle_in_plane([0,1,0], [1,0,0]) == π/2
    @test angle_in_plane([1,0,0], [0,1,0]) == 3π/2
    @test angle_in_plane([1,1,0], [1,0,0]) == π/4
    @test angle_in_plane([1,0,0], [1,1,0]) == 7π/4


    
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # phase angle between Kerbin and Duna for a typical transfer, from from alexmoon's app: https://alexmoon.github.io/ksp/
    @test isapprox(angle_in_plane(duna.orbit, kerbin.orbit, 5091552), deg2rad(36.49); rtol=0.001)

end

@testset "Lambert Solver" begin
    # test from Problems 5.3 and 5.4, http://www.braeunig.us/
    #  au = 1.495978707e11     # number of meters per astronomical unit
    μ = 3.964016e-14        # the real sun AU^3/m^2
    Δt = 207*24*60*60
    r̄ₛ = [0.473265, -0.899215, 0.] # .* au
    r̄ₑ = [0.066842, 1.561256, 0.030948] # .* au
    v̄ = KT.p_lambert_velocity(r̄ₛ, r̄ₑ, Δt, μ)
    @test isapprox(v̄, [0.000000193828, 0.000000101824, 0.00000000861759]; rtol=1e-3)

    # transfer from Kerbin to Duna, from alexmoon's app: https://alexmoon.github.io/ksp/
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    torb = KT.p_lambert(kerbin.orbit, duna.orbit, 5091522., 5588238. +5091522.)
    @test isapprox(torb.a, 1.68e10, rtol=0.01)            
    @test isapprox(torb.e, 0.194, rtol=0.01)

    # generate random elliptical orbits for testing, and recover them with the Lambert solver
    μ = sun.μ
    for i in 1:50
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
        (spos, svel), (epos, evel) = time_state_vector(starttime, orb), time_state_vector(endtime, orb)
        orb2 = p_lambert(orb, orb, starttime, endtime)
        @test orb2 ≈ orb
        @test isapprox(time_orbital_position(starttime, orb2), spos)
        @test isapprox(time_orbital_velocity(starttime, orb2), svel)
        @test isapprox(time_orbital_position(endtime, orb2), epos)
        @test isapprox(time_orbital_velocity(endtime, orb2), evel)
        d = orb.i<π/2 ? 1 : -1
        @test isapprox(KT.p_lambert_velocity(spos, epos, endtime-starttime, μ; dir=d), svel, rtol=1e-6)
    end

    # generate random hyperbolic orbits for testing, and recover them with the Lambert solver
    μ = sun.μ
    for i in 1:50
        e = 1+4*rand()
        a = -(moho.orbit.a + 3*duna.orbit.a*rand())/(e-1)
        orb = Orbit(a,
                    e,
                    π*rand(),
                    2π*rand(),
                    2π*rand(),
                    2π*rand(),
                    1e6*rand(),
                    sun
                    )
        starttime = orb.epoch + 1e6*rand()
        endtime = starttime + orb.period*rand()
        (spos, svel), (epos, evel) = time_state_vector(starttime, orb), time_state_vector(endtime, orb)
        orb2 = p_lambert(orb, orb, starttime, endtime)
        @test orb2 ≈ orb
        @test isapprox(time_orbital_position(starttime, orb2), spos)
        @test isapprox(time_orbital_velocity(starttime, orb2), svel)
        @test isapprox(time_orbital_position(endtime, orb2), epos)
        @test isapprox(time_orbital_velocity(endtime, orb2), evel)
        d = orb.i<π/2 ? 1 : -1
        @test isapprox(KT.p_lambert_velocity(spos, epos, endtime-starttime, μ; dir=d), svel, rtol=1e-6)
    end
end

@testset "Departure and arrival trajectory calculation" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # define transfer trajectory
    starttime = 5091522.
    endtime = starttime + 5588238.
    torb = p_lambert(kerbin.orbit, duna.orbit, starttime, endtime)
    v̄rel = time_orbital_velocity(starttime, torb) - time_orbital_velocity(starttime, kerbin.orbit)

    # define parking orbit
    pkorb = Orbit(kerbin.eqradius+100000, kerbin)

    # get fast departure orbit (valid for circular parking orbits)
    dorb, Δv̄ = fast_departure_orbit(pkorb, v̄rel, starttime)
    @test isapprox(time_orbital_velocity(KT.ejection_time(dorb), dorb), v̄rel; rtol=1e-6)

    # recover (something similar to) the fast departure orbit from the optimizer
    θₒ = time_to_true(dorb.epoch, pkorb)
    d2orb, eccobj, Δv̄2 = departure_orbit(pkorb, v̄rel, starttime)
    @test isapprox(d2orb, dorb; atol=0.01, rtol=0.01)
end

@testset "Instantaneous Burns" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end

@testset "Patched Conic Transfer Calculations" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # Transfer from Kerbin to Duna
    startorb = Orbit(kerbin.eqradius+100000, kerbin)
    endorb   = Orbit(duna.eqradius  +100000, duna)
    starttime = 5091522.
    endtime = starttime + 5588238.;
    tfer = Transfer(startorb, endorb, starttime, endtime)
    ftfer = fastTransfer(startorb, endorb, starttime, endtime)
    @test isapprox(tfer.Δv, ftfer.Δv, rtol=0.01)

    # Porkchop plot containing time ranges around above transfer
    startrange = [0.,               852*6*3600.]
    flightrange =[151 *6 *3600.,    453*6. *3600.]
    pc = Porkchop(startorb, endorb, startrange..., flightrange...)
    @test isapprox(minimum(pc.Δv), tfer.Δv; rtol=0.01)
    @test Transfer(startorb, endorb, KT.best_start_end_times(pc)...).Δv == minimum(pc.Δv)

    fpc = fastPorkchop(startorb, endorb, startrange..., flightrange...)
    @test isapprox(minimum(fpc.Δv), tfer.Δv; rtol=0.01)
    @test fastTransfer(startorb, endorb, KT.best_start_end_times(fpc)...).Δv == minimum(fpc.Δv)

    # optimize patch times and positions
    @test sum(norm, KT.patch_position_errors(tfer)) > 0
    tfer = match_patch_positions(tfer)

    @test sum(norm, KT.patch_position_errors(tfer)) < 1

end

@testset "SoI patch propagation" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # TODO
end


