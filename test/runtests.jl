using Base: byte_string_classify
using KeplerTools
using LinearAlgebra
using Test

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
    for i=1:50
        e = rand()
        θ = 2π*rand()
        M = KeplerTools.true_to_mean(θ, e)
        @test KeplerTools.mean_to_true(M, e) ≈ θ
    end 
    # hyperbolic case
    for i=1:50
        e = 1+9*rand()
        θ = π*rand() - π/2
        M = KeplerTools.true_to_mean(θ, e)
        @test KeplerTools.mean_to_true(M, e) ≈ θ
    end 
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
        svec = time_StateVector(t, bd.orbit)
        orb = Orbit(t, svec, bd.primary, bd.orbit.epoch)
        @test orb ≈ bd.orbit
        svec2 = time_StateVector(t, orb)
        @test isapprox(svec.position, svec2.position, atol=norm(svec.position)*1e-6)
        @test isapprox(svec.velocity, svec2.velocity, atol=norm(svec.velocity)*1e-6)
        orb2 = Orbit(t, svec2, bd.primary)
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
        stvec = time_StateVector(time, orb)
        @test Orbit(time, stvec, sun) ≈ orb
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
        stvec = time_StateVector(time, orb)
        @test Orbit(time, stvec, sun) ≈ orb
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
    v̄ = p_lambert(r̄ₛ, r̄ₑ, Δt, μ)
    @test isapprox(v̄, [0.000000193828, 0.000000101824, 0.00000000861759]; rtol=1e-3)

    # transfer from Kerbin to Duna, from alexmoon's app: https://alexmoon.github.io/ksp/
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    torb = KeplerTools.p_lambert(kerbin.orbit, duna.orbit, 5091522., 5588238. +5091522.)
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
        startstate, endstate = time_StateVector(starttime, orb), time_StateVector(endtime, orb)
        orb2 = p_lambert(orb, orb, starttime, endtime)
        @test orb2 ≈ orb
        @test isapprox(time_StateVector(starttime, orb2), startstate)
        @test isapprox(time_StateVector(endtime, orb2), endstate)
        d = orb.i<π/2 ? 1 : -1
        @test isapprox(p_lambert(startstate.position, endstate.position, endtime-starttime, μ; dir=d), startstate.velocity, rtol=1e-6)
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
        startstate, endstate = time_StateVector(starttime, orb), time_StateVector(endtime, orb)
        orb2 = p_lambert(orb, orb, starttime, endtime)
        @test orb2 ≈ orb
        @test isapprox(time_StateVector(starttime, orb2), startstate)
        @test isapprox(time_StateVector(endtime, orb2), endstate)
        d = orb.i<π/2 ? 1 : -1
        @test isapprox(p_lambert(startstate.position, endstate.position, endtime-starttime, μ; dir=d), startstate.velocity, rtol=1e-6)
    end
end

@testset "Departure and arrival trajectory calculation" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # define transfer trajectory
    starttime = 5091522.
    endtime = starttime + 5588238.
    torb = p_lambert(kerbin.orbit, duna.orbit, starttime, endtime)
    v̄rel = time_StateVector(starttime, torb).velocity - time_StateVector(starttime, kerbin.orbit).velocity

    # define parking orbit
    pkorb = Orbit(kerbin.eqradius+100000, kerbin)

    # get quick departure orbit (valid for circular parking orbits)
    dorb, Δv̄ = quick_departarrive_orbit(pkorb, v̄rel, starttime)
    @test time_StateVector(starttime, dorb).velocity ≈ v̄rel
    # TODO: match position of start orbit 
    # @test StateVector(dorb.epoch, dorb).position ≈ StateVector(dorb.epoch, pkorb).position

    # recover (something similar to) the quick departure orbit from the optimizer
    θₒ = time_to_true(dorb.epoch, pkorb)
    d2orb, eccobj, Δv̄2 = departarrive_orbit(pkorb, v̄rel, starttime)
    @test isapprox(d2orb.a, dorb.a; rtol=0.01, atol=0.01)
    @test isapprox(d2orb.e, dorb.e; rtol=0.01, atol=0.01)
    @test KT.isapproxangle(d2orb.i, dorb.i ; rtol=0.01, atol=0.01)
    @test KT.isapproxangle(d2orb.Ω, dorb.Ω; rtol=0.01, atol=0.01)
    @test KT.isapproxangle(d2orb.ω, dorb.ω; rtol=0.01, atol=0.01)

end

@testset "Patched Conic Transfer Calculations" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end
