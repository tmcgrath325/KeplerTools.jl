using KeplerToolbox
using LinearAlgebra
using Test

@testset "Orbits and Bodies" begin
    sun = Star("Kerbol", 261600000., 1.1723328e18);
    @test length(sun.satellite_bodies) == 0
    kerbin =    CelestialBody("Kerbin", 1, 600000,  3.5316000e12,       21549.425,      0.,     Orbit(13599840256., 0., 0., 0., 0., 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 1
    moho =      CelestialBody("Moho",   4, 250000,  168.609378655e9,    1210000.,       190.,   Orbit(5263138304.,0.200000003, deg2rad(7.), deg2rad(70.), deg2rad(15.), 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 2
    eve =       CelestialBody("Eve",    5, 700000,  8171.730229211e9,   80500.,         0.,     Orbit(9832684544.,0.01, deg2rad(2.099999905), deg2rad(15.000000000), deg2rad(0.), 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 3
    duna =      CelestialBody("Duna",   6, 320000,  301.363211975e9,    65517.859375,   90.,    Orbit(20726155264., 0.051, deg2rad(0.06), deg2rad(135.5), 0., 3.14, 0., sun));
    @test length(sun.satellite_bodies) == 4
end

@testset "Anomalies" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    # period
    for bd in sun.satellite_bodies
        @test_throws KeyError bd.orbit.attrs[:period]
        @test period!(bd) == period(bd.orbit.a, bd.primary.μ)
        @test bd.orbit.attrs[:period] == period(bd.orbit)
    end
    # mean anomaly
    for bd in sun.satellite_bodies
        t = rand()*10*period!(bd)  - 5*period!(bd) + bd.orbit.epoch
        M = time_to_mean(t, bd.orbit)
        @test isapprox(t, mean_to_time(M, bd.orbit, t-1), atol=t*1e-12)
    end
    # true anomaly
    for bd in sun.satellite_bodies
        t = rand()*10*period!(bd)  - 5*period!(bd) + bd.orbit.epoch
        θ = time_to_true(t, bd.orbit)
        @test isapprox(t, true_to_time(θ, bd.orbit, t-1), atol=t*1e-12)
    end
end

@testset "Bases" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end

@testset "State vectors" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    for bd in sun.satellite_bodies
        t = rand()*10*period!(bd)  - 5*period!(bd) + bd.orbit.epoch
        svec = StateVector(t, bd.orbit)
        orb = Orbit(t, svec, bd.primary, bd.orbit.epoch)
        @test orb ≈ bd.orbit
        svec2 = StateVector(t, orb)
        @test isapprox(svec.position, svec2.position, atol=norm(svec.position)*1e-6)
        @test isapprox(svec.velocity, svec2.velocity, atol=norm(svec.velocity)*1e-6)
        orb2 = Orbit(t, svec2.position, svec2.velocity, bd.primary)
        @test orb2 ≈ bd.orbit
    end
end

@testset "Plane angles" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end

@testset "Lambert Solvers" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end

@testset "Patched Conic Transfer Calculations" begin
    include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
    #TODO
end
