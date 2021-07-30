# KeplerTools.jl

Tools for 2-body orbit propagation, Lambert problem solving, and ejection/insertion trajectory optimization. Put together as a project to get comfortable with [Julia](https://julialang.org/). Intended for use to plan simple missions in Kerbal Space Program. 

## Orbit Propagation

Get the position, velocity, or both (i.e. a state vector) for an orbit as a function of either true anomaly or time:

```julia
julia> using KeplerTools

julia> include(joinpath(dirname(pathof(KeplerTools)), "..", "data", "kerbol_system.jl"));

julia> position, velocity = state_vector(3.14, kerbin.orbit)
([-1.3599823007697136e10, 2.1659825247473374e7, 0.0], [-14.786987156243779, -9284.488948751467, 0.0])

julia> position, velocity = time_state_vector(0, kerbin.orbit)
([-1.3599823007697136e10, 2.1659825247473374e7, 0.0], [-14.786987156243779, -9284.488948751467, 0.0])
```

## Lambert Solver

Solve Lambert's problem given two orbits, a start time, and an end time.

```julia
julia> starttime = 5091522.; 

julia> endtime = starttime + 5588238.;

julia> torb = p_lambert(kerbin.orbit, duna.orbit, starttime, endtime)
Orbit above Kerbol:
  semi-major axis:              1.6833240192579033e10 m
  eccentricity:                 0.1926164653078007
  inclination:                  0.1434866354451145°
  RA of ascending node:         19.065485668978525°
  argument of periapsis:        354.82322931318316°
  mean anomaly at epoch:        3.4404814375472044°
  epoch:                        5.091522e6 s
```

## Ejection/Insertion Optimization

Obtain an optimal ejection or insertion trajectory given a parking orbit and a relative velocity at escape from/encounter with the sphere of influence. 

```julia
julia> parking_orb = Orbit(kerbin.eqradius + 100000, kerbin)
Orbit above Kerbin:
  semi-major axis:              700000.0 m
  eccentricity:                 0.0
  inclination:                  0.0°
  RA of ascending node:         0.0°
  argument of periapsis:        0.0°
  mean anomaly at epoch:        0.0°
  epoch:                        0.0 s

julia> v̄_ejection = time_orbital_velocity(starttime, torb) - time_orbital_velocity(starttime, kerbin.orbit)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 -138.4274524214611
  853.0248945699786
   25.38365880099945

julia> ejection_orb = departure_orbit(parking_orb, v̄_ejection, starttime)[1]
Orbit above Kerbin:
  semi-major axis:              -5.322430914063159e6 m
  eccentricity:                 1.1315005119137438
  inclination:                  3.4432174905445865°
  RA of ascending node:         308.4387644550581°
  argument of periapsis:        358.686708883049°
  mean anomaly at epoch:        0.04289926269835575°
  epoch:                        5.004695828269423e6 s

julia> target_orb = Orbit(duna.eqradius + 100000, duna)
Orbit above Duna:
  semi-major axis:              420000.0 m
  eccentricity:                 0.0
  inclination:                  0.0°
  RA of ascending node:         0.0°
  argument of periapsis:        0.0°
  mean anomaly at epoch:        0.0°
  epoch:                        0.0 s

julia> v̄_insertion = time_orbital_velocity(endtime, torb) - time_orbital_velocity(endtime, duna.orbit)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 -97.7579049517824
 890.6147177047887
 -22.10024628653403

julia> insertion_orb = arrival_orbit(target_orb, v̄_insertion, endtime)[1]
Orbit above Duna:
  semi-major axis:              -381152.85373850045 m
  eccentricity:                 2.1019197583378046
  inclination:                  1.605586783222423°
  RA of ascending node:         214.6218218467353°
  argument of periapsis:        180.05673308987429°
  mean anomaly at epoch:        -0.037260288183548824°
  epoch:                        1.0732725844524838e7 s
```

## Transfers between two orbits

Obtain all ejection, insertion, and transfer orbits for a ballistic transfer between any two orbits within a star system:

```julia
julia> tfer = Transfer(parking_orb, target_orb, starttime, endtime)
Transfer from orbit around Kerbin to orbit around Duna
  departure time:       5.091522e6 s
  arrival time:         1.067976e7 s
  Δv:                   1692.0562845940308 m/s
  
julia> tfer.transfer_orbits
1-element Vector{Orbit}:
 Orbit above Kerbol:
  semi-major axis:              1.6833240192579033e10 m
  eccentricity:                 0.1926164653078007
  inclination:                  0.1434866354451145°
  RA of ascending node:         19.065485668978525°
  argument of periapsis:        354.82322931318316°
  mean anomaly at epoch:        3.4404814375472044°
  epoch:                        5.091522e6 s
  
julia> tfer.ejection_orbits
1-element Vector{Orbit}:
 Orbit above Kerbin:
  semi-major axis:              -5.322430914063159e6 m
  eccentricity:                 1.1315005119137438
  inclination:                  3.4432174905445865°
  RA of ascending node:         308.4387644550581°
  argument of periapsis:        358.686708883049°
  mean anomaly at epoch:        0.04289926269835575°
  epoch:                        5.004695828269423e6 s
  
julia> tfer.insertion_orbits
1-element Vector{Orbit}:
 Orbit above Duna:
  semi-major axis:              -381152.85373850045 m
  eccentricity:                 2.1019197583378046
  inclination:                  1.605586783222423°
  RA of ascending node:         214.6218218467353°
  argument of periapsis:        180.05673308987429°
  mean anomaly at epoch:        -0.037260288183548824°
  epoch:   
  
julia> tfer.burns
2-element Vector{Burn}:
 Instantaneous burn in orbit above Kerbin:
  UT:           5.004695828269423e6 s
  total Δv:     1046.4726500752242 m/s
  prograde:     1026.9996921227305 m/s
  normal:       196.93801584416403 m/s
  radial:       39.90059712139015 m/s


 Instantaneous burn in orbit above Duna:
  UT:           1.0732725844524838e7 s
  total Δv:     645.5836345188064 m/s
  prograde:     -644.2281248462269 m/s
  normal:       41.80131947390173 m/s
  radial:       1.001002430789867 m/s
```
