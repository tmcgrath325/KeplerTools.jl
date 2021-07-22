sun = Star("Kerbol", 261600000., 1.1723328e18);
kerbin =    CelestialBody("Kerbin", 600000.,  3.5316000e12,       21549.425,      0.,     Orbit(13599840256., 0., 0., 0., 0., 3.14, 0., sun));
moho =      CelestialBody("Moho",   250000.,  168.609378655e9,    1210000.,       190.,   Orbit(5263138304.,0.200000003, deg2rad(7.), deg2rad(70.), deg2rad(15.), 3.14, 0., sun));
eve =       CelestialBody("Eve",    700000.,  8171.730229211e9,   80500.,         0.,     Orbit(9832684544.,0.01, deg2rad(2.099999905), deg2rad(15.000000000), deg2rad(0.), 3.14, 0., sun));
duna =      CelestialBody("Duna",   320000.,  301.363211975e9,    65517.859375,   90.,    Orbit(20726155264., 0.051, deg2rad(0.06), deg2rad(135.5), 0., 3.14, 0., sun));