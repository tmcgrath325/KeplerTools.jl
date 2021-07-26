using KeplerTools
using ProfileView

include(joinpath(@__DIR__, "..", "data", "kerbol_system.jl"))
stime1, stime2 = 0.,                    852. * 6. * 3600.
ftime1, ftime2 = 151. * 6. * 3600.,     453. * 6. * 3600.
startorb = Orbit(kerbin.eqradius + 100000, kerbin)
endorb   = Orbit(duna.eqradius   + 100000, duna)

ProfileView.@profview pc = Porkchop(startorb, endorb, stime1, stime2, ftime1, ftime2; npts=100)
ProfileView.@profview qpc = quickPorkchop(startorb, endorb, stime1, stime2, ftime1, ftime2; npts=100)
