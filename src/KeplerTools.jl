module KeplerTools

using StaticArrays
using Rotations
using LinearAlgebra
using Optim

export Orbit, StateVector, OrbitalState, state_vector
export CelestialBody, Star, SpaceObject
export period, period!, time_to_mean, mean_to_time, time_to_true, true_to_time
export basis_matrix, basis_matrix!
export angle_in_plane
export lambert
export departarrive_orbit

include("orbit.jl")
include("body.jl")
include("kepler.jl")
include("anomaly.jl")
include("basis.jl")
include("angles.jl")
include("utils.jl")
include("transfer/lambert.jl")
include("transfer/departarrive.jl")

end
