module KeplerTools

using LinearAlgebra: parent
using StaticArrays
using Rotations
using LinearAlgebra
using Optim
using PlotlyJS

export orbital_position, orbital_velocity, state_vector
export time_orbital_position, time_orbital_velocity, time_state_vector
export Orbit, StateVector, time_StateVector, OrbitalState
export CelestialBody, Star, SpaceObject
export period, period!, time_to_mean, mean_to_time, time_to_true, true_to_time
export basis_matrix, basis_matrix!
export angle_in_plane
export p_lambert
export departarrive_orbit, quick_departarrive_orbit
export draw_orbit, draw_central_body, draw_orbiting_body, draw_system
export Transfer

include("orbit.jl")
include("body.jl")
include("kepler.jl")
include("anomaly.jl")
include("basis.jl")
include("angles.jl")
include("utils.jl")
include("transfer/lambert.jl")
include("transfer/burn.jl")
include("transfer/departarrive.jl")
include("transfer/transfer.jl")
include("draw/draworbit.jl")
include("draw/drawobject.jl")
include("draw/drawsystem.jl")
include("draw/drawmarkers.jl")

end
