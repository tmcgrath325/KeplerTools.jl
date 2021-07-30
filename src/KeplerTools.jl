module KeplerTools

using StaticArrays
using Rotations
using LinearAlgebra
using Optim
using PlotlyJS

export Orbit
export CelestialBody, Star, SpaceObject
export period, orbital_angle, orbital_distance, flight_path_angle, orbital_direction
export orbital_position, orbital_velocity, state_vector
export time_orbital_position, time_orbital_velocity, time_state_vector
export kepler, time_to_mean, mean_to_time, time_to_true, true_to_time
export basis_rotation, orientation_rotation
export angle_in_plane

export p_lambert, plane_change_p_lambert
export departure_orbit, arrival_orbit, fast_departure_orbit, fast_arrival_orbit
export ejection_angle, ejection_time, ejection_position, ejection_velocity, ejection_state_vector
export insertion_angle, insertion_time, insertion_position, insertion_velocity, insertion_state_vector
export propagate_across_soi
export Burn, orientation_components, apply_burn
export Transfer, fastTransfer, match_transfer_patch_times, match_patch_positions
export Porkchop, fastPorkchop

export draw_orbit, draw_central_body, draw_orbiting_body, draw_system
export draw_burn_arrow, draw_orbit_position, draw_angle_in_plane
export draw_transfer, draw_ejection, draw_insertion
export draw_porkchop

include("orbit.jl")
include("body.jl")
include("kepler.jl")
include("anomaly.jl")
include("basis.jl")
include("angles.jl")
include("utils.jl")

include("transfer/lambert.jl")
include("transfer/departarrive.jl")
include("transfer/soipropagate.jl")
include("transfer/burn.jl")
include("transfer/transfer.jl")
include("transfer/porkchop.jl")

include("draw/draworbit.jl")
include("draw/drawobject.jl")
include("draw/drawsystem.jl")
include("draw/drawmarkers.jl")
include("draw/drawtransfer.jl")
include("draw/drawporkchop.jl")

end
