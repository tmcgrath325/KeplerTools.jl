"""
    basis = basis_rotation(Ω, i, ω)
    basis = basis_rotation(xvec, yvec)
    basis_rotation() = one(RotMatrix{3, Float64})

From the RAAN `Ω`, inclination `i`, and argument of the periapsis `ω`, computes the rotation matrix that transforms 
an orbit's perifocal plane to the inertial frame. For this matrix, the first column is a unit vector that points toward 
the periapsis of the orbit, and the third column is a unit vector that points in the normal direction. 

If two vectors `xvec` and `yvec` are provided, a rotation matrix is constructed that describes the plane the two vectors
define, where `normalize(xvec)` is the first basis vector.
"""
function basis_rotation(Ω, i, ω)
    # Ω, ω, i = Float64(Ω), Float64(ω), Float64(i)
    # R₁ = [cos(Ω) -sin(Ω)  0.;    # rotate around z-axis to match Ω
    #       sin(Ω)  cos(Ω)  0.;
    #       0.      0.      1.]
    # R₂ = [1.      0.      0.;    # rotate around x-axis to match i
    #       0.  cos(i) -sin(i);
    #       0.  sin(i)  cos(i)]
    # R₃ = [cos(ω) -sin(ω)  0.;    # rotate around z-axis to match ω
    #       sin(ω)  cos(ω)  0.;
    #       0.      0.      1.]
    # return R₃*R₂*R₁
    return RotZXZ(Ω, i, ω)
end

function basis_rotation(xvec::AbstractVector{<:Real}, yvec::AbstractVector{<:Real})
    @assert length(xvec) == length(yvec) == 3
    xbasisvec = normalize(xvec)
    ybasisvec = normalize(yvec - xbasisvec .* dot(xbasisvec, yvec))
    zbasisvec = cross(xbasisvec, ybasisvec)
    return RotMatrix(hcat(xbasisvec,ybasisvec,zbasisvec))
end

basis_rotation() = one(RotMatrix{3, Float64})

# methods for switching orbital reference frames

"""
    pvec = inertial_to_perifocal_bases(vec, basis)
    pvec = inertial_to_perifocal_bases(vec, orb)

Transforms a vector `vec` from the inertial frame to an orbit's perifocal frame,
defined by `basis`.
"""
inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, basis::Rotation{3,Float64}) = basis\vec
inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, orb::Orbit) = inertial_to_perifocal_bases(vec, orb.basis)

"""
    ivec = perifocal_to_inertial_bases(vec, basis)
    ivec = perifocal_to_inertial_bases(vec, orb)

Transforms a vector `vec` from an orbit's perifocal frame, defined by `basis`, 
to the inertial frame.
"""
perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, basis::Rotation{3,Float64}) = basis*vec
perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, orb::Orbit) = perifocal_to_inertial_bases(vec, orb.basis)

"""
    rotmat = time_orientation_rotation(t, orb)

Computes the rotation matrix which defines, the prograde, radial, and normal directions
for the orbital state vector of Orbit `orb` at time `t`.
"""
time_orientation_rotation(t, orb::Orbit) = orientation_rotation(time_to_true(t, orb), orb)

"""
    rotmat = align_vectors(vec, vecref)

Computes the rotation matrix to align `vec` with `vecref`.
"""
function align_vectors(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real})
    nvec, nvecref = normalize(vec), normalize(vecref)
    axis = normalize(cross(nvec, nvecref))
    angle = acos(dot(nvec, nvecref))

    return AngleAxis(angle, axis...)
end
