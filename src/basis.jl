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
    xbasisvec = xvec ./ norm(xvec)
    ybasisvec = yvec - xbasisvec .* dot(xbasisvec, yvec)
    ybasisvec = ybasisvec ./ norm(ybasisvec)
    zbasisvec = cross(xbasisvec, ybasisvec)
    return RotMatrix(hcat(xbasisvec,ybasisvec,zbasisvec))
end

basis_rotation() = one(RotMatrix{3, Float64})
# basis_rotation(orb::Orbit) = basis_rotation(orb.Ω, orb.ω, orb.i)
# basis_rotation!(orb::Orbit) = get!(orb.attrs, :basis_rotation, basis_rotation(orb))

# methods for switching orbital reference frames

inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, basis::Rotation{3,Float64}) = basis\vec
inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, orb::Orbit) = inertial_to_perifocal_bases(vec, orb.basis)

perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, basis::Rotation{3,Float64}) = basis*vec
perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, orb::Orbit) = perifocal_to_inertial_bases(vec, orb.basis)

# methods for instantaneous orientation (prograde, radial, normal) 
function orientation_rotation(pos, vel, h=cross(pos, vel))
    prograde = vel/norm(vel)
    radial = cross(prograde, h)
    radial = radial/norm(radial)
    return basis_rotation(prograde, radial)
end

orientation_rotation(stvec::StateVector) = orientation_rotation(stvec.position, stvec.velocity)
orientation_rotation(θ, orb::Orbit) = orientation_rotation(StateVector(θ, orb))

time_orientation_rotation(t, orb::Orbit) = orientation_rotation(time_to_true(t, orb), orb)

function align_vectors(vec::AbstractVector{<:Real}, vecref::AbstractVector{<:Real})
    nvec, nvecref = normalize(vec), normalize(vecref)
    axis = normalize(cross(nvec, nvecref))
    angle = acos(dot(nvec, nvecref))

    return AngleAxis(angle, axis...)
end
