function basis_MRP(Ω, ω, i)
    Ω, ω, i = Float64(Ω), Float64(ω), Float64(i)
    # R₁ = MRP([cos(Ω)  sin(Ω)  0.;    # rotate around z-axis to match Ω
    #          -sin(Ω)  cos(Ω)  0.;
    #           0.      0.      1.])
    # R₂ = MRP([1.      0.      0.;    # rotate around x-axis to match i
    #           0.  cos(i)  sin(i);
    #           0. -sin(i)  cos(i)])
    # R₃ = MRP([cos(ω)  sin(ω)  0.;    # rotate around z-axis to match ω
    #          -sin(ω)  cos(ω)  0.;
    #           0.      0.      1.])
    R₁ = MRP((cos(Ω), -sin(Ω), 0., sin(Ω), cos(Ω), 0., 0., 0., 1.))
    R₂ = MRP((1., 0., 0., 0., cos(i), -sin(i), 0., sin(i), cos(i)))
    R₃ = MRP((cos(ω), -sin(ω), 0., sin(ω), cos(ω), 0., 0., 0., 1.))
    return R₃*R₂*R₁
end

function basis_MRP(xvec::AbstractVector{<:Real}, yvec::AbstractVector{<:Real})
    @assert length(xvec) == length(yvec) == 3
    xbasisvec = xvec ./ norm(xvec)
    ybasisvec = yvec - xbasisvec .* dot(xbasisvec, yvec)
    ybasisvec = ybasisvec ./ norm(ybasisvec)
    zbasisvec = cross(xbasisvec, ybasisvec)
    return MRP(hcat(xbasisvec,ybasisvec,zbasisvec))
end

basis_MRP() = MRP(SMatrix{3,3}(1I))
# basis_MRP(orb::Orbit) = basis_MRP(orb.Ω, orb.ω, orb.i)
# basis_MRP!(orb::Orbit) = get!(orb.attrs, :basis_MRP, basis_MRP(orb))

# methods for switching orbital reference frames

inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, basis::MRP{Float64}) = basis*vec
inertial_to_perifocal_bases(vec::AbstractVector{<:Real}, orb::Orbit) = inertial_to_perifocal_bases(vec, orb.basis)

perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, basis::MRP{Float64}) = basis\vec
perifocal_to_inertial_bases(vec::AbstractVector{<:Real}, orb::Orbit) = perifocal_to_inertial_bases(vec, orb.basis)

# methods for instantaneous orientation (prograde, radial, normal) 
function orientation_MRP(pos, vel, h=cross(pos, vel))
    prograde = vel/norm(vel)
    radial = cross(prograde, h)
    radial = radial/norm(radial)
    return basis_MRP(prograde, radial)
end

orientation_MRP(stvec::StateVector) = orientation_MRP(stvec.position, stvec.velocity)
orientation_MRP(θ, orb::Orbit) = orientation_MRP(StateVector(θ, orb))

time_orientation_MRP(t, orb::Orbit) = orientation_MRP(time_to_true(t, orb), orb)
