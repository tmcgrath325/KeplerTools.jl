function basis_MRP(Ω::Real, ω::Real, i::Real)
    R1 = MRP([cos(Ω)  sin(Ω)  0;    # rotate around z-axis to match Ω
             -sin(Ω)  cos(Ω)  0;
              0       0       1])
    R2 = MRP([1       0       0;    # rotate around x-axis to match i
              0  cos(i)  sin(i);
              0 -sin(i)  cos(i)])
    R3 = MRP([cos(ω)  sin(ω)  0;    # rotate around z-axis to match ω
             -sin(ω)  cos(ω)  0;
              0       0       1])
    return R3*R2*R1
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
basis_MRP(orb::Orbit) = basis_MRP(orb.Ω, orb.ω, orb.i)
basis_MRP!(orb::Orbit) = get!(orb.attrs, :basis_MRP, basis_MRP(orb))

inertial_to_perifocal_bases(vec::SVector{3}, orb::Orbit) = basis_MRP!(orb)*vec
perifocal_to_inertial_bases(vec::SVector{3}, orb::Orbit) = basis_MRP!(orb)\vec
