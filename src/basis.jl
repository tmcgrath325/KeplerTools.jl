function basis_matrix(Ω, ω, i)
    R1 = MRP(@SMatrix([cos(Ω)  sin(Ω)  0;                   # rotate around z-axis to match Ω
                      -sin(Ω)  cos(Ω)  0;
                       0       0       1]))
    R2 = MRP(@SMatrix([1       0       0;                   # rotate around x-axis to match i
                       0  cos(i)  sin(i);
                       0 -sin(i)  cos(i)]))
    R3 = MRP(@SMatrix([cos(ω)  sin(ω)  0;                   # rotate around z-axis to match ω
                      -sin(ω)  cos(ω)  0;
                       0       0       1]))
    return R3*R2*R1
end

basis_matrix() = MRP(SMatrix{3,3}(1I))
basis_matrix(orb::Orbit) = basis_matrix(orb.Ω, orb.ω, orb.i)
basis_matrix!(orb::Orbit) = get!(orb.attrs, :bases, basis_matrix(orb))

inert_to_perif_bases(vec::SVector{3}, orb::Orbit,) = basis_matrix!(orb)*vec
perif_to_inert_bases(vec::SVector{3}, orb::Orbit) = basis_matrix!(orb)\vec
