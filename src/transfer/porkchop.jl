struct Porkchop
    startorbit::Orbit
    endorbit::Orbit
    starttimes::Vector{Float64}
    flighttimes::Vector{Float64}
    # transfers::Matrix{Transfer}
    Δv::Matrix{Float64}
end

function Porkchop(startorb::Orbit, endorb::Orbit, starttimes::AbstractVector{<:Real}, flighttimes::AbstractVector{<:Real}; tfer_fun = Transfer)
    Δv = fill(NaN, length(starttimes), length(flighttimes))
    tferpath=path_to_body(startorb.primary, endorb.primary)
    cmnparent=closest_common_parent(startorb, endorb)
    for (i,st) in enumerate(starttimes)
        for (j,ft) in enumerate(flighttimes)
            Δv[i,j] = tfer_fun(startorb, endorb, st, st+ft; transferpath=tferpath, commonparent=cmnparent).Δv
        end
    end
    return Porkchop(startorb, endorb, starttimes, flighttimes, Δv)
end

quickPorkchop(startorb::Orbit, endorb::Orbit, starttimes::AbstractVector{<:Real}, flighttimes::AbstractVector{<:Real} 
    ) = Porkchop(startorb, endorb, starttimes, flighttimes; tfer_fun = quickTransfer)

Porkchop(startorb::Orbit, endorb::Orbit, stime1, stime2, ftime1, ftime2; npts=100, tfer_fun = Transfer
    ) = Porkchop(startorb, endorb, collect(range(stime1, stop=stime2, length=npts)), collect(range(ftime1, stop=ftime2, length=npts)); tfer_fun = tfer_fun)
quickPorkchop(startorb::Orbit, endorb::Orbit, stime1, stime2, ftime1, ftime2; npts=250
    ) = Porkchop(startorb, endorb, stime1, stime2, ftime1, ftime2; npts=npts, tfer_fun = quickTransfer)

### descriptive display ###

Base.show(io::IO, pc::Porkchop) = println(io,
    "Porkchop for transfer from ",
    "orbit around $(pc.startorbit.primary.name) to ",
    "orbit around $(pc.endorbit.primary.name)\n",
    "  Best Δv:\t$(minimum(pc.Δv)) m/s",
)
