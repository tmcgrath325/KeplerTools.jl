struct Porkchop
    startorbit::Orbit
    endorbit::Orbit
    starttimes::Vector{Float64}
    flighttimes::Vector{Float64}
    # transfers::Matrix{Transfer}
    departure_Δv::Matrix{Float64}
    arrival_Δv::Matrix{Float64}
    Δv::Matrix{Float64}
end

function Porkchop(startorb::Orbit, endorb::Orbit, starttimes::AbstractVector{<:Real}, flighttimes::AbstractVector{<:Real};
                  tfer_fun = Transfer, match_departure_plane=false, match_arrival_plane=false, kwargs...)
    departure_Δv = fill(NaN, length(starttimes), length(flighttimes))
    arrival_Δv = fill(NaN, length(starttimes), length(flighttimes))
    Δv = fill(NaN, length(starttimes), length(flighttimes))
    tferpath=path_to_body(startorb.primary, endorb.primary)
    cmnparent=closest_common_parent(startorb, endorb)
    for (i,st) in enumerate(starttimes)
        for (j,ft) in enumerate(flighttimes)
            tfer = tfer_fun(startorb, endorb, st, st+ft; transferpath=tferpath, commonparent=cmnparent, kwargs...)
            if match_departure_plane
                tfer = match_departure_inclination(tfer)
            end
            if match_arrival_plane
                tfer = match_arrival_inclination(tfer)
            end
            departure_Δv[i,j] = norm(tfer.burns[1].Δv̄)
            arrival_Δv[i,j] = norm(tfer.burns[end].Δv̄)
            Δv[i,j] = tfer.Δv
        end
    end
    # remove NaN and replace with highes value for Δv
    # this is done to avoid JSON serialization problems, and probably should be moved elsewhere
    departure_Δv[isnan.(departure_Δv)] .= maximum(x->isnan(x) ? -Inf : x, departure_Δv)
    arrival_Δv[isnan.(arrival_Δv)] .= maximum(x->isnan(x) ? -Inf : x, arrival_Δv)
    Δv[isnan.(Δv)] .= maximum(x->isnan(x) ? -Inf : x, Δv)
    return Porkchop(startorb, endorb, starttimes, flighttimes, departure_Δv, arrival_Δv, Δv)
end

Porkchop(startorb::Orbit, endorb::Orbit, stime1::Real, stime2::Real, ftime1::Real, ftime2::Real; npts=100, kwargs...
    ) = Porkchop(startorb, endorb, collect(range(stime1, stop=stime2, length=npts)), collect(range(ftime1, stop=ftime2, length=npts)); tfer_fun = Transfer, kwargs...)

function Porkchop(startorb::Orbit, endorb::Orbit, stime1::Real=0., stime2=nothing, ftime1=nothing, ftime2=nothing; kwargs...)
    transferbody = closest_common_parent(startorb, endorb)
    if isnothing(stime2) || isnothing(ftime1) || isnothing(ftime2)
        sorb = startorb
        while sorb.primary!=transferbody
            sorb = sorb.primary.orbit
        end
        eorb = endorb
        while eorb.primary!=transferbody
            eorb = eorb.primary.orbit
        end
    end
    # complete a full rotation of one orbit relative to the other to get start time range
    if isnothing(stime2)
        if isapprox(1/sorb.period, 1/eorb.period)
            stime2 = stime1 + sorb.period
        else
            stime2 = stime1 + 1.25*abs(1/(1/sorb.period - 1/eorb.period))   # a full relative rotation + 25%
        end
    end
    # choose flight times based on the apses of sorb and eorb
    if isnothing(ftime1) || isnothing(ftime2)
        hohmannaxis = max(abs(sorb.a*(1-sorb.e) + eorb.a*(1+eorb.e)), abs(sorb.a*(1+sorb.e) + eorb.a*(1-eorb.e)))
        ftimemid = period(hohmannaxis/2, transferbody.μ)/2
        if isnothing(ftime1) && isnothing(ftime2)
            ftime1 = ftimemid/2
            ftime2 = ftimemid*2
        elseif isnothing(ftime1)
            ftime1 = min(ftimemid/2, ftime2/2)
        else
            ftime2 = max(ftimemid*2, ftime1*2)
        end
    end
    return Porkchop(startorb, endorb, stime1, stime2, ftime1, ftime2; kwargs...)
end

fastPorkchop(args...; kwargs...) = Porkchop(args...; tfer_fun=fastTransfer, kwargs...)

### other methods ###

# returns best start time and end time
function best_start_end_times(pc::Porkchop; field=:Δv)
    (sidx, fidx) = Tuple(findmin(getfield(pc, field))[2])
    return pc.starttimes[sidx], pc.starttimes[sidx]+pc.flighttimes[fidx]
end

function best_transfer(pc::Porkchop; field=:Δv, tfer_fun=Transfer, kwargs...)
    stime, etime = best_start_end_times(pc; field=field)
    return tfer_fun(pc.startorbit, pc.endorbit, stime, etime; kwargs...)
end

### descriptive display ###

Base.show(io::IO, pc::Porkchop) = println(io,
    "Porkchop for transfer from ",
    "orbit around $(pc.startorbit.primary.name) to ",
    "orbit around $(pc.endorbit.primary.name)\n",
    "  Best Δv:\t$(minimum(pc.Δv)) m/s",
)