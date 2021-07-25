### type definition ###

struct Transfer
    startorbit::Orbit
    endorbit::Orbit
    departure_time::Float64
    arrival_time::Float64
    transfer_orbit::Orbit
    ejection_orbits::Vector{Orbit}
    insertion_orbits::Vector{Orbit}
    burns::Vector{Burn}
    Δv::Float64
end

### alternate constructors ###

function Transfer(startorb::Orbit, endorb::Orbit, starttime, endtime; dorb_fun = departure_orbit, aorb_fun = arrival_orbit)
    transferpath = path_to_body(startorb.primary, endorb.primary)
    commonparent = closest_common_parent(startorb, endorb)
    tidx = findfirst(x->x==commonparent, transferpath)
    
    # transfer
    sorb = tidx == 1 ? startorb : transferpath[tidx-1].orbit
    eorb = tidx == length(transferpath) ? endorb : transferpath[tidx+1].orbit
    torb = p_lambert(sorb, eorb, starttime, endtime)
    v̄dep = time_state_vector(starttime,torb)[2] .- time_state_vector(starttime,sorb)[2]
    v̄arr = time_state_vector(endtime,  eorb)[2] .- time_state_vector(endtime,  torb)[2]
    
    # backward to start
    dorbs = Orbit[]
    stime = starttime
    didx = tidx-1
    while didx > 0
        sorb = didx == 1 ? startorb : transferpath[didx-1].orbit
        dorb = dorb_fun(sorb, v̄dep, stime)[1]
        push!(dorbs, dorb)
        stime = dorb.epoch
        v̄dep = time_state_vector(stime,dorb)[2] .- time_state_vector(stime,sorb)[2]
        didx = didx - 1
    end

    # forward to end
    aorbs = Orbit[]
    etime = endtime
    aidx = tidx+1
    while aidx < length(transferpath)+1
        eorb = aidx == length(transferpath) ? endorb : transferpath[aidx+1].orbit
        aorb = aorb_fun(eorb, v̄arr, etime)[1]
        push!(aorbs, aorb)
        etime = aorb.epoch
        v̄arr = time_state_vector(etime,eorb)[2] .- time_state_vector(etime,aorb)[2]
        aidx = aidx + 1
    end

    # get burns from beggining and end
    dep_orb = !isempty(dorbs) ? dorbs[end] : torb
    dΔv̄ = time_state_vector(dep_orb.epoch, dep_orb)[2] .- time_state_vector(dep_orb.epoch, startorb)[2]
    departure_burn = Burn(startorb, dΔv̄, dep_orb.epoch)
    
    arr_orb = !isempty(aorbs) ? aorbs[end] : torb
    aΔv̄ = time_state_vector(arr_orb.epoch, endorb)[2] .- time_state_vector(arr_orb.epoch, arr_orb)[2]
    arrival_burn = Burn(endorb, aΔv̄, arr_orb.epoch)

    burns = [departure_burn, arrival_burn]
    Δv = sum([norm(b.Δv̄) for b in burns])

    # construct Transfer 
    return Transfer(startorb,
                    endorb,
                    starttime,
                    endtime,
                    torb,
                    reverse(dorbs),
                    aorbs,
                    burns,
                    Δv,
                    )
end

quickTransfer(startorb::Orbit, endorb::Orbit, starttime, endtime) = Transfer(startorb, endorb, starttime, endtime; dorb_fun=quick_departarrive_orbit, aorb_fun=quick_arrival_orbit)

### helper methods ###

function patch_time(daorb::Orbit, c)
    θₛₚ = c*orbital_angle(daorb.primary.SoI, daorb)
    return true_to_time(θₛₚ, daorb)
end
ejection_time(dorb::Orbit)  = patch_time(dorb,  1)
insertion_time(aorb::Orbit) = patch_time(aorb, -1)

function match_patch_times(tfer::Transfer)
    time_mismatches = Float64[]
    for dorb in tfer.ejection_orbits
        t1 = dorb.epoch
        t2 = ejection_time(dorb)
        append!(time_mismatches(t2-t1))
    end
    for aorb in tfer.insertion_orbits
        t1=insertion_time(aorb)
        t2 = aorb.epoch
        append!(time_mismatches, t2-t1)
    end
    return time_mismatches
end

### descriptive display ###

Base.show(io::IO, tfer::Transfer) = println(io,
    "Transfer from ",
    "orbit around $(tfer.startorbit.primary.name) to ",
    "orbit around $(tfer.endorbit.primary.name)\n",
    "  departure time:\t$(tfer.departure_time) s\n",
    "  arrival time:\t\t$(norm(tfer.arrival_time)) s\n",
    "  Δv:\t\t\t$(tfer.Δv) m/s",
)

