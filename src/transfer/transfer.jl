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

function Transfer(startorb::Orbit, endorb::Orbit, starttime, endtime)
    transferpath = path_to_body(startorb.primary, endorb.primary)
    commonparent = closest_common_parent(startorb, endorb)
    tidx = findfirst(x->x==commonparent, transferpath)
    
    # transfer
    sorb = tidx == 1 ? startorb : transferpath[tidx-1].orbit
    eorb = tidx == length(transferpath) ? endorb : transferpath[tidx+1].orbit
    torb = p_lambert(sorb, eorb, starttime, endtime)
    v̄dep = time_state_vector(starttime,torb)[2] .- time_state_vector(starttime,sorb)[2]
    v̄arr = time_state_vector(endtime,eorb)[2] .- time_state_vector(endtime,torb)[2]
    
    # backward to start
    dorbs = Orbit[]
    stime = starttime
    didx = tidx-1
    while didx > 0
        sorb = didx == 1 ? startorb : transferpath[didx-1].orbit
        dorb = departure_orbit(sorb, v̄dep, stime)[1]
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
        aorb = arrival_orbit(eorb, v̄arr, etime)[1]
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

function math_patch_times()
end