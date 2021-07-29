### type definition ###

"""
A container with information about a series of patched-conic trajectories for a ballistic transfer from a starting orbit
to an ending orbit, where the two orbits are in the same solar system.
"""
struct Transfer
    startorbit::Orbit
    endorbit::Orbit
    departure_time::Float64
    arrival_time::Float64
    patch_positions::Vector{SVector{3,Float64}}
    transfer_orbits::Vector{Orbit}
    ejection_orbits::Vector{Orbit}
    insertion_orbits::Vector{Orbit}
    burns::Vector{Burn}
    Δv::Float64
end

### alternate constructors ###

function Transfer(startorb::Orbit, endorb::Orbit, starttime, endtime,
                  patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}} = nothing;
                  transferpath=path_to_body(startorb.primary, endorb.primary), commonparent=closest_common_parent(startorb, endorb), 
                  torb_fun = p_lambert, dorb_fun = departure_orbit, aorb_fun = arrival_orbit)
    #
    if isnothing(patch_posns)
        patch_posns = fill(@SVector(zeros(3)), length(transferpath)-1)
    end
    tidx = findfirst(x->x==commonparent, transferpath)
    
    # transfer
    sorb = tidx == 1 ? startorb : transferpath[tidx-1].orbit
    eorb = tidx == length(transferpath) ? endorb : transferpath[tidx+1].orbit
    torbs = torb_fun(sorb, eorb, starttime, endtime, patch_posns[tidx-1], patch_posns[tidx])
    if typeof(torbs)<:Orbit
        torbs = [torbs]
    else
        torbs = [torbs...]
    end
    v̄dep = time_orbital_velocity(starttime,torbs[1])   - time_orbital_velocity(starttime,sorb)
    v̄arr = time_orbital_velocity(endtime,  torbs[end]) - time_orbital_velocity(endtime,  eorb)
    
    # backward to start
    dorbs = Orbit[]
    stime = starttime
    didx = tidx-1
    while didx > 0
        sorb = didx == 1 ? startorb : transferpath[didx-1].orbit
        dorb = dorb_fun(sorb, v̄dep, stime)[1]
        if didx>1
            if norm(patch_posns[didx-1])>0
                dorb = departarrive_true_anomaly(time_orbital_position(dorb.epoch, dorb) + patch_posns[didx-1], 
                                                 dorb.primary, 
                                                 v̄dep,
                                                 stime,
                                                 1)[1]
            end
        end
        push!(dorbs, dorb)
        stime = dorb.epoch
        v̄dep = time_orbital_velocity(stime,dorb) - time_orbital_velocity(stime,sorb)
        didx = didx - 1
    end

    # forward to end
    aorbs = Orbit[]
    etime = endtime
    aidx = tidx+1
    while aidx < length(transferpath)+1
        eorb = aidx == length(transferpath) ? endorb : transferpath[aidx+1].orbit
        aorb = aorb_fun(eorb, v̄arr, etime)[1]
        if aidx<=length(patch_posns)
            if norm(patch_posns[aidx])>0
                aorb = departarrive_true_anomaly(time_orbital_position(aorb.epoch, aorb) + patch_posns[aidx], 
                                                 aorb.primary, 
                                                 v̄arr,
                                                 etime,
                                                 -1)[1]
            end
        end
        push!(aorbs, aorb)
        etime = aorb.epoch
        v̄arr = time_orbital_velocity(etime,eorb) - time_orbital_velocity(etime,aorb)
        aidx = aidx + 1
    end

    # get burns from beggining and end
    dep_orb = !isempty(dorbs) ? dorbs[end] : torbs[1]
    dΔv̄ = time_orbital_velocity(dep_orb.epoch, dep_orb) - time_orbital_velocity(dep_orb.epoch, startorb)
    departure_burn = Burn(startorb, dΔv̄, dep_orb.epoch)
    
    arr_orb = !isempty(aorbs) ? aorbs[end] : torbs[end]
    aΔv̄ = time_orbital_velocity(arr_orb.epoch, endorb) - time_orbital_velocity(arr_orb.epoch, arr_orb)
    arrival_burn = Burn(endorb, aΔv̄, arr_orb.epoch)

    burns = [departure_burn]
    for i=1:length(torbs)-1
        Δv̄ = time_orbital_velocity(torbs[i+1].epoch, torbs[i+1]) - time_orbital_velocity(torbs[i+1].epoch, torbs[i])
        push!(burns, Burn(torbs[i], Δv̄, torbs[i+1].epoch))
    end
    push!(burns, arrival_burn)
    
    Δv = sum([norm(b.Δv̄) for b in burns])

    # construct Transfer 
    return Transfer(startorb,
                    endorb,
                    starttime,
                    endtime,
                    patch_posns,
                    torbs,
                    reverse(dorbs),
                    aorbs,
                    burns,
                    Δv)
end

Transfer(tfer::Transfer, patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing; kwargs...
    ) = Transfer(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time, patch_posns; kwargs...)

fastTransfer(startorb::Orbit, endorb::Orbit, starttime, endtime, patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}} = nothing; kwargs...
    ) = Transfer(startorb, endorb, starttime, endtime, patch_posns; dorb_fun=fast_departarrive_orbit, aorb_fun=fast_arrival_orbit, kwargs...)

### helper methods ###

# ensure continuous time across SoI patches

function patch_time_errors(tfer::Transfer)
    up_time_mismatches   = Float64[]
    down_time_mismatches = Float64[]
    for (i,dorb) in enumerate(reverse(tfer.ejection_orbits))
        if i==1
            t1 = tfer.departure_time
        else
            t1 = tfer.ejection_orbits[end-i+2].epoch
        end
        t2 = ejection_time(dorb)
        prepend!(up_time_mismatches, t2-t1)
    end
    for (j,aorb) in enumerate(tfer.insertion_orbits)
        if j==1
            t1 = tfer.arrival_time
        else
            t1 = tfer.insertion_orbits[j-1].epoch
        end
        t2 = insertion_time(aorb)
        append!(down_time_mismatches, t2-t1)
    end
    return up_time_mismatches, down_time_mismatches
end

# function optimize_patch_times(startorb, endorb, stime, etime)
#     function objectivefun(stet) 
#         up_errs, down_errs = patch_time_errors(
#                                 Transfer(startorb,
#                                 endorb,
#                                 stet[1],
#                                 stet[2])
#         )
#         return abs(up_errs[end]) + abs(down_errs[1])
#     end
#     res = optimize(objectivefun, [stime, etime], LBFGS())
#     (stime, etime) = Optim.minimizer(res)
#     return Transfer(startorb, endorb, stime, etime)
# end
# optimize_patch_times(tfer::Transfer
#     ) = optimize_patch_times(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time)

"""
    tfer = match_transfer_patch_times(startorb, endorb, stime, etime, patch_posns; atol=1e-6, maxit=100)
    tfer = match_transfer_patch_times(tfer; kwargs...)

Attempts to optimize transfer start and end times to ensure that patch times are continuous
on either side of a patch, and returns the optimized Transfer.
"""
function match_transfer_patch_times(startorb, endorb, stime, etime, patch_posns; atol=1e-6, maxit=100)
    tfer = Transfer(startorb, endorb, stime, etime, patch_posns)
    up_errs, down_errs = patch_time_errors(tfer)
    tup_err = up_errs[end]
    tdown_err = down_errs[1]
    it=0
    err = abs(tup_err) + abs(tdown_err)
    while err>atol && it<maxit
        it+=1
        stime = stime + tup_err
        etime = etime + tdown_err

        tfer = Transfer(startorb, endorb, stime, etime, patch_posns)
        up_errs, down_errs = patch_time_errors(tfer)
        tup_err = up_errs[end]
        tdown_err = down_errs[1]
        err = abs(tup_err) + abs(tdown_err)
    end
    return tfer
end
match_transfer_patch_times(tfer::Transfer; kwargs...
    ) = match_transfer_patch_times(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time, tfer.patch_positions; kwargs...)

# ensure continuous position across SoI patches

new_patch_positions(tfer::Transfer
    ) = vcat(
            [ejection_position(dorb) for dorb in tfer.ejection_orbits], 
            [insertion_position(aorb) for aorb in tfer.insertion_orbits]
)

function patch_position_errors(tfer::Transfer)
    up_time_distances   = Float64[]
    down_time_distances = Float64[]
    for (i,dorb) in enumerate(reverse(tfer.ejection_orbits))
        if i==1
            r̄1 = time_orbital_position(tfer.departure_time, tfer.transfer_orbits[1])
        else
            r̄1 = time_orbital_position(tfer.ejection_orbits[end-i+2].epoch, tfer.ejection_orbits[end-i+2])
        end
        r̄2 = ejection_position(dorb) + time_orbital_position(ejection_time(dorb), dorb.primary.orbit)
        prepend!(up_time_distances, norm(r̄2.-r̄1))
    end
    for (j,aorb) in enumerate(tfer.insertion_orbits)
        if j==1
            r̄1 = time_orbital_position(tfer.arrival_time, tfer.transfer_orbits[end])
        else
            r̄1 = time_orbital_position(tfer.insertion_orbits[j-1].epoch, tfer.insertion_orbits[j-1])
        end
        r̄2 = insertion_position(aorb) + time_orbital_position(insertion_time(aorb), aorb.primary.orbit)
        append!(down_time_distances, norm(r̄2.-r̄1))
    end
    return up_time_distances, down_time_distances
end

"""
    tfer = match_transfer_patch_positions(tfer; atol=1., maxit=100)

Attempts to optimize transfer start and end times to ensure that patch times and positions are continuous
on either side of a patch, and returns the optimized Transfer.
"""
function match_patch_positions(tfer::Transfer; atol=1., maxit=100)
    err = sum(norm, patch_position_errors(tfer))
    it = 0
    while err > atol && it<maxit
        it += 1
        tfer = match_transfer_patch_times(tfer)
        tfer = Transfer(tfer, new_patch_positions(tfer))
        err = sum(x->sum(abs,x), patch_position_errors(tfer))
    end
    return tfer
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

