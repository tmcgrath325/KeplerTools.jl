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

function Transfer(startorb::Orbit, endorb::Orbit, starttime, endtime;
                  transferpath=path_to_body(startorb.primary, endorb.primary), commonparent=closest_common_parent(startorb, endorb), 
                  dorb_fun = departure_orbit, aorb_fun = arrival_orbit)
    
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

quickTransfer(startorb::Orbit, endorb::Orbit, starttime, endtime; kwargs...
    ) = Transfer(startorb, endorb, starttime, endtime; dorb_fun=quick_departarrive_orbit, aorb_fun=quick_arrival_orbit, kwargs...)

### helper methods ###

function patch_time(daorb::Orbit, c)
    θₛₚ = c*orbital_angle(daorb.primary.SoI, daorb)
    return true_to_time(θₛₚ, daorb)
end
ejection_time(dorb::Orbit)  = patch_time(dorb,  1)
insertion_time(aorb::Orbit) = patch_time(aorb, -1)

function patch_time_errors(tfer::Transfer)
    up_time_mismatches   = Float64[]
    down_time_mismatches = Float64[]
    for (i,dorb) in enumerate(reverse(tfer.ejection_orbits))
        if i==1
            t1 = tfer.departure_time
            t2 = ejection_time(dorb)
            prepend!(up_time_mismatches, t2-t1)
        else
            t1 = tfer.ejection_orbits[end-i+2].epoch
            t2 = ejection_time(dorb)
            prepend!(up_time_mismatches, t2-t1 ) # - up_time_mismatches[1])
        end
    end
    for (j,aorb) in enumerate(tfer.insertion_orbits)
        if j==1
            t1 = tfer.arrival_time
            t2 = insertion_time(aorb)
            append!(down_time_mismatches, t2-t1)
        else
            t1 = tfer.insertion_orbits[j-1].epoch
            t2 = insertion_time(aorb)
            append!(down_time_mismatches, t2-t1 ) # - down_time_mismatches[j-1])
        end
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
#         @show up_errs[end], down_errs[1]
#         return abs(up_errs[end]) + abs(down_errs[1])
#     end
#     res = optimize(objectivefun, [stime, etime], LBFGS())
#     (stime, etime) = Optim.minimizer(res)
#     @show Optim.minimum(res)
#     return Transfer(startorb, endorb, stime, etime)
# end
# optimize_patch_times(tfer::Transfer
#     ) = optimize_patch_times(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time)


function match_transfer_patch_times(startorb, endorb, stime, etime; atol=1e-6, maxit=100)
    err = 1+atol
    it=0
    tfer = quickTransfer(startorb, endorb, stime, etime)
    up_errs, down_errs = patch_time_errors(tfer)
    tup_err = up_errs[end]
    tdown_err = down_errs[1]
    err = abs(tup_err) + abs(tdown_err)
    while err>atol && it<maxit
        it+=1
        stime = stime + tup_err
        etime = etime + tdown_err

        tfer = quickTransfer(startorb, endorb, stime, etime)
        up_errs, down_errs = patch_time_errors(tfer)
        tup_err = up_errs[end]
        tdown_err = down_errs[1]
        err = abs(tup_err) + abs(tdown_err)
        @show it, up_errs, down_errs
    end
    return tfer
end
match_transfer_patch_times(tfer::Transfer
    ) = match_transfer_patch_times(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time)

### descriptive display ###

Base.show(io::IO, tfer::Transfer) = println(io,
    "Transfer from ",
    "orbit around $(tfer.startorbit.primary.name) to ",
    "orbit around $(tfer.endorbit.primary.name)\n",
    "  departure time:\t$(tfer.departure_time) s\n",
    "  arrival time:\t\t$(norm(tfer.arrival_time)) s\n",
    "  Δv:\t\t\t$(tfer.Δv) m/s",
)

