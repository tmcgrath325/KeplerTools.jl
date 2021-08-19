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
    planechange::Bool
    match_sorb_Mo::Bool
    match_eorb_Mo::Bool
end

### alternate constructors ###

function Transfer(startorb::Orbit, endorb::Orbit, starttime, endtime,
                  patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}} = nothing;
                  transferpath = path_to_body(startorb.primary, endorb.primary), 
                  commonparent = closest_common_parent(startorb, endorb), 
                  torb_fun::Union{typeof(p_lambert), typeof(plane_change_p_lambert)} = p_lambert, 
                  dorb_fun::Union{typeof(departure_orbit), typeof(fast_departure_orbit)} = departure_orbit, 
                  aorb_fun::Union{typeof(arrival_orbit), typeof(fast_arrival_orbit)} = arrival_orbit,
                  match_sorb_Mo=true, match_eorb_Mo=true,)
    
    tidx = findfirst(x->x==commonparent, transferpath)
    n_ejc = tidx - 1
    n_ins = length(transferpath) - tidx
    if isnothing(patch_posns)
        patch_posns = fill(@SVector(zeros(3)), length(transferpath)-1)
    end
    
    # transfer
    sorb = tidx == 1 ? startorb : transferpath[tidx-1].orbit
    eorb = tidx == length(transferpath) ? endorb : transferpath[tidx+1].orbit
    start_pos = n_ejc > 0 ? patch_posns[tidx-1] : @SVector(zeros(3))
    end_pos = n_ins > 0 ? patch_posns[tidx] : @SVector(zeros(3))
    torbs = torb_fun(sorb, eorb, starttime, endtime, start_pos, end_pos)
    if typeof(torbs)<:Orbit
        torbs = Orbit[torbs]
    else
        torbs = Orbit[torbs...]
    end
    v̄dep = time_orbital_velocity(starttime,torbs[1])   - time_orbital_velocity(starttime,sorb)
    v̄arr = time_orbital_velocity(endtime,  torbs[end]) - time_orbital_velocity(endtime,  eorb)
    
    # backward to start
    dorbs = Orbit[]
    stime = starttime
    didx = tidx-1
    while didx > 0
        sorb = didx == 1 ? startorb : transferpath[didx-1].orbit
        if didx==1 && !match_sorb_Mo
            dorb = dorb_fun(startorb, v̄dep, stime; match_pkorb_Mo=false)[1]
            r̄ₒ = time_orbital_position(dorb.epoch, dorb)
            θₚₖ = angle_in_plane(r̄ₒ, startorb)
            startorb = Orbit(startorb.a, startorb.e, startorb.i, startorb.Ω, startorb.ω, 
                             true_to_mean(θₚₖ, startorb), 
                             dorb.epoch, startorb.primary,
            )
        else
            dorb = dorb_fun(sorb, v̄dep, stime)[1]
        end
        if didx>1
            if norm(patch_posns[didx-1])>0
                depoch = dorb.epoch
                dorb = departarrive_true_anomaly(time_orbital_position(dorb.epoch, dorb) + patch_posns[didx-1], 
                                                 dorb.primary, 
                                                 v̄dep,
                                                 stime,
                                                 1)[1]
                dorb = Orbit(dorb.a, dorb.e, dorb.i, dorb.Ω, dorb.ω, dorb.Mo, depoch, dorb.primary)
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
        if aidx==length(transferpath) && !match_eorb_Mo
            aorb = aorb_fun(endorb, v̄arr, etime; match_pkorb_Mo=false)[1]
            r̄ₒ = time_orbital_position(aorb.epoch, aorb)
            θₚₖ = angle_in_plane(r̄ₒ, endorb)
            endorb = Orbit(endorb.a, endorb.e, endorb.i, endorb.Ω, endorb.ω, 
                           true_to_mean(θₚₖ, endorb), 
                           aorb.epoch, endorb.primary,
            )
        else
            aorb = aorb_fun(eorb, v̄arr, etime)[1]
        end
        if aidx<=length(patch_posns)
            if norm(patch_posns[aidx])>0
                aepoch = aorb.epoch
                aorb = departarrive_true_anomaly(time_orbital_position(aorb.epoch, aorb) + patch_posns[aidx], 
                                                 aorb.primary, 
                                                 v̄arr,
                                                 etime,
                                                 -1)[1]
                aorb = Orbit(aorb.a, aorb.e, aorb.i, aorb.Ω, aorb.ω, aorb.Mo, aepoch, aorb.primary)
            end
        end
        push!(aorbs, aorb)
        etime = aorb.epoch
        v̄arr = time_orbital_velocity(etime,aorb) - time_orbital_velocity(etime,eorb)
        aidx = aidx + 1
    end

    # get burns from beggining and end
    dep_orb, dep_time = !isempty(dorbs) ? (dorbs[end], dorbs[end].epoch) : (torbs[1], starttime)
    dΔv̄ = time_orbital_velocity(dep_orb.epoch, dep_orb) - time_orbital_velocity(dep_orb.epoch, startorb)
    departure_burn = Burn(startorb, dΔv̄, dep_orb.epoch)
    
    arr_orb, arr_time = !isempty(aorbs) ? (aorbs[end], aorbs[end].epoch) : (torbs[end], endtime)
    aΔv̄ = time_orbital_velocity(arr_time, endorb) - time_orbital_velocity(arr_time, arr_orb)
    arrival_burn = Burn(endorb, aΔv̄, arr_time)

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
                    Δv,
                    torb_fun==plane_change_p_lambert,
                    match_sorb_Mo,
                    match_eorb_Mo)
end

Transfer(tfer::Transfer, patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing; kwargs...
    ) = Transfer(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time, patch_posns; 
                 torb_fun = tfer.planechange ? plane_change_p_lambert : p_lambert,
                 match_sorb_Mo=tfer.match_sorb_Mo, match_eorb_Mo=tfer.match_eorb_Mo, 
                 # match_dorb_plane=tfer.match_dorb_plane, match_aorb_plane=tfer.match_aorb_plane,
                 kwargs...)

fastTransfer(startorb::Orbit, endorb::Orbit, starttime, endtime, patch_posns::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}} = nothing; kwargs...
    ) = Transfer(startorb, endorb, starttime, endtime, patch_posns; dorb_fun=fast_departure_orbit, aorb_fun=fast_arrival_orbit, kwargs...)

### helper methods ###

# ensure continuous time across SoI patches

function patch_time_errors(tfer::Transfer)
    up_time_mismatches   = Float64[]
    down_time_mismatches = Float64[]
    for (i,dorb) in enumerate(tfer.ejection_orbits)
        if i==length(tfer.ejection_orbits)
            t1 = tfer.departure_time
        else
            t1 = tfer.ejection_orbits[i+1].epoch
        end
        t2 = ejection_time(dorb)
        push!(up_time_mismatches, t2-t1)
    end
    for (j,aorb) in enumerate(tfer.insertion_orbits)
        if j==1
            t1 = tfer.arrival_time
        else
            t1 = tfer.insertion_orbits[j-1].epoch
        end
        t2 = insertion_time(aorb)
        push!(down_time_mismatches, t2-t1)
    end
    return up_time_mismatches, down_time_mismatches
end

function transfer_patch_time_errors(tfer::Transfer)
    up_time_mismatches, down_time_mismatches = patch_time_errors(tfer)
    start_err = isempty(up_time_mismatches)   ? 0. : up_time_mismatches[end]
    end_err =   isempty(down_time_mismatches) ? 0. : down_time_mismatches[1]
    return up_time_mismatches[end], down_time_mismatches[1]
end

"""
    tfer = match_transfer_patch_times(startorb, endorb, stime, etime, patch_posns; atol=1e-6, maxit=100)
    tfer = match_transfer_patch_times(tfer; kwargs...)

Attempts to optimize transfer start and end times to ensure that patch times are continuous
on either side of a patch, and returns the optimized Transfer.
"""
function match_transfer_patch_times(startorb, endorb, stime, etime, patch_posns; atol=1e-3, maxit=100, kwargs...)
    start_err, end_err = transfer_patch_time_errors(Transfer(startorb, endorb, stime, etime, patch_posns; kwargs...))
    err = sum(abs, [start_err, end_err])
    it=0
    # naive iteration
    while err>atol && it<maxit
        it+=1
        stime = stime + start_err
        etime = etime + end_err

        start_err, end_err = transfer_patch_time_errors(Transfer(startorb, endorb, stime, etime, patch_posns; kwargs...))
        err = sum(abs, [start_err, end_err])
    end
    # # attempt at optimization 
    # function err_fun(se, startorb, endorb, patch_posns; kwargs...)
    #     start_err, end_err = transfer_patch_time_errors(Transfer(startorb, endorb, se[1], se[2], patch_posns; kwargs...))
    #     return sum(abs, [start_err, end_err])
    # end
    # obj_fun(se) = err_fun(se, startorb, endorb, patch_posns; kwargs...)
    # res = optimize(obj_fun, [stime, etime], SimulatedAnnealing())
    # @show res.minimum
    # stime, etime = res.minimizer
    return Transfer(startorb, endorb, stime, etime, patch_posns; kwargs...)
end
match_transfer_patch_times(tfer::Transfer; kwargs...
    ) = match_transfer_patch_times(tfer.startorbit, tfer.endorbit, tfer.departure_time, tfer.arrival_time, tfer.patch_positions; 
                                   torb_fun = tfer.planechange ? plane_change_p_lambert : p_lambert,
                                   match_sorb_Mo=tfer.match_sorb_Mo, match_eorb_Mo=tfer.match_eorb_Mo, 
                                   # match_dorb_plane=tfer.match_dorb_plane, match_aorb_plane=tfer.match_aorb_plane,
                                   kwargs...)

# ensure continuous position across SoI patches

function new_patch_positions(tfer::Transfer, coef=1)
    if tfer.startorbit.primary == tfer.endorbit.primary
       return nothing
    end
    if coef == 1
        return vcat(SVector{3,Float64}[ejection_position(dorb) for dorb in tfer.ejection_orbits], 
                    SVector{3,Float64}[insertion_position(aorb) for aorb in tfer.insertion_orbits])
    else 
        return vcat(SVector{3,Float64}[normalize(coef*ejection_position(dorb)+(1-coef)*tfer.patch_positions[i])*norm(ejection_position(dorb)) for (i,dorb) in enumerate(tfer.ejection_orbits)], 
                    SVector{3,Float64}[normalize(coef*insertion_position(aorb)+(1-coef)*tfer.patch_positions[length(tfer.ejection_orbits)+j])*norm(insertion_position(aorb)) for (j,aorb) in enumerate(tfer.insertion_orbits)])
    end
end

function patch_position_errors(tfer::Transfer)
    if tfer.startorbit.primary == tfer.endorbit.primary
        return 0.
     end
    return tfer.patch_positions - new_patch_positions(tfer)
end

"""
    tfer = match_transfer_patch_positions(tfer; atol=1., maxit=100)

Attempts to optimize transfer start and end positions to ensure that patch times and positions are continuous
on either side of a patch, and returns the optimized Transfer.
"""
function match_patch_positions(tfer::Transfer; atol=1., maxit=100, kwargs...)
    err = sum(norm, patch_position_errors(tfer))
    best_err = err
    best_tfer = tfer
    coef = 1
    it = 0
    while best_err > atol && it<maxit
        it += 1
        tfer = match_transfer_patch_times(tfer; kwargs...)
        tfer = Transfer(tfer, new_patch_positions(tfer, coef); kwargs...)
        err = sum(norm, patch_position_errors(tfer))
        if err < best_err
            best_err = err
            best_tfer = tfer
        end
    end
    return best_tfer
end

# align start orbit with departure orbit plane
function match_departure_inclination(tfer::Transfer)
    # align start orbit with departure orbit plane
    if isempty(tfer.ejection_orbits)
        return tfer
    end
    rot = align_vectors(tfer.startorbit.basis[:,3], tfer.ejection_orbits[1].basis[:,3])
    r̄ₒ, v̄ₒ = time_state_vector(tfer.ejection_orbits[1].epoch, tfer.startorbit)
    v̄ₒ = rot * v̄ₒ
    startorb = Orbit(tfer.ejection_orbits[1].epoch, r̄ₒ, v̄ₒ, tfer.startorbit.primary, tfer.startorbit.epoch)

    dΔv̄ = time_orbital_velocity(tfer.ejection_orbits[1].epoch, tfer.ejection_orbits[1]) - time_orbital_velocity(tfer.ejection_orbits[1].epoch, startorb)
    departure_burn = Burn(startorb, dΔv̄, tfer.ejection_orbits[1].epoch)
    burns = [departure_burn, tfer.burns[2:end]...]

    Δv = sum([norm(b.Δv̄) for b in burns])

    return Transfer(startorb, tfer.endorbit, tfer.departure_time, tfer.arrival_time, 
                    tfer.patch_positions, tfer.transfer_orbits, tfer.ejection_orbits, tfer.insertion_orbits,
                    burns, Δv,
                    tfer.planechange, tfer.match_sorb_Mo, tfer.match_eorb_Mo)
end

# align end orbit with arrival orbit plane
function match_arrival_inclination(tfer::Transfer)
    # align start orbit with departure orbit plane
    if isempty(tfer.insertion_orbits)
        return tfer
    end
    rot = align_vectors(tfer.endorbit.basis[:,3], tfer.insertion_orbits[end].basis[:,3])
    r̄ₒ, v̄ₒ = time_state_vector(tfer.insertion_orbits[end].epoch, tfer.endorbit)
    v̄ₒ = rot * v̄ₒ
    endorb = Orbit(tfer.insertion_orbits[end].epoch, r̄ₒ, v̄ₒ, tfer.endorbit.primary, tfer.endorbit.epoch)

    aΔv̄ = time_orbital_velocity(tfer.insertion_orbits[end].epoch, endorb) - time_orbital_velocity(tfer.insertion_orbits[end].epoch, tfer.insertion_orbits[end])
    arrival_burn = Burn(tfer.insertion_orbits[end], aΔv̄, tfer.insertion_orbits[end].epoch)
    burns = [tfer.burns[1:end-1]..., arrival_burn]

    Δv = sum([norm(b.Δv̄) for b in burns])

    return Transfer(tfer.startorbit, endorb, tfer.departure_time, tfer.arrival_time, 
                    tfer.patch_positions, tfer.transfer_orbits, tfer.ejection_orbits, tfer.insertion_orbits,
                    burns, Δv,
                    tfer.planechange, tfer.match_sorb_Mo, tfer.match_eorb_Mo)
end

## Get additional inforation about transfers that won't be stored as part of the type ##

function phase_angle(tfer::Transfer)
    sorb = isempty(tfer.ejection_orbits)  ? tfer.startorbit : tfer.ejection_orbits[end].primary.orbit
    eorb = isempty(tfer.insertion_orbits) ? tfer.endorbit   : tfer.insertion_orbits[1].primary.orbit
    r̄ₛ = time_orbital_position(tfer.departure_time, sorb)
    r̄ₑ = time_orbital_position(tfer.departure_time, eorb)
    return wrap_angle(angle_in_plane(r̄ₑ, r̄ₛ, sorb.basis), -π)
end

function ejection_angle(tfer::Transfer)
    if isempty(tfer.ejection_orbits)
        return nothing
    end
    v̂ₚ = normalize(time_orbital_velocity(tfer.ejection_orbits[1].epoch, tfer.startorbit.primary.orbit))
    r̄ₒ = time_orbital_position(tfer.ejection_orbits[1].epoch, tfer.ejection_orbits[1])
    wrap_angle(angle_in_plane(r̄ₒ, v̂ₚ, tfer.startorbit.basis), -π)
end

function transfer_angle(tfer::Transfer)
    r̄ₛ = time_orbital_position(tfer.departure_time, tfer.transfer_orbits[1])
    r̄ₑ = time_orbital_position(tfer.arrival_time, tfer.transfer_orbits[end])
    return wrap_angle(angle_in_plane(r̄ₑ, r̄ₛ, basis_rotation(r̄ₛ, r̄ₑ)))
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

