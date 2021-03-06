function draw_burn_arrow(brn::Burn; normalize=nothing, scale=nothing, color=(255,0,0), name="Burn")
    if isnothing(normalize)
        normalize = norm(orbital_velocity(0, Orbit(brn.orbit.primary.eqradius+100000, brn.orbit.primary)))
    end
    if isnothing(scale)
        scale = brn.orbit.primary.eqradius
    end

    length = norm(brn.Δv̄)/normalize*scale
    arrowdir = length*brn.Δv̄/norm(brn.Δv̄)
    tail_pos = time_orbital_position(brn.time, brn.orbit)
    head_pos = tail_pos .+ arrowdir
    components = orientation_components(brn)

    tail_trace = scatter3d(;x=[tail_pos[1], head_pos[1]], 
                            y=[tail_pos[2], head_pos[2]], 
                            z=[tail_pos[3], head_pos[3]],
                            mode="lines", 
                            line=attr(color="rgb$color", ),
                            showlegend=false, hoverinfo="skip",
    )
    head_trace = cone(;x=[head_pos[1]],
                       y=[head_pos[2]],
                       z=[head_pos[3]],
                       u=[arrowdir[1]],
                       v=[arrowdir[2]],
                       w=[arrowdir[3]],
                       colorscale = [[0, "rgb$color"],
                                     [1, "rgb$color"]],
                       showscale = false, 
                       showlegend = false,
                       name = name,
                       hovertemplate = join(["$(round(components[1],digits=2)) m/s prograde<br>",
                                             "$(round(components[3],digits=2)) m/s normal<br>",
                                             "$(round(components[2],digits=2)) m/s radial<br><br>",
                                             "UT: $(round(brn.time,digits=2)) s"]),
                       hoverlabel = orbit_hoverlabel(color),
                       
    )
    return [tail_trace, head_trace]
end


function draw_ref_line(r; color=(255,255,255), dash="longdash", width=1)
    trace = scatter3d(;x=[0.,r], y=[0.,0.], z=[0.,0.,],
                       mode="lines",
                       line=attr(color="rgb$color",
                                 dash=dash,
                                 width=width,
                       ),
                       showlegend=false, 
                       hoverinfo="skip",
    )
    return trace
end


function draw_angle_in_plane(vec::AbstractArray{<:Real}, vecref::AbstractArray{<:Real}, basis; npts=100, color=(255,255,255), linedash="dot", arcdash="dash", width=1, rscale=1.25, wrap_min=-π)
    r = rscale*max(norm(vec), norm(vecref))
    θ = angle_in_plane(vecref, basis)
    Δθ = wrap_angle(angle_in_plane(vec, vecref, basis), wrap_min)
    vecplane, vecrefplane = basis\vec, basis\vecref
    arcplane = fill(NaN, 3, npts)
    angles = collect(range(θ, stop=θ+Δθ, length=npts))
    for i=1:npts
        arcplane[:,i] = r*[cos(angles[i]), sin(angles[i]), 0.]
    end
    vecproj, vecrefproj = basis*@SVector([vecplane[1], vecplane[2], 0.]), basis*@SVector([vecrefplane[1], vecrefplane[2], 0.]) 
    arc = basis*arcplane
    
    traces = AbstractTrace[]
    # dashed lines
    xs = arc[1,:]
    prepend!(xs, [vecrefproj[1], arc[1,1]])
    append!(xs, vecproj[1], arc[1,end])
    ys = arc[2,:]
    prepend!(ys, [vecrefproj[2], arc[2,1]])
    append!(ys, vecproj[2], arc[2,end])
    zs = arc[3,:]
    prepend!(zs, [vecrefproj[3], arc[3,1]])
    append!(zs, vecproj[3], arc[3,end])
    push!(traces, scatter3d(;x=xs, 
                             y=ys, 
                             z=zs,
                             mode="lines",
                             line=attr(color="rgb$color",
                                       dash=linedash,
                                       width=width,
                             ),   
                             hoverinfo="skip",           
                             showlegend=false,    
        )
    )
    # text
    sqrtrscale = sqrt(rscale)
    xmid, ymid, zmid = arc[:,Int(round(npts/2))]
    push!(traces, scatter3d(;x=[sqrtrscale*xmid],
                             y=[sqrtrscale*ymid],
                             z=[sqrtrscale*zmid],
                             mode="text",
                             text="$(round(rad2deg(Δθ), digits=2))°",
                             textfont = attr(color="rgb$color",
                                             size=12,
                                             family="Courier New, monospace",
                             ),
                             textposition="middle center",
                             hoverinfo="skip", 
                             showlegend=false,     
    )
)
    return traces
end
draw_angle_in_plane(vec::AbstractArray{<:Real}, orbref::Orbit, t) = draw_angle_in_plane(vec, time_orbital_position(t,orbref), orbref.basis)
draw_angle_in_plane(orb::Orbit, orbref::Orbit, t) = draw_angle_in_plane(time_orbital_position(t, orb), time_orbital_position(t,orbref), orbref.basis)


function draw_orbit_position(orb::Orbit, time; color=(255,255,255), name="", symbol="diamond")
    pos = time_orbital_position(time, orb)
    trace = scatter3d(;x=[pos[1]], y=[pos[2]], z=[pos[3]],
                       mode="markers", 
                       marker=attr(color="rgb$color", 
                                   symbol=symbol,
                                   size = 3,
                       ),
                       showlegend=false, name=name,
                       hovertemplate=name,
                       hoverlabel = orbit_hoverlabel(color),
    )
    return trace
end
