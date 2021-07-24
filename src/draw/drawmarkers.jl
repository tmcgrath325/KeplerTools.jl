function draw_burn_arrow(brn::Burn; normalize=nothing, scale=nothing, color=(255,0,0), name="Burn")
    if isnothing(normalize)
        normalize = norm(orbital_velocity(0, Orbit(brn.orbit.primary.eqradius+100000, brn.orbit.primary)))
    end
    if isnothing(scale)
        scale = brn.orbit.primary.eqradius
    end

    length = norm(brn.Δv̄)/normalize*scale
    @show length
    arrowdir = length*brn.Δv̄/norm(brn.Δv̄)
    tail_pos = time_orbital_position(brn.time, brn.orbit)
    head_pos = tail_pos .+ arrowdir
    @show head_pos
    @show arrowdir
    components = time_orientation_MRP(brn.time, brn.orbit) * brn.Δv̄

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
                       hoverlabel = standard_hoverlabel(color),
                       
    )
    return [tail_trace, head_trace]
end