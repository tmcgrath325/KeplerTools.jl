function draw_wireframe_sphere(pos, r; npts=11, color=(255,255,255), name=nothing, opacity=1.)
    ϕs = range(0, stop=2π, length=npts)
    θs = range(-π/2, stop=π/2, length=npts)
    # vertical
    v_xs = [[r * cos(θ) * sin(ϕ) + pos[1] for θ in θs] for ϕ in ϕs[1:end-1]]
    v_ys = [[r * cos(θ) * cos(ϕ) + pos[2] for θ in θs] for ϕ in ϕs[1:end-1]]
    v_zs = [[r * sin(θ) + pos[3] for θ in θs] for ϕ in ϕs[1:end-1]]
    # horizontal
    h_xs = [[r * cos(θ) * sin(ϕ) + pos[1] for ϕ in ϕs] for θ in θs]
    h_ys = [[r * cos(θ) * cos(ϕ) + pos[2] for ϕ in ϕs] for θ in θs]
    h_zs = [[r * sin(θ) + pos[3] for ϕ in ϕs] for θ in θs] 

    hover = isnothing(name) ? "skip" : nothing

    xs, ys, zs = Float64[], Float64[], Float64[]
    for i=1:length(v_xs)
        append!(xs, v_xs[i])
        push!(xs, NaN)
        append!(ys, v_ys[i])
        push!(ys, NaN)
        append!(zs, v_zs[i])
        push!(zs, NaN)
    end
    for j=1:length(h_xs)
        append!(xs, h_xs[j])
        push!(xs, NaN)
        append!(ys, h_ys[j])
        push!(ys, NaN)
        append!(zs, h_zs[j])
        push!(zs, NaN)
    end

    return scatter3d(;x=xs, y=ys, z=zs, 
                      mode="lines",
                      line=attr(color="rgb$color"),
                      opacity=opacity,
                      showlegend=false, showscale=false, name=name,
                      hovertemplate=name, hoverinfo=hover,
            )
end

draw_central_body(bd::Union{CelestialBody,Star}; kwargs...
    ) = draw_wireframe_sphere(SVector{3,Float64}(0.,0.,0.), bd.eqradius; color=bd.color, name=bd.name, kwargs...)

draw_orbiting_body(bd::CelestialBody, t; kwargs...
    ) = draw_wireframe_sphere(time_orbital_position(t, bd.orbit), bd.eqradius; color=bd.color, name=bd.name, kwargs...)

draw_soi(bd::CelestialBody, t; kwargs...
    ) = draw_wireframe_sphere(time_orbital_position(t, bd.orbit), bd.SoI; color=bd.color, name=bd.name, opacity=0.25, kwargs...)