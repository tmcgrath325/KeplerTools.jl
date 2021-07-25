# Divides each element of the tuple by the specified number
fade_color(color, div=2) = color.÷div

standard_layout(bd::CelestialObject, 
                size=isempty(bd.satellite_bodies) ? bd.SoI : 1.1*maximum(bd.satellite_bodies).orbit.a
                ) = Layout(autosize=false, width=800, height=600,
                           paper_bgcolor = "rgb(30, 30, 30)",
                           margin=attr(l=0, r=0, b=0, t=0),
                           scene=attr(
                                aspectmode="cube",
                                xaxis=attr(visible=false, range=[-size, size]),
                                yaxis=attr(visible=false, range=[-size, size]),
                                zaxis=attr(visible=false, range=[-size, size]))
                            )


standard_hoverlabel(color) = attr(bgcolor = "rgb$color",
                                  font = attr(family = "Courier New, monospace",
                                              size = 10,
                                  ),
                                  align="left",
)

function draw_orbit(orb::Orbit, angles::AbstractVector{<:Real}, starttime; color=(255,255,255), fadedcolor=fade_color(color), name="", timedef=KerbalTime)
    pos_mat = fill(NaN, 3, length(angles))
    vel_mat = fill(NaN, 3, length(angles))
    times = fill(NaN, length(angles))
    tmin = starttime - orb.period/2
    for (i,θ) in enumerate(angles)
        pos_mat[:,i] = orbital_position(θ, orb)
        vel_mat[:,i] = orbital_velocity(θ, orb)
        tmin = true_to_time(θ, orb, tmin)
        times[i] = tmin
    end

    cdata = [[norm(pos_mat[:,i])/1000,
              norm(vel_mat[:,i]),
              seconds_to_datetime(times[i], timedef)...,
              times[i],
              wrap_angle(angles[i]),] for i=1:length(angles)]

    hovertemp = join([  "<b>Position and Speed</b><br>",
                        "   r  = %{customdata[0]:.3e} km<br>",
                        "   v  = %{customdata[1]:.3e} m/s<br><br>",
                        "<b>Date and Time</b><br>",
                        "   Year %{customdata[2]:.0f}, " ,
                        "Day %{customdata[3]:.0f}, " ,
                        "%{customdata[4]:0>2d}:" ,
                        "%{customdata[5]:0>2d}:" ,
                        "%{customdata[6]:0>2d}<br>" ,
                        "   UT:  %{customdata[7]:.3f} s<br>" ,
                        "   θ:   %{customdata[8]:.5f} rad<br><br>" ,
                        "<b>Orbit Parameters</b><br>",
                        "   a  = $(round(orb.a)) m<br>" ,
                        "   e  = $(round(orb.e, digits=4))<br>" ,
                        "   i  = $(round(rad2deg(orb.i), digits=4))°<br>" ,
                        "   ω  = $(round(rad2deg(orb.ω), digits=4))°<br>",
                        "   Ω  = $(round(rad2deg(orb.Ω), digits=4))°<br>",
                        "   Mₒ = $(round(rad2deg(orb.Mo), digits=4)) rad<br>",
                        "   tₒ = $(round(orb.epoch, digits=2)) s"
    ])

    trace = scatter3d(;x=pos_mat[1,:], y=pos_mat[2,:], z=pos_mat[3,:],
                       customdata = cdata,
                       mode="lines", 
                       line=attr(color=angles, 
                                 colorscale = [[0, "rgb$fadedcolor"], 
                                               [1, "rgb$color"]]
                       ),
                       showlegend=false, name=name,
                       hovertemplate=hovertemp,
                       hoverlabel = standard_hoverlabel(color),
    )
    return trace
end

function draw_orbit(orb::Orbit, starttime, endtime=nothing; npts=200, kwargs...)
    # choose end time if not decided (full revolution or until SoI departure)
    if isnothing(endtime)
        if orb.e < 1
            endtime = starttime + orb.period
        else
            if typeof(orb.prim) == star
                furthest_body = maximum(orb.prim.satellitebodies)
                distlim = 1.05 * furthest_body.orb.a * (1+furthest_body.orb.e)
            else
                distlim = orb.prim.SoI
            end
            endtime = true_to_time(orbital_angle(distlim, orb), orb) 
        end
    end

    # trim endtime if the specified times span more than a full revolution
    if (orb.e < 1) && (endtime-starttime > orb.period)
        endtime = starttime + orb.period
    end

    startangle = time_to_true(starttime, orb)
    endangle = wrap_angle(time_to_true(endtime, orb), startangle+1e-3)
    @show startangle, endangle
    @show true_to_time(startangle, orb, starttime-eps(Float64)), true_to_time(endangle, orb, starttime)
    angles = collect(range(startangle, stop=endangle, length=npts))
    return draw_orbit(orb, angles, starttime; kwargs...)
end

draw_orbit(bd::CelestialBody, starttime, endtime=nothing; kwargs...
    ) = draw_orbit(bd.orbit, starttime, endtime; name=bd.name, color=bd.color, kwargs...)
draw_orbit(bd::SpaceObject, starttime, endtime=nothing; kwargs...
    ) = draw_orbit(bd.orbit, starttime, endtime; name=bd.name, kwargs...)

function draw_central_body_orbit(bd::CelestialObject, time; npts=200)
    if typeof(bd) == Star
        return AbstractTrace[]
    end
    θ = time_to_true(time, bd.orbit)
    r = orbital_distance(θ, bd.orbit)
    dlim = isempty(bd.satellite_bodies) ? bd.SoI : 1.1*maximum(bd.satellite_bodies).orbit.a
    arcwidth = dlim / r              # angle in radians
    angles = collect(range(θ-arcwidth, stop=θ+arcwidth, length=npts))

    mid_pos = orbital_position(θ, bd.orbit)
    pos_mat = fill(NaN, 3, length(angles))
    for (i,a) in enumerate(angles)
        pos_mat[:,i] = orbital_position(a, bd.orbit) .- mid_pos
    end

    trace = scatter3d(;x=pos_mat[1,:], y=pos_mat[2,:], z=pos_mat[3,:],
                       mode="lines", 
                       line=attr(color=angles, 
                                 colorscale = [[0, "rgb$(fade_color(bd.color,4))"], 
                                               [1, "rgb$(bd.color)"]]
                       ),
                       showlegend=false,
                       hoverinfo="skip",
    )
    return trace
end