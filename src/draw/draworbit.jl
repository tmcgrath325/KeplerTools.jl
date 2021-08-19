# Divides each element of the tuple by the specified number
fade_color(color, div=3) = color.÷div

function get_orbit_end_time(orb::Orbit, starttime)
    if orb.e < 1
        return starttime + orb.period
    else
        if typeof(orb.primary) == Star
            furthest_body = maximum(orb.primary.satellite_bodies)
            distlim = 1.05 * furthest_body.orbit.a * (1+furthest_body.orbit.e)
        else
            distlim = orb.primary.SoI
        end
        return true_to_time(orbital_angle(distlim, orb), orb) 
    end
end

orbit_plot_size(bd::CelestialObject; rscale=1.5) = isempty(bd.satellite_bodies) ? bd.SoI : rscale*(maximum(bd.satellite_bodies).orbit.a*(1+maximum(bd.satellite_bodies).orbit.e) + maximum(bd.satellite_bodies).SoI)

orbit_layout(bd::CelestialObject, size=orbit_plot_size(bd)
                ) = Layout(# autosize=false, width=800, height=600,
                           autosize=true,
                           paper_bgcolor = "rgb(30, 30, 30)",
                           plot_bgcolor="rgb(30, 30, 30)",
                           margin=attr(l=0, r=0, b=0, t=0),
                           scene=attr(
                                aspectmode="cube",
                                xaxis=attr(visible=false, range=[-size, size], showgrid=false),
                                yaxis=attr(visible=false, range=[-size, size], showgrid=false),
                                zaxis=attr(visible=false, range=[-size, size], showgrid=false)
                            )
                        )


orbit_hoverlabel(color) = attr(bgcolor = "rgb$color",
                                  font = attr(family = "Courier New, monospace",
                                              size = 10,
                                  ),
                                  align="left",
)

# draw the apoapsis and periapsis of the orbit if they fall within the specified true anomaly
function draw_orbit_apses(orb::Orbit, startangle, endangle; color=(0,125,255), symbol="circle-open")
    if orb.e == 0 return AbstractTrace[] end

    periapsis = startangle <= wrap_angle(0., startangle) <= endangle ? orbital_position(0, orb) : fill(NaN,3)
    apoapsis =  startangle <= wrap_angle(π,  startangle) <= endangle ? orbital_position(π, orb) : fill(NaN,3)
    if count(isnan, [periapsis..., apoapsis...]) > 3
        return AbstractTrace[]
    else
        return [scatter3d(;x=[periapsis[1], apoapsis[1]], 
                           y=[periapsis[2], apoapsis[2]], 
                           z=[periapsis[3], apoapsis[3]],
                           mode="markers", 
                           marker=attr(color="rgb$color", 
                                       symbol=symbol,
                                       size = 4,
                           ),
                           showlegend=false,
                           hoverinfo="skip",
        )]
    end
end

# draw the ascending and descending nodes of the orbit if they fall within the specified true anomaly range
function draw_orbit_nodes(orb::Orbit, startangle, endangle; color=(0,255,125), symbol="circle-open")
    if orb.i == 0 || orb.e > 1 return AbstractTrace[] end

    ascending  = startangle <= wrap_angle(-orb.ω,  startangle) <= endangle ? orbital_position(-orb.ω,  orb) : fill(NaN,3)
    descending = startangle <= wrap_angle(π-orb.ω, startangle) <= endangle ? orbital_position(π-orb.ω, orb) : fill(NaN,3)
    if count(isnan, [ascending..., descending...]) > 3
        return AbstractTrace[]
    else
        return [scatter3d(;x=[ascending[1], descending[1]], 
                           y=[ascending[2], descending[2]], 
                           z=[ascending[3], descending[3]],
                           mode="markers", 
                           marker=attr(color="rgb$color", 
                                       symbol=symbol,
                                       size = 3,
                           ),
                           showlegend=false,
                           hoverinfo="skip",
            )]
    end
end

# draw the orbit positions at each angle, connected by a line
function draw_orbit_path(orb::Orbit, angles::AbstractVector{<:Real}, starttime; color=(255,255,255), name="", timedef=EarthTime)
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
                        "   Mₒ = $(round(orb.Mo, digits=4)) rad<br>",
                        "   tₒ = $(round(orb.epoch, digits=2)) s"
    ])
    return scatter3d(;x=pos_mat[1,:], y=pos_mat[2,:], z=pos_mat[3,:],
                      customdata = cdata,
                      mode="lines", 
                      line=attr(color=angles, 
                                colorscale = [[0, "rgb$(fade_color(color))"], 
                                              [1, "rgb$color"]]
                      ),
                      showlegend=false, name=name,
                      hovertemplate=hovertemp,
                      hoverlabel = orbit_hoverlabel(color),
        )
end

"""
    trace = draw_orbit(orb, starttime, endtime=nothing; kwargs...)

Creates a vector of traces of an orbit `orb` at each true anomaly in `angles`.

When both a start time and end time are provided, the orbit is drawn between the
two times. If the end time is not specified, a full period is drawn. 
"""
function draw_orbit(orb::Orbit, starttime, endtime=nothing; maxpts=100, apses=false, nodes=false, kwargs...)
    # choose end time if not decided (full revolution or until SoI departure)
    if isnothing(endtime)
        endtime = get_orbit_end_time(orb, starttime)
        if orb.e > 1
            if endtime < starttime 
                return AbstractTrace[]
            end
            starttime = true_to_time(-time_to_true(endtime, orb), orb) 
        end
    end

    # trim endtime if the specified times span more than a full revolution
    if (orb.e < 1) && (endtime-starttime > orb.period)
        endtime = starttime + orb.period
    end

    traces = AbstractTrace[]

    startangle = time_to_true(starttime, orb)
    endangle = wrap_angle(time_to_true(endtime, orb), startangle+1e-3)
    npts = max(2, Int(round(maxpts/(2π/(endangle-startangle)))))
    angles = collect(range(startangle, stop=endangle, length=npts))
    push!(traces, draw_orbit_path(orb, angles, starttime; kwargs...))
    if apses
        append!(traces, draw_orbit_apses(orb, startangle, endangle))
    end
    if nodes
        append!(traces, draw_orbit_nodes(orb, startangle, endangle))
    end
    return traces
end

draw_orbit(bd::CelestialBody, starttime, endtime=nothing; kwargs...
    ) = draw_orbit(bd.orbit, starttime, endtime; name=bd.name, color=bd.color, kwargs...)
draw_orbit(bd::SpaceObject, starttime, endtime=nothing; kwargs...
    ) = draw_orbit(bd.orbit, starttime, endtime; name=bd.name, kwargs...)

function draw_central_body_orbit(bd::CelestialObject, time; npts=50, kwargs...)
    if typeof(bd) == Star
        return AbstractTrace[]
    end
    θ = time_to_true(time, bd.orbit)
    r = orbital_distance(θ, bd.orbit)
    dlim = isempty(bd.satellite_bodies) ? bd.SoI : orbit_plot_size(bd; kwargs...)
    arcwidth = dlim / r              # angle in radians
    angles = collect(range(θ-arcwidth, stop=θ+arcwidth, length=npts))

    mid_pos = orbital_position(θ, bd.orbit)
    pos_mat = fill(NaN, 3, length(angles))
    for (i,a) in enumerate(angles)
        pos_mat[:,i] = orbital_position(a, bd.orbit) .- mid_pos
    end
    
    angle_intensity = [wrap_angle(a, θ) for a in angles]
    return  scatter3d(;x=pos_mat[1,:], y=pos_mat[2,:], z=pos_mat[3,:],
                       mode="lines", 
                       line=attr(color=angle_intensity, 
                                 colorscale = [[0, "rgb$(fade_color(bd.color,3))"], 
                                               [1, "rgb$(bd.color)"]]
                       ),
                       showlegend=false,
                       hoverinfo="skip",
        )
end
