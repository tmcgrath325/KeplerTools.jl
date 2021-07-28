standard_porkchop_layout(
    ) = Layout(autosize=false, width=800, height=600,
            paper_bgcolor = "rgb(30, 30, 30)",
            font = attr(family = "Courier New, monospace",
                        size = 10,
                        color = "rgb(200, 200, 200)"
            ), 
            xaxis=attr(title_text="Transfer start (day #)",
                       showspikes=true,
                       spikemode="across",
                       spikecolor="rgb(200, 200, 200)",
                       spikedash="solid",
                       spikethickness=-1,
            ),
            yaxis=attr(title_text="Transfer duration (days)",
                       showspikes=true,
                       spikemode="across",
                       spikecolor="rgb(200, 200, 200)",
                       spikedash="solid",
                       spikethickness=-1,
            ),
)

"""
    trace = draw_porkchop(pc; kwargs...)

Creates a trace of a porkchop plot.
"""
function draw_porkchop(pc::Porkchop; levelscale=1.1, nlevels=16, timedef=KerbalTime, plottype=:Δv)
    minΔv = minimum(getfield(pc, plottype))
    lvls = minΔv * [levelscale .^ i for i=0:nlevels]
    logΔv = fill(NaN, size(getfield(pc, plottype))...) 
    for (i,el) in enumerate(getfield(pc, plottype))
        logΔv[i] = log(el)/log(levelscale)
    end
    loglvls = [log(lvl)/log(levelscale) for lvl in lvls]
    @show loglvls
    colorlabels = ["$(Int(floor(lvl)))" for lvl in lvls]

    day_to_sec = timedef["day"] * timedef["hour"] * timedef["minute"]

    trace = contour(
        z = logΔv,
        x = pc.starttimes ./ day_to_sec,
        y = pc.flighttimes ./ day_to_sec,
        contours = Dict(
            :start => loglvls[1]+0.05,
            :end => loglvls[end],
            :size => (loglvls[end]-loglvls[1])/nlevels,
            # :showlines => false,
        ),
        colorscale = "Viridis",
        reversescale = true,
        colorbar = attr(
            tickvals = loglvls,
            ticktext = colorlabels,
            title = attr(
                text = "Δv (m/s)",
                font = attr(
                    family = "Courier New, monospace",
                ),  
            ), 
            tickfont = attr(
                family = "Courier New, monospace",
            ),
        ),
        customdata = getfield(pc, plottype),
        hovertemplate = "Δv = %{customdata:.2f} m/s<extra></extra>",
    )
    return trace
end