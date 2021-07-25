function draw_porkchop(pc::Porkchop; levelscale=1.1, nlevels=16, timedef=KerbalTime)
    minΔv = minimum(pc.Δv)
    lvls = minΔv * [levelscale .^ i for i=0:nlevels]
    logΔv = fill(NaN, size(pc.Δv)...) 
    for (i,el) in enumerate(pc.Δv)
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
            # tickcolor = ,
            tickfont = attr(
                family = "Courier New, monospace",
            ),
        ),
        customdata = pc.Δv,
        hovertemplate = "Δv = %{customdata:.2f} m/s<extra></extra>",
        # xaxis = attr(
        #     title_text = "Transfer start (day #)",
        # ),
        # yaxis = attr(
        #     title_text="Transfer duration (days)",
        # ),
    )
    return trace
end