function plot_survival_curves(
    times,
    observed,
    simulated;
    savepath=nothing
)
    colors = (
        low = :blue,
        mid = :orange,
        high = :green
    )

    p = plot(
        times,
        observed.low,
        label = "Low",
        color = colors.low,
        linestyle = :solid,
        linewidth = 3,
        xlabel = "Days after Dox removal",
        ylabel = "CML-free survival",
        yscale = :log10,
        ylim = (0.025, 1.05),
        xlim = (0, 270)
    )

    plot!(
        p,
        times,
        observed.mid,
        label = "Mid",
        color = colors.mid,
        linewidth = 3
    )

    plot!(
        p,
        times,
        observed.high,
        label = "High",
        color = colors.high,
        linewidth = 3
    )

    plot!(
        p,
        times,
        simulated.low,
        label = false,
        color = colors.low,
        linestyle = :dot,
        linewidth = 2
    )

    plot!(
        p,
        times,
        simulated.mid,
        label = false,
        color = colors.mid,
        linestyle = :dot,
        linewidth = 2
    )

    plot!(
        p,
        times,
        simulated.high,
        label = false,
        color = colors.high,
        linestyle = :dot,
        linewidth = 2
    )
    # add the simulated points
    scatter!(times, simulated.low, color = colors.low, label = false)
    scatter!(times, simulated.mid, color = colors.mid, label = false)
    scatter!(times, simulated.high, color = colors.high, label = false)

    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end


function plot_chimerism_curves(
    times,
    prog,
    nonprog;
    savepath=nothing
)
    colors = (
        low = :blue,
        mid = :orange,
        high = :green
    )

    p = plot(
        times,
        nonprog.low,
        label = "Low",
        color = colors.low,
        linestyle = :dot,
        linewidth = 3,
        xlabel = "Days after Dox removal",
        ylabel = "Mean Chimerism",
        yscale = :log10,
        ylim = (0.001, 1.05),
        xlim = (0, 270)
    )

    plot!(
        p,
        times,
        prog.low,
        linestyle = :solid,
        label = false,
        color = colors.low,
        linewidth = 3
    )

    plot!(
        p,
        times,
        nonprog.mid,
        label = "mid",
        linestyle = :dot,
        color = colors.mid,
        linewidth = 3
    )

    plot!(
        p,
        times,
        prog.mid,
        label = false,
        linestyle = :solid,
        color = colors.mid,
        linewidth = 3
    )

    plot!(
        p,
        times,
        nonprog.high,
        label = "high",
        linestyle = :dot,
        color = colors.high,
        linewidth = 3
    )

    plot!(
        p,
        times,
        prog.high,
        label = false,
        linestyle = :solid,
        color = colors.high,
        linewidth = 3
    )


    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end
