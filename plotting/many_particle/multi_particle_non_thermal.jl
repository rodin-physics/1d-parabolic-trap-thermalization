include("../../src/main.jl")

att = load_object(
    "data/Multi_Non_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ1000.jld2",
)
rep = load_object(
    "data/Multi_Non_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωTnothing_τ1000.jld2",
)
σsAtt = att.σs
ρsAtt = att.ρs
τsAtt = att.τs

σsRep = rep.σs
ρsRep = rep.ρs
τsRep = rep.τs



print(σsAtt[200000, :])
print(σsRep[200000, :])

step_size = 10
let
    fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = "Displacement")
    ax3 = Axis(fig[1, 2], xlabel = L"\tau")
    ax5 = Axis(fig[1, 3], xlabel = L"\tau")

    ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = "Displacement")
    ax4 = Axis(fig[2, 2], xlabel = L"\tau")
    ax6 = Axis(fig[2, 3], xlabel = L"\tau")

    ax_x_lims = [(0, 1000), (0, 1000), (250, 255), (250, 255), (600, 601), (600, 601)]
    # ax_y_lims = [(-150, 150), (-150, 150), (-35, 35), (-35, 35), (-5, 0), (-5, 0)]
    ax_y_lims = [(-150, 150), (-150, 150), (-150, 150), (-150, 150), (-15, 15), (-15, 15)]
    ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    txt = ["(a)", "(d)", "(b)", "(e)", "(c)", "(f)"]

    for n = 1:length(ax)
        xlims!(ax[n], ax_x_lims[n])
        ylims!(ax[n], ax_y_lims[n])
        text!(
            ax[n],
            txt[n],
            position = (
                ax_x_lims[n][1] + 5 / 6 * (ax_x_lims[n][2] - ax_x_lims[n][1]),
                ax_y_lims[n][2] / 2,
            ),
            textsize = 16,
            font = "CMU Serif",
        )
    end
    for ii = 1:size(σsAtt)[2]
        idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), τsAtt)
        lines!(
            ax[1],
            τsAtt[idx][1:step_size*10:end],
            σsAtt[:, ii][idx][1:step_size*10:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )
        lines!(
            ax[2],
            τsRep[idx][1:step_size*10:end],
            σsRep[:, ii][idx][1:step_size*10:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )

        idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), τsAtt)
        lines!(
            ax[3],
            τsAtt[idx][1:step_size:end],
            σsAtt[:, ii][idx][1:step_size:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )
        lines!(
            ax[4],
            τsRep[idx][1:step_size:end],
            σsRep[:, ii][idx][1:step_size:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )


        idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), τsAtt)
        lines!(
            ax[5],
            τsAtt[idx][1:step_size:end],
            σsAtt[:, ii][idx][1:step_size:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )
        lines!(
            ax[6],
            τsRep[idx][1:step_size:end],
            σsRep[:, ii][idx][1:step_size:end],
            color = (ii == 22  ? my_green : (ii == 2 ? my_blue : my_black)),
            linewidth = 1,
        )
    end
    idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), τsAtt)
    lines!(
        ax[1],
        τsAtt[idx][1:step_size*10:end],
        ρsAtt[idx][1:step_size*10:end],
        color = my_vermillion,
        linewidth = 1,
    )
    lines!(
        ax[2],
        τsRep[idx][1:step_size*10:end],
        ρsRep[idx][1:step_size*10:end],
        color = my_vermillion,
        linewidth = 1,
    )

    idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), τsAtt)
    lines!(
        ax[3],
        τsAtt[idx][1:step_size:end],
        ρsAtt[idx][1:step_size:end],
        color = my_vermillion,
        linewidth = 2,
    )
    lines!(
        ax[4],
        τsRep[idx][1:step_size:end],
        ρsRep[idx][1:step_size:end],
        color = my_vermillion,
        linewidth = 2,
    )

    idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), τsAtt)
    lines!(
        ax[5],
        τsAtt[idx][1:1:end],
        ρsAtt[idx][1:1:end],
        color = my_vermillion,
        linewidth = 2,
    )
    lines!(
        ax[6],
        τsRep[idx][1:1:end],
        ρsRep[idx][1:1:end],
        color = my_vermillion,
        linewidth = 2,
    )

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    save("Multi_NonThermal.pdf", fig)
end


# save("test_A.pdf",fig)
