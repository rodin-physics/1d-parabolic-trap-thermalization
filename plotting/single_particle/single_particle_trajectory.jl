include("../../src/main.jl")
## Single Trajectory Analysis
data = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
σs = data.σs |> vec
ρs = data.ρs
τs = data.τs
step_size = 10
let
    fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"τ", ylabel = L"σ")
    ax3 = Axis(fig[1, 2], xlabel = L"τ", ylabel = L"σ")
    ax5 = Axis(fig[1, 3], xlabel = L"τ", ylabel = L"σ")

    ax2 = Axis(fig[2, 1], xlabel = L"τ", ylabel = L"ρ")
    ax4 = Axis(fig[2, 2], xlabel = L"τ", ylabel = L"ρ")
    ax6 = Axis(fig[2, 3], xlabel = L"τ", ylabel = L"ρ")

    ax_x_lims = [(0, 250), (0, 250), (170, 175), (170, 175), (230, 231), (230, 231)]
    ax_y_lims =
        [(-110, 110), (-1.5, 1.5), (-35, 35), (-1.5, 1.5), (-0.2, 0.2), (-0.15, 0.15)]
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
    idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), τs)
    lines!(ax[1], τs[idx][1:step_size:end], σs[idx][1:step_size:end], color = my_blue, linewidth = 1)
    lines!(ax[2], τs[idx][1:step_size:end], ρs[idx][1:step_size:end], color = my_vermillion, linewidth = 1)

    idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), τs)
    lines!(ax[3], τs[idx], σs[idx], color = my_blue, linewidth = 1)
    lines!(ax[4], τs[idx], ρs[idx], color = my_vermillion, linewidth = 1)


    idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), τs)
    lines!(ax[5], τs[idx], σs[idx], color = my_blue, linewidth = 1)
    lines!(ax[6], τs[idx], ρs[idx], color = my_vermillion, linewidth = 1)

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    # fig
    save("General_Example.pdf", fig)
end
