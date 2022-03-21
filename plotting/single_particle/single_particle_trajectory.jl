include("../../src/main.jl")
## Single Trajectory Analysis
data = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)

m = data.m                      # Mass of the chain atoms
k = data.k                      # Spring force constant
K = data.K                      # Confining potential force constant
Ωmin = sqrt(K / m)              # Smallest chain frequency
Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

M = data.M                      # Mass of the mobile particles
K_M = data.K_M                  # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass
Rs = data.Rs |> vec
rs = data.rs
ts = data.ts ./ t_M
step_size = 100
let
    fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"t/t_M", ylabel = L"R")
    ax3 = Axis(fig[1, 2], xlabel = L"t/t_M", ylabel = L"R")
    ax5 = Axis(fig[1, 3], xlabel = L"t/t_M", ylabel = L"R")

    ax2 = Axis(fig[2, 1], xlabel = L"t/t_M", ylabel = L"r")
    ax4 = Axis(fig[2, 2], xlabel = L"t/t_M", ylabel = L"r")
    ax6 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"r")

    ax_x_lims = [(0, 120), (0, 120), (88, 92), (88, 92), (116, 117), (116, 117)]
    ax_y_lims =
        [(-11, 11), (-0.15, 0.15), (-5.5, 5.5), (-0.15, 0.15), (-0.04, 0.04), (-0.04, 0.04)]
    ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    txt = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

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
    idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), ts)
    lines!(ax[1], ts[idx][1:step_size:end], Rs[idx][1:step_size:end], color = my_red, linewidth = 1)
    lines!(ax[2], ts[idx][1:step_size:end], rs[idx][1:step_size:end], color = my_red, linewidth = 1)

    idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), ts)
    lines!(ax[3], ts[idx], Rs[idx], color = my_red, linewidth = 1)
    lines!(ax[4], ts[idx], rs[idx], color = my_red, linewidth = 1)


    idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), ts)
    lines!(ax[5], ts[idx], Rs[idx], color = my_red, linewidth = 1)
    lines!(ax[6], ts[idx], rs[idx], color = my_red, linewidth = 1)

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    save("General_Example.pdf", fig)
end
