include("../../src/main.jl")

let
    fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax5 = Axis(fig[1, 2], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax9 = Axis(fig[1, 3], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))

    ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax6 = Axis(fig[2, 2], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax10 = Axis(fig[2, 3], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))

    ax3 = Axis(fig[3, 1], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax7 = Axis(fig[3, 2], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax11 = Axis(fig[3, 3], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))

    ax4 = Axis(fig[4, 1], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax8 = Axis(fig[4, 2], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    ax12 = Axis(fig[4, 3], xlabel = L"\tau", ylabel = L"\sigma", limits = (0, 1000, -140, 140))
    step_size = 100
    files = [
        "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT1.0e-5_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT1.0e-5_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT1.0e-5_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT1.0e-5_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
        "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
    ]
    ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    txt =
        ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"]
    colors = repeat([my_vermillion, my_green, my_sky, my_blue], 4)
    for n = 1:length(files)
        data = load_object(files[n])
        lines!(
            ax[n],
            data.τs[1:step_size:end] ,
            data.σs[1:step_size:end] |> vec,
            color = colors[n],
            linewidth = 0.1,
        )
        # text!(
        #     ax[n],
        #     txt[n],
        #     position = (5.5 * 1000 / 6, 20 / 2),
        #     textsize = 16,
        #     font = "CMU Serif",
        # )
    end

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    fig
    save("Thermal_Example.pdf", fig)
end