include("../../src/main.jl")

let
    fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax5 = Axis(fig[1, 2], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax9 = Axis(fig[1, 3], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))

    ax2 = Axis(fig[2, 1], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax6 = Axis(fig[2, 2], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax10 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))

    ax3 = Axis(fig[3, 1], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax7 = Axis(fig[3, 2], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax11 = Axis(fig[3, 3], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))

    ax4 = Axis(fig[4, 1], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax8 = Axis(fig[4, 2], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))
    ax12 = Axis(fig[4, 3], xlabel = L"t/t_M", ylabel = L"R", limits = (0, 1000, -14, 14))

    files = [
        "data/Single_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem50.0TM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem1.0TM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem50.0TM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem1.0TM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem50.0TM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem1.0TM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
        "data/Single_Thermal/Single_R0[10]_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
    ]
    ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    txt =
        ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"]
    colors = repeat([my_blue, my_red, my_green, my_violet], 4)
    for n = 1:length(files)
        data = load_object(files[n])
        M = data.M                      # Mass of the mobile particles
        K_M = data.K_M                  # Trap force constant
        ΩM = √(K_M / M)                 # Trap frequency
        t_M = 2 * π / ΩM                # Period of the trapped mass
        idx = findall(x -> x < 1000, data.ts ./ t_M)
        lines!(
            ax[n],
            data.ts[idx] ./ t_M,
            data.Rs[idx] |> vec,
            color = colors[n],
            linewidth = 2,
        )
        text!(
            ax[n],
            txt[n],
            position = (5.5 * 1000 / 6, 20 / 2),
            textsize = 16,
            font = "CMU Serif",
        )
    end

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    save("Thermal_Example.pdf", fig)
end