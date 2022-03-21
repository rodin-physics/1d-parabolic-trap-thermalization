include("../../src/main.jl")

fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 32)
ax1 = fig[1, 1] = Axis(fig, xlabel = L"E/\hbar\Omega_T", ylabel = L"\ln[P(E)]")
ax2 = fig[1, 2] = Axis(fig, xlabel = L"E/\hbar\Omega_T", ylabel = L"\ln[P(E)]")
ax3 = fig[2, 1] = Axis(fig, xlabel = L"E/\hbar\Omega_T", ylabel = L"\ln[P(E)]")
ax4 = fig[2, 2] = Axis(fig, xlabel = L"E/\hbar\Omega_T", ylabel = L"\ln[P(E)]")

colors = [my_red, my_orange, my_green, my_blue, my_violet, colorant"rgba(0, 0, 0, 0.35)"]
labs = [L"Ω_T = 50", L"Ω_T = 100", L"Ω_T = 250", L"Ω_T = 500", L"Ω_T = 1000"]


function mkFigure(ax, filename, clr, lab, drop)
    data = load_object(filename)

    ΩT = data.ΩT
    F = data.F
    s = data.s
    Rs = data.Rs
    rs = data.rs
    ts = data.ts
    δ = ts[2] - ts[1]
    K_M = data.K_M
    M = data.M
    ħ = data.ħ

    Rs = Rs[drop:end, :]
    rs = rs[drop:end]

    prt = size(Rs)[2]

    kin_en = vcat(
        zeros(1, prt),
        ((Rs[2:end, :] - Rs[1:(end-1), :]) ./ δ) .^ 2 ./ 2 .* M ./ (ħ * ΩT),
    )
    pot_en = (Rs .^ 2 / 2 * K_M + F .* exp.(-(rs .- Rs) .^ 2 ./ (2 * s^2))) ./ (ħ * ΩT)
    tot_en = (pot_en + kin_en) |> vec

    hist_fit = fit(Histogram, tot_en, -10:0.05:25)
    hist_fit = normalize(hist_fit, mode = :pdf)
    scatter!(
        ax,
        (hist_fit.edges[1])[1:end-1],
        log.(hist_fit.weights),
        color = clr,
        label = lab,
        markersize = 12,
    )

end

files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT1000.0_τ1000.jld2",
]
for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax1, filename, clr, lab, 100000)
end
lines!(ax1, collect(0:0.1:7), -collect(0:0.1:7), linewidth = 2, color = :black)
axislegend(ax1, position = :lb, labelsize = 32)


xlims!(ax1, (-1, 10))
# ylims!(ax1, (-6, 1))

files = [
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F1_m1.0_d60_ΩT1000.0_τ1000.jld2",
]
for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax2, filename, clr, lab, 100000)
end
lines!(ax2, collect(0:0.1:7), -collect(0:0.1:7), linewidth = 2, color = :black)
axislegend(ax2, position = :lb, labelsize = 32)


xlims!(ax2, (-1, 10))
# ylims!(ax1, (-6, 1))


files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT1000.0_τ1000.jld2",
]
for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax3, filename, clr, lab, 100000)
end
lines!(ax3, collect(0:0.1:7), -collect(0:0.1:7), linewidth = 2, color = :black)
axislegend(ax3, position = :lb, labelsize = 32)


xlims!(ax3, (-5, 10))
# ylims!(ax1, (-6, 1))

files = [
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩT1000.0_τ1000.jld2",
]
for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax4, filename, clr, lab, 100000)
end
lines!(ax4, collect(0:0.1:7), -collect(0:0.1:7), linewidth = 2, color = :black)
axislegend(ax4, position = :lb, labelsize = 32)


xlims!(ax4, (-5, 10))
# ylims!(ax1, (-6, 1))

save("Energy_distribution.pdf", fig)