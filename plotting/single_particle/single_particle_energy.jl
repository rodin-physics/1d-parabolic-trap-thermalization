include("../../src/main.jl")

fig = Figure(resolution = (2400, 1200), font = "CMU Serif", fontsize = 32)
ax1 = Axis(fig[1, 1], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax2 = Axis(fig[2, 1], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax3 = Axis(fig[1, 2], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax4 = Axis(fig[2, 2], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax5 = Axis(fig[1, 3], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax6 = Axis(fig[2, 3], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax7 = Axis(fig[1, 4], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")
ax8 = Axis(fig[2, 4], xlabel = L"E/\omega_T", ylabel = L"\ln[P(E)]")

colors = [my_vermillion, my_orange, my_green, my_sky, my_blue, my_black]
labs = [
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

function mkFigure(ax, filename, clr, lab, drop)
    data = load_object(filename)

    ωT = data.ωT
    Φ0 = data.Φ
    λ = data.λ
    σs = data.σs
    ρs = data.ρs
    τs = data.τs
    δ = τs[2] - τs[1]

    σs = σs[drop:end, :]
    ρs = ρs[drop:end]

    prt = size(σs)[2]

    kin_en =
        vcat(zeros(1, prt), ((σs[2:end, :] - σs[1:(end-1), :]) ./ δ) .^ 2 ./ 2 ./ ωT) ./
        (2 * pi)^2
    pot_en = (σs .^ 2 / 2 + Φ0 .* exp.(-(ρs .- σs) .^ 2 ./ (2 * λ^2))) ./ ωT
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
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax1, filename, clr, lab, 100000)
end
lines!(ax1, collect(0:0.1:7), -collect(0:0.1:7), linewidth = 2, color = :black)
axislegend(ax1, position = :lb, labelsize = 32)

files = [
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]


for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax3, filename, clr, lab, 100000)
end

lines!(ax3, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax3, position = :lb, labelsize = 32)

files = [
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax5, filename, clr, lab, 100000)
end

lines!(ax5, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax5, position = :lb, labelsize = 32)

files = [
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax7, filename, clr, lab, 100000)
end

lines!(ax7, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax7, position = :lb, labelsize = 32)

files = [
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax2, filename, clr, lab, 150000)
end

lines!(ax2, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax2, position = :lb, labelsize = 32)


files = [
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax4, filename, clr, lab, 150000)
end

lines!(ax4, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax4, position = :lb, labelsize = 32)

files = [
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax6, filename, clr, lab, 150000)
end

lines!(ax6, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax6, position = :lb, labelsize = 32)


files = [
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Single_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for ii = 1:length(files)
    filename = files[ii]
    clr = colors[ii]
    lab = labs[ii]
    mkFigure(ax8, filename, clr, lab, 150000)
end

lines!(ax8, collect(0:0.1:6), -collect(0:0.1:6), linewidth = 2, color = :black)
axislegend(ax8, position = :lb, labelsize = 32)

fig
# save("Single_Energy2.pdf", fig)
