include("../../src/main.jl")

colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]|> reverse
labs = [
    # L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]|> reverse

fig = Figure(resolution = (600, 1600), font = "CMU Serif", fontsize = 16)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\tau",
    ylabel = L"\mathrm{corr}\left[\sigma_i(t)\sigma_i(t+\tau)\right]_\mathrm{avg}",
)
ax2 = Axis(
    fig[2, 1],
    xlabel = L"\tau",
    ylabel = L"\mathrm{corr}\left[\sigma_i(t)\sigma_j(t+\tau)\right]_\mathrm{avg}",
)
ax3 = Axis(
    fig[3, 1],
    xlabel = L"\tau",
    ylabel = L"\mathrm{corr}\left[\sigma_i(t)\sigma_i(t+\tau)\right]_\mathrm{avg}",
)
ax4 = Axis(
    fig[4, 1],
    xlabel = L"\tau",
    ylabel = L"\mathrm{corr}\left[\sigma_j(t)\sigma_i(t+\tau)\right]_\mathrm{avg}",
)

files = [
    "data/Correlations/Corr_Φ0500.0_λ4.0_μ2.0_ωT100.0.jld2",
    "data/Correlations/Corr_Φ0500.0_λ4.0_μ2.0_ωT250.0.jld2",
    "data/Correlations/Corr_Φ0500.0_λ4.0_μ2.0_ωT500.0.jld2",
    "data/Correlations/Corr_Φ0500.0_λ4.0_μ2.0_ωT1000.0.jld2",
    "data/Correlations/Corr_Φ0500.0_λ4.0_μ2.0_ωT2500.0.jld2",
]|> reverse

for ii = 1:length(files)
    data = load_object(files[ii])
    autocorr = map(x -> x |> diag |> mean, data[1])
    crosscorr = map(x -> (sum(x) - (x |> diag |> sum)) / 600, data[1])
    lines!(ax1, data[2], autocorr, color = colors[ii], label = labs[ii])
    lines!(ax2, data[2], crosscorr, color = colors[ii], label = labs[ii])
end
axislegend(ax1, position = :rt, labelsize = 14, nbanks = 2)
axislegend(ax2, position = :rt, labelsize = 14, nbanks = 2)
xlims!(ax1, (0, 100))
xlims!(ax2, (0, 100))


files = [
    "data/Correlations/Corr_Φ0-500.0_λ4.0_μ2.0_ωT100.0.jld2",
    "data/Correlations/Corr_Φ0-500.0_λ4.0_μ2.0_ωT250.0.jld2",
    "data/Correlations/Corr_Φ0-500.0_λ4.0_μ2.0_ωT500.0.jld2",
    "data/Correlations/Corr_Φ0-500.0_λ4.0_μ2.0_ωT1000.0.jld2",
    "data/Correlations/Corr_Φ0-500.0_λ4.0_μ2.0_ωT2500.0.jld2",
]|> reverse

for ii = 1:length(files)
    data = load_object(files[ii])
    autocorr = map(x -> x |> diag |> mean, data[1])
    crosscorr = map(x -> (sum(x) - (x |> diag |> sum)) / 600, data[1])
    lines!(ax3, data[2], autocorr, color = colors[ii], label = labs[ii])
    lines!(ax4, data[2], crosscorr, color = colors[ii], label = labs[ii])
end
axislegend(ax3, position = :rt, labelsize = 14, nbanks = 2)
axislegend(ax4, position = :rt, labelsize = 14, nbanks = 2)
xlims!(ax3, (0, 100))
xlims!(ax4, (0, 100))
fig


save("Correlation.pdf", fig)

