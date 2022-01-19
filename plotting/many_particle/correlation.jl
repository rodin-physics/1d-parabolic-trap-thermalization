include("../../src/main.jl")

colors = [my_violet, my_blue, my_green, my_orange, my_red, colorant"rgba(0, 0, 0, 0.35)"]
labs = [
    L"$Ω_T = 1000$",
    L"$Ω_T = 500$",
    L"$Ω_T = 250$",
    L"$Ω_T = 100$",
    L"$Ω_T = 50$",
    L"$Ω_T = 10$",
]

fig = Figure(resolution = (600, 1600), font = "CMU Serif", fontsize = 16)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle R_i(t)R_i(t+\tau)\rangle_\mathrm{avg}",
)
ax2 = Axis(
    fig[2, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle R_i(t)R_j(t+\tau)\rangle_\mathrm{avg}",
)
ax3 = Axis(
    fig[3, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle R_i(t)R_i(t+\tau)\rangle_\mathrm{avg}",
)
ax4 = Axis(
    fig[4, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle R_i(t)R_j(t+\tau)\rangle_\mathrm{avg}",
)

files = [
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT1000.0_hbar0.005.jld2",
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT500.0_hbar0.005.jld2",
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT250.0_hbar0.005.jld2",
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT100.0_hbar0.005.jld2",
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT50.0_hbar0.005.jld2",
    "data/Correlations/Corr_F1.0_K1.0_k20.0_m1.0_ΩT10.0_hbar0.005.jld2",
]

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
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT1000.0_hbar0.005.jld2",
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT500.0_hbar0.005.jld2",
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT250.0_hbar0.005.jld2",
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT100.0_hbar0.005.jld2",
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT50.0_hbar0.005.jld2",
    "data/Correlations/Corr_F-1.0_K1.0_k20.0_m1.0_ΩT10.0_hbar0.005.jld2",
]

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

save("Correlation.pdf", fig)


colors = [my_violet, my_blue, my_green, my_orange, my_red, colorant"rgba(0, 0, 0, 0.35)"]
labs = [
    L"$Ω_T = 1000$",
    L"$Ω_T = 500$",
    L"$Ω_T = 250$",
    L"$Ω_T = 100$",
    L"$Ω_T = 50$",
    L"$Ω_T = 10$",
]

fig = Figure(resolution = (600, 1600), font = "CMU Serif", fontsize = 16)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle E_i(t)E_i(t+\tau)\rangle_\mathrm{avg}",
)
ax2 = Axis(
    fig[2, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle E_i(t)E_j(t+\tau)\rangle_\mathrm{avg}",
)
ax3 = Axis(
    fig[3, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle E_i(t)E_i(t+\tau)\rangle_\mathrm{avg}",
)
ax4 = Axis(
    fig[4, 1],
    xlabel = L"τ/t_M",
    ylabel = L"\langle E_i(t)E_j(t+\tau)\rangle_\mathrm{avg}",
)

files = [
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT1000.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT500.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT250.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT100.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT50.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F1.0_K1.0_k20.0_m1.0_ΩT10.0_hbar0.005.jld2",
]

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
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT1000.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT500.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT250.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT100.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT50.0_hbar0.005.jld2",
    "data/Correlations/Energy_Corr_F-1.0_K1.0_k20.0_m1.0_ΩT10.0_hbar0.005.jld2",
]

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

save("Energy_Correlation.pdf", fig)