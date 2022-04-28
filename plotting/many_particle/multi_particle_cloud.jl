include("../../src/main.jl")

fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 32)
ax1 = fig[1, 1] = Axis(fig, xlabel = L"\tau", ylabel = L"\sigma_\mathrm{rms}")
ax2 = fig[2, 1] = Axis(fig, xlabel = L"\tau", ylabel = L"\sigma_\mathrm{rms}")

files = [
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]
colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]
step_size = 250
labs = [
    L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

data = [load_object(f) for f in files]

times = [d.τs for d in data]
σs = [d.σs for d in data]
σ_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in σs]

for ii = 1:length(files)
    lines!(
        ax1,
        times[ii][1:step_size:end],
        σ_rms[ii][1:step_size:end] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax1, position = :lt, labelsize = 32, orientation = :horizontal)

files = [
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]
colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]
step_size = 250
labs = [
    L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

data = [load_object(f) for f in files]

times = [d.τs for d in data]
σs = [d.σs for d in data]
σ_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in σs]

for ii = 1:length(files)
    lines!(
        ax2,
        times[ii][1:step_size:end],
        σ_rms[ii][1:step_size:end] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax2, position = :lt, labelsize = 32, orientation = :horizontal)
ylims!(ax1, (-5, 80))
ylims!(ax2, (-5, 80))
save("Cloud_Size.pdf", fig)

# WARMING UP
fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 32)
ax1 = fig[1, 1] = Axis(fig, xlabel = L"\tau", ylabel = L"\sigma_\mathrm{rms}")
ax2 = fig[2, 1] = Axis(fig, xlabel = L"\tau", ylabel = L"\sigma_\mathrm{rms}")

files = [
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT50.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]
colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]
step_size = 250
labs = [
    L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

data = [load_object(f) for f in files]

times = [d.τs for d in data]
σs = [d.σs for d in data]
σ_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in σs]

for ii = 1:length(files)
    lines!(
        ax1,
        times[ii][1:step_size:end],
        σ_rms[ii][1:step_size:end] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax1, position = :lt, labelsize = 32, orientation = :horizontal)

files = [
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT50.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/WarmingUp_Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]
colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]
step_size = 250
labs = [
    L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

data = [load_object(f) for f in files]

times = [d.τs for d in data]
σs = [d.σs for d in data]
σ_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in σs]

for ii = 1:length(files)
    lines!(
        ax2,
        times[ii][1:step_size:end],
        σ_rms[ii][1:step_size:end] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax2, position = :lt, labelsize = 32, orientation = :horizontal)
ylims!(ax1, (-5, 80))
ylims!(ax2, (-5, 80))

save("WarmingUp_Cloud_Size.pdf", fig)