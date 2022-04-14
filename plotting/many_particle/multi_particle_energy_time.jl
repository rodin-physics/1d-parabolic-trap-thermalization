include("../../src/main.jl")

files = [
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]
colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]
step_size = 100
labs = [
    L"\omega_T = 50",
    L"\omega_T = 100",
    L"\omega_T = 250",
    L"\omega_T = 500",
    L"\omega_T = 1000",
    L"\omega_T = 2500",
]

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 24)
ax1 = Axis(fig[1, 1], ylabel = L"E/\omega_T", xlabel = L"\tau")

for ii = 1:length(files)

    data = load_object(files[ii])

    ωT = data.ωT
    Φ0 = data.Φ
    λ = data.λ
    σs = data.σs
    ρs = data.ρs
    τs = data.τs
    δ = τs[2] - τs[1]

    σs = σs
    ρs = ρs

    prt = size(σs)[2]

    kin_en =
        vcat(zeros(1, prt), ((σs[2:end, :] - σs[1:(end-1), :]) ./ δ) .^ 2 ./ 2 ./ ωT) ./
        (2 * pi)^2
    pot_en = (σs .^ 2 / 2 + Φ0 .* exp.(-(ρs .- σs) .^ 2 ./ (2 * λ^2))) ./ ωT
    tot_en = (pot_en + kin_en)

    tot_en = sum(tot_en, dims = 2) |> vec
    lines!(
        ax1,
        τs[1:step_size:end],
        tot_en[1:step_size:end] ./ prt,
        color = colors[ii],
        linewidth = 2,
        label = labs[ii],
    )

end

# ylims!(ax1, (-1 / 4, 4))
axislegend(ax1, position = :rt, labelsize = 24)

save("Energy_over_time.pdf", fig)
