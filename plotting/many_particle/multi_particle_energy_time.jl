include("../../src/main.jl")

files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT10.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT1000.0_τ1000.jld2",
]

colors = [colorant"rgba(0, 0, 0, 0.35)", my_red, my_orange, my_green, my_blue, my_violet]
labs = [L"Ω_T = 10", L"Ω_T = 50", L"Ω_T = 100", L"Ω_T = 250", L"Ω_T = 500", L"Ω_T = 1000"]

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 24)
ax1 = Axis(fig[1, 1], ylabel = L"E/\hbar\Omega_T", xlabel = L"t/t_M")

for ii = 1:length(files)

    data = load_object(files[ii])

    ΩT = data.ΩT
    F = data.F
    s = data.s
    Rs = data.Rs
    rs = data.rs
    k = data.k
    K = data.K
    ts = data.ts
    δ = ts[2] - ts[1]
    K_M = data.K_M
    M = data.M
    m = data.m
    ħ = data.ħ
    ΩM = √(K_M / M)                 # Trap frequency
    t_M = 2 * π / ΩM                # Period of the trapped mass
    Ωmin = sqrt(K / m)                  # Smallest chain frequency
    Ωmax = sqrt(4 * k / m + K / m)      # Largest chain frequency

    prt = size(Rs)[2]

    kin_en = vcat(
        zeros(1, prt),
        ((Rs[2:end, :] - Rs[1:(end-1), :]) ./ δ) .^ 2 ./ 2 .* M ./ (ħ * ΩT),
    )
    pot_en = (Rs .^ 2 / 2 * K_M + F .* exp.(-(rs .- Rs) .^ 2 ./ (2 * s^2))) ./ (ħ * ΩT)
    tot_en = (pot_en + kin_en)
    tot_en = sum(tot_en, dims = 2) |> vec
    lines!(
        ax1,
        ts ./ t_M,
        tot_en ./ prt,
        color = colors[ii],
        linewidth = 2,
        label = labs[ii],
    )


end

ylims!(ax1, (-1 / 4, 4))
axislegend(ax1, position = :rt, labelsize = 24)

fig

save("Energy_over_time.pdf", fig)
