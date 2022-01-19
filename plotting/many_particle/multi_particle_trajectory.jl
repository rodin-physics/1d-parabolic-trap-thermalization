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
files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT10.0_τ1000.jld2",
]

data = [load_object(f) for f in files]
M = data[1].M                   # Mass of the mobile particles
K_M = data[1].K_M               # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass

times = [d.ts for d in data] ./ t_M
Rs = [d.Rs for d in data]
R_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in Rs]

fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 32)
ax1 = fig[1, 1] = Axis(fig, xlabel = L"t/t_M", ylabel = L"R_\mathrm{rms}")
ax2 = fig[2, 1] = Axis(fig, xlabel = L"t/t_M", ylabel = L"R_\mathrm{rms}")

for ii = 1:length(files)
    lines!(
        ax1,
        times[ii],
        R_rms[ii] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax1, position = :lt, labelsize = 32, orientation = :horizontal)

files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
]

data = [load_object(f) for f in files]
M = data[1].M                   # Mass of the mobile particles
K_M = data[1].K_M               # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass

times = [d.ts for d in data] ./ t_M
Rs = [d.Rs for d in data]
R_rms = [sqrt.(mean(x .^ 2, dims = 2)) for x in Rs]
for ii = 1:length(files)
    lines!(
        ax2,
        times[ii],
        R_rms[ii] |> vec,
        color = colors[ii],
        label = labs[ii],
        linewidth = 2,
    )
end
axislegend(ax2, position = :lt, labelsize = 32, orientation = :horizontal)

save("Cloud_Size.pdf", fig)

