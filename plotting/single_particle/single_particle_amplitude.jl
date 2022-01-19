include("../../src/main.jl")

function amp_idx(R)
    return filter(x -> (R[x] > R[x-1]) && (R[x] > R[x+1]), 2:(length(R)-1))
end

fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

ax1 = Axis(fig[1, 1], xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")
ax2 = Axis(fig[1, 2], xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")
ax3 = Axis(fig[1, 3], xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")

ax4 = Axis(fig[2, 1], xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")
ax5 = Axis(fig[2, 2], xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")

ax6 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"-dU/dR")
ax7 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"r", yaxisposition = :right)

ts_fit = exp.(range(-1, 6, length = 20))

colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)

# Mass dependence
data = [
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.0625_F-1_m0.25_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.0625_F-1_m0.5_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.0625_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
]
labs = [L"m = 1/4", L"m = 1/2", L"m = 1"]
colors = [my_red, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    m = d.m
    M = d.M          # Mass of the mobile particles
    K_M = d.K_M      # Trap force constant
    ΩM = √(K_M / M)  # Trap frequency
    t_M = 2 * π / ΩM # Period of the trapped mass
    ts = d.ts ./ t_M
    Rs = abs.(d.Rs |> vec)
    idx = amp_idx(Rs)
    scatter!(
        ax1,
        log.(ts[idx]),
        log.(abs.(Rs[idx] .^ 6 .- Rs[1] .^ 6)),
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax1,
        log.(ts_fit),
        log.(900 * (1 / d.m)^2 * ts_fit),
        color = colors[n],
        linewidth = 2,
    )

end

## Width Scaling

data = [
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.0625_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.125_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.5_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
]
labs = [L"s = 1/16", L"s = 1/8", L"s = 1/4", L"s = 1/2"]
colors = [my_red, my_green, my_blue, my_violet]
for n = 1:length(data)
    d = data[n]
    M = d.M          # Mass of the mobile particles
    K_M = d.K_M      # Trap force constant
    ΩM = √(K_M / M)  # Trap frequency
    t_M = 2 * π / ΩM # Period of the trapped mass
    ts = d.ts ./ t_M
    Rs = abs.(d.Rs |> vec)
    idx = amp_idx(Rs)
    scatter!(
        ax2,
        log.(ts[idx]),
        log.(abs.(Rs[idx] .^ 6 .- Rs[1] .^ 6)),
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax2,
        log.(ts_fit),
        log.(225000 * (d.s)^2 * ts_fit),
        color = colors[n],
        linewidth = 2,
    )

end

data = [
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.0625_F1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.125_F1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.5_F1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
]
labs = [L"s = 1/16", L"s = 1/8", L"s = 1/4", L"s = 1/2"]
colors = [my_red, my_green, my_blue, my_violet]
for n = 1:length(data)
    d = data[n]
    M = d.M          # Mass of the mobile particles
    K_M = d.K_M      # Trap force constant
    ΩM = √(K_M / M)  # Trap frequency
    t_M = 2 * π / ΩM # Period of the trapped mass
    ts = d.ts ./ t_M
    Rs = abs.(d.Rs |> vec)
    idx = amp_idx(Rs)
    scatter!(
        ax3,
        log.(ts[idx]),
        log.(abs.(Rs[idx] .^ 6 .- Rs[1] .^ 6)),
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax3,
        log.(ts_fit),
        log.(450000 * (d.s)^2 * ts_fit),
        color = colors[n],
        linewidth = 2,
    )
end

## Depth Scaling
data = [
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-0.5_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-2.0_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-4.0_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
]
labs = [L"F = -1/2", L"F = -1", L"F = -2", L"F = -4"]
colors = [my_red, my_green, my_blue, my_violet]
for n = 1:length(data)
    d = data[n]
    M = d.M          # Mass of the mobile particles
    K_M = d.K_M      # Trap force constant
    ΩM = √(K_M / M)  # Trap frequency
    t_M = 2 * π / ΩM # Period of the trapped mass
    ts = d.ts ./ t_M
    Rs = abs.(d.Rs |> vec)
    idx = amp_idx(Rs)
    scatter!(
        ax4,
        log.(ts[idx]),
        log.(abs.(Rs[idx] .^ 6 .- Rs[1] .^ 6)),
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax4,
        log.(ts_fit),
        log.(15000 * (d.F)^2 * ts_fit),
        color = colors[n],
        linewidth = 2,
    )

end

data = [
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F0.5_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F1_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F2.0_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F4.0_m1.0_d60_ΩTnothing_τ250.jld2",
    ),
]
labs = [L"F = 1/2", L"F = 1", L"F = 2", L"F = 4"]
colors = [my_red, my_green, my_blue, my_violet]
for n = 1:length(data)
    d = data[n]
    M = d.M          # Mass of the mobile particles
    K_M = d.K_M      # Trap force constant
    ΩM = √(K_M / M)  # Trap frequency
    t_M = 2 * π / ΩM # Period of the trapped mass
    ts = d.ts ./ t_M
    Rs = abs.(d.Rs |> vec)
    idx = amp_idx(Rs)
    scatter!(
        ax5,
        log.(ts[idx]),
        log.(abs.(Rs[idx] .^ 6 .- Rs[1] .^ 6)),
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax5,
        log.(ts_fit),
        log.(18000 * (d.F)^2 * ts_fit),
        color = colors[n],
        linewidth = 2,
    )
end

txt = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
axs = [ax1, ax2, ax3, ax4, ax5]
for n = 1:length(axs)
    ax = axs[n]
    xlims!(ax, (-1, 6))
    ylims!(ax, (6, 14))
    text!(
        ax,
        txt[n],
        position = (-1 + 1 / 10 * 7, 6 + 5 / 6 * 8),
        textsize = 16,
        font = "CMU Serif",
    )
    ax.yticks = 6:2:14
    axislegend(ax, position = :rb)

end

## Attractive vs repulsive
repulsive = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F1_m1.0_d60_ΩTnothing_τ250.jld2",
)
attractive = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
m = attractive.m                # Mass of the chain atoms
k = attractive.k                # Spring force constant
K = attractive.K                # Confining potential force constant
Ωmin = sqrt(K / m)              # Smallest chain frequency
Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

M = attractive.M                # Mass of the mobile particles
K_M = attractive.K_M            # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass

Rs_attractive = attractive.Rs |> vec
Rs_repulsive = repulsive.Rs |> vec

rs_attractive = attractive.rs |> vec
rs_repulsive = repulsive.rs |> vec

Us_attractive = attractive.U_pr_mob |> vec
Us_repulsive = repulsive.U_pr_mob |> vec

n_pts = length(Rs_attractive)
ts = attractive.ts ./ t_M
idx = findall(x -> (x >= 0) && (x <= 1), ts)

lines!(ax6, ts[idx], -Us_attractive[idx], color = my_red, linewidth = 2, label = L"F = -1")
lines!(ax6, ts[idx], -Us_repulsive[idx], color = my_blue, linewidth = 2, label = L"F = 1")
lines!(ax7, ts[idx], rs_attractive[idx], color = my_red, linewidth = 4, linestyle = "-")
lines!(ax7, ts[idx], rs_repulsive[idx], color = my_blue, linewidth = 4, linestyle = "-")

xlims!(ax6, (0.23, 0.27))
xlims!(ax7, (0.23, 0.27))
ylims!(ax6, (-2.5, 2.5))
ylims!(ax7, (-0.03, 0.03))
l = axislegend(ax6, position = :lb)
text!(
    ax6,
    "(f)",
    position = (0.23 + 1 / 10 * 0.04, -2.5 + 5 / 6 * 5),
    textsize = 16,
    font = "CMU Serif",
)

save("Amplitude_Decay.pdf", fig)
