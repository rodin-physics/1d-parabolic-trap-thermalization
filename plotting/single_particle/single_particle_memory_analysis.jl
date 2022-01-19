include("../../src/main.jl")
## Different Starts
data = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)

m = data.m                      # Mass of the chain atoms
k = data.k                      # Spring force constant
K = data.K                      # Confining potential force constant
Ωmin = sqrt(K / m)              # Smallest chain frequency
Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

M = data.M                      # Mass of the mobile particles
K_M = data.K_M                  # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass
Rs10 = data.Rs |> vec
rs10 = data.rs
ts10 = data.ts ./ t_M

data = load_object(
    "data/Single_Non_Thermal/Single_R0[7.422582732835284]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
Rs7_5 = data.Rs |> vec
rs7_5 = data.rs
ts7_5 = data.ts ./ t_M

let
    fig = Figure(resolution = (1200, 300), font = "CMU Serif", fontsize = 14)

    ax1 = fig[1, 1] = Axis(fig, xlabel = L"t/t_M", ylabel = L"R")
    ax2 = fig[1, 2] = Axis(fig, xlabel = L"t/t_M", ylabel = L"R")
    ax3 = fig[1, 3] = Axis(fig, xlabel = L"t/t_M", ylabel = L"R")

    ln1 = lines!(ax1, ts10, Rs10, color = my_blue, linewidth = 2)
    ln2 = lines!(ax1, ts7_5 .+ ts7_5[84317], Rs7_5, color = my_red, linewidth = 2)
    xlims!(ax1, (0, 120))
    ylims!(ax1, (-11, 11))
    text!(ax1, "(a)", position = (5 * 120 / 6, 11 / 2), textsize = 16, font = "CMU Serif")

    idx = findall(x -> (x <= 74) && (x >= 69), ts10)
    ln3 = lines!(ax2, ts10[idx], Rs10[idx], color = my_blue, linewidth = 2)
    idx = findall(x -> (x <= 74) && (x >= 69), ts7_5 .+ ts7_5[84317])
    ln4 = lines!(ax2, ts7_5[idx] .+ ts7_5[84317], Rs7_5[idx], color = my_red, linewidth = 2)
    xlims!(ax2, (69, 74))
    ylims!(ax2, (-10, 10))
    text!(
        ax2,
        "(b)",
        position = (69 + 5 * 5 / 6, 10 / 2),
        textsize = 16,
        font = "CMU Serif",
    )

    idx = findall(x -> (x <= 111) && (x >= 110), ts10)
    ln5 = lines!(ax3, ts10[idx], Rs10[idx], color = my_blue, linewidth = 2)
    idx = findall(x -> (x <= 111) && (x >= 110), ts7_5 .+ ts7_5[84318])
    ln6 = lines!(ax3, ts7_5[idx] .+ ts7_5[84317], Rs7_5[idx], color = my_red, linewidth = 2)

    xlims!(ax3, (110, 111))
    ylims!(ax3, (-0.04, 0.04))
    text!(ax3, "(c)", position = (110 + 5 / 6, 0.04 / 2), textsize = 16, font = "CMU Serif")

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    ax3.xticks = 110:1:111
    # fig

    save("Different_Starts.pdf", fig)
end

## Memory length
Full = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_MemInfTM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T50 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem50.0TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T10 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem10.0TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T2 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem2.0TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T155 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem1.55TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T105 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem1.05TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T101 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem1.01TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T1 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem1.0TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T055 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.55TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T051 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.51TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T05 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.5TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T025 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.25TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T01 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.1TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)
T005 = load_object(
    "data/Single_Non_Thermal/Single_R0[10]_Mem0.05TM_s0.25_F-1_m1.0_d60_ΩTnothing_τ250.jld2",
)

m = Full.m                      # Mass of the chain atoms
k = Full.k                      # Spring force constant
K = Full.K                      # Confining potential force constant
Ωmin = sqrt(K / m)              # Smallest chain frequency
Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

M = Full.M                      # Mass of the mobile particles
K_M = Full.K_M                  # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass
ts = Full.ts ./ t_M

function amp_idx(R)
    return filter(
        x -> (abs(R[x]) > abs(R[x-1])) && (abs(R[x]) > abs(R[x+1]) && abs(R[x]) > Full.s),
        2:(length(R)-1),
    )
end

fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 32)
ax1 = fig[1, 1] = Axis(fig, xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")
ax2 = fig[2, 1] = Axis(fig, xlabel = L"\ln(t/t_M)", ylabel = L"\ln\left(R_0^6 - R^6\right)")

Rs = [(x.Rs |> vec) for x in [Full, T50, T10, T2, T025, T005]]
idx = [amp_idx(x) for x in Rs]
markers = [:circle, :+, :circle, :rect, :+, :circle]
marker_size = [10, 20, 10, 10, 20, 10] .* 2
colors = [:black, my_red, my_orange, my_green, my_blue, my_violet]
labs = [L"\infty", L"50t_M", L"10t_M", L"2t_M", L"t_M/4", L"t_M/20"]

for ii = 2:length(idx)
    p = Rs[ii]
    ind = idx[ii]
    scatter!(
        ax1,
        log.(ts[ind]),
        log.(abs.(p[1] .^ 6 .- abs.(p[ind]) .^ 6)),
        marker = markers[ii],
        markersize = marker_size[ii],
        color = colors[ii],
        label = labs[ii],
    )
end

pos_Inf = Rs[1]
indices_Inf = idx[1]
lines!(
    ax1,
    log.(ts[indices_Inf]),
    log.(abs.(pos_Inf[1] .^ 6 .- abs.(pos_Inf[indices_Inf]) .^ 6)),
    linewidth = 2,
    color = :black,
)
axislegend(ax1, position = :rb, labelsize = 32)

text!(ax1, "(a)", position = (-1 / 2, 13.5), textsize = 32, font = "CMU Serif")
xlims!(ax1, (-1, 5))
ax1.yticks = 8:1:14

# Second panel


Rs = [x.Rs |> vec for x in [Full, T155, T101, T1, T055, T05]]
idx = [amp_idx(x) for x in Rs]

markers = [:circle, :+, :circle, :rect, :+, :circle]
marker_size = [10, 20, 10, 10, 20, 10] .* 2
colors = [:black, my_red, my_orange, my_green, my_blue, my_violet]
labs = [L"\infty", L"1.55t_M", L"1.01t_M", L"1t_M", L"0.55t_M", L"0.5t_M"]

for ii = 2:length(idx)
    p = Rs[ii]
    ind = idx[ii]
    scatter!(
        ax2,
        log.(ts[ind]),
        log.((p[1] .^ 6 .- abs.(p[ind]) .^ 6)),
        marker = markers[ii],
        markersize = marker_size[ii],
        color = colors[ii],
        label = labs[ii],
    )
end

pos_Inf = Rs[1]
indices_Inf = idx[1]
lines!(
    ax2,
    log.(ts[indices_Inf]),
    log.(abs.(pos_Inf[1] .^ 6 .- abs.(pos_Inf[indices_Inf]) .^ 6)),
    linewidth = 2,
    color = :black,
)
axislegend(ax2, position = :rb, labelsize = 32)

text!(ax2, "(b)", position = (-1 / 2, 13.5), textsize = 32, font = "CMU Serif")

xlims!(ax2, (-1, 5))
ax2.yticks = 8:1:14

save("Memory_Amp.pdf", fig)
