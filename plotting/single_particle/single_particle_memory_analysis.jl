include("../../src/main.jl")
## Different Starts
data = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)

σs100 = data.σs |> vec
ρs100 = data.ρs
τs100 = data.τs

data = load_object(
    "data/Single_Non_Thermal/Single_σ0[73.4561560502402]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
σs73 = data.σs |> vec
ρs73 = data.ρs
τs73 = data.τs
step_size = 10
let
    fig = Figure(resolution = (1200, 300), font = "CMU Serif", fontsize = 14)

    ax1 = fig[1, 1] = Axis(fig, xlabel = L"\tau", ylabel = L"\rho")
    ax2 = fig[1, 2] = Axis(fig, xlabel = L"\tau", ylabel = L"\rho")
    ax3 = fig[1, 3] = Axis(fig, xlabel = L"\tau", ylabel = L"\rho")

    ln1 = lines!(
        ax1,
        τs100[1:step_size:end],
        σs100[1:step_size:end],
        color = my_blue,
        linewidth = 2,
    )
    ln2 = lines!(
        ax1,
        (τs73.+τs73[150869])[1:step_size:end],
        σs73[1:step_size:end],
        color = my_vermillion,
        linewidth = 2,
    )
    xlims!(ax1, (0, 250))
    ylims!(ax1, (-110, 110))
    text!(
        ax1,
        "(a)",
        position = (5 * 250 / 6, 2 * 110 / 3),
        textsize = 16,
        font = "CMU Serif",
    )

    idx = findall(x -> (x <= 130) && (x >= 125), τs100)
    ln3 = lines!(ax2, τs100[idx], σs100[idx], color = my_blue, linewidth = 2)
    idx = findall(x -> (x <= 130) && (x >= 125), τs73 .+ τs73[150870])
    ln4 = lines!(
        ax2,
        τs73[idx] .+ τs73[150870],
        σs73[idx],
        color = my_vermillion,
        linewidth = 2,
    )
    xlims!(ax2, (125, 130))
    ylims!(ax2, (-110, 110))
    text!(
        ax2,
        "(b)",
        position = (125 + 5 * 5 / 6, 2 * 110 / 3),
        textsize = 16,
        font = "CMU Serif",
    )

    idx = findall(x -> (x <= 181) && (x >= 180), τs100)
    ln5 = lines!(ax3, τs100[idx], σs100[idx], color = my_blue, linewidth = 2)
    idx = findall(x -> (x <= 181) && (x >= 180), τs73 .+ τs73[150870])
    ln6 = lines!(
        ax3,
        τs73[idx] .+ τs73[150870],
        σs73[idx],
        color = my_vermillion,
        linewidth = 2,
    )

    xlims!(ax3, (180, 181))
    ylims!(ax3, (-0.25, 0.25))
    text!(
        ax3,
        "(c)",
        position = (180 + 5 / 6, 2 * 0.25 / 3),
        textsize = 16,
        font = "CMU Serif",
    )

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    ax3.xticks = 180:1:181
    # fig

    save("Different_Starts.pdf", fig)
end

## Memory length
Full = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T50 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ050.0_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T2 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ02.0_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T155 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ01.55_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T105 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ01.05_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T101 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ01.01_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T1 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ01.0_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T055 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ00.55_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T051 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ00.51_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T05 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ00.5_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T01 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ00.1_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)
T005 = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ00.05_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)

function amp_idx(σ)
    return filter(x -> (σ[x] > σ[x-1]) && (σ[x] > σ[x+1]), 2:(length(σ)-1))
end

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 32)
ax1 =
    fig[1, 1] = Axis(
        fig,
        xlabel = L"\tau",
        ylabel = L"1 - (\sigma / \sigma_0)^6",
        yscale = log10,
        xscale = log10,
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorticks = IntervalsBetween(10),
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
    )

τs = Full.τs
σs = [(x.σs |> vec) for x in [T50, T2, T101, T1, T051, T05, T005, Full]]
idx = [amp_idx(x) for x in σs]
colors = [my_red, my_vermillion, my_orange, my_yellow, my_green, my_sky, my_blue, my_black]
labs = [L"\tau_0 = 50", L"\tau_0 = 2", L"\tau_0 = 1.01", L"\tau_0 = 1", L"\tau_0 = 0.51", L"\tau_0 = 0.50", L"\tau_0 = 0.05", L"\tau_0 = \infty"]

for ii = 1:length(idx)
    p = σs[ii]
    ind = idx[ii]
    if ii == length(idx)
        mk = :+
    else
        mk = :circle
    end
    scatter!(
        ax1,
        τs[ind],
        abs.(p[1] .^ 6 .- abs.(p[ind]) .^ 6) ./ p[1] .^ 6,
        markersize = 15,
        color = colors[ii],
        label = labs[ii],
        marker = mk,
    )
end

axislegend(ax1, position = :rb, labelsize = 32)
ax1.xticks = [1, 10, 100, 1000]
fig
save("Memory_Amp.pdf", fig)
