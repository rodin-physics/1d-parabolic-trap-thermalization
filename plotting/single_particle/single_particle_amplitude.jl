include("../../src/main.jl")

function amp_idx(σ)
    return filter(x -> (σ[x] > σ[x-1]) && (σ[x] > σ[x+1]), 2:(length(σ)-1))
end
C_att = 5e-9
C_rep = 7e-9

fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)

ax2 = Axis(
    fig[2, 1],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)

ax3 = Axis(
    fig[1, 2],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)

ax4 = Axis(
    fig[2, 2],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)

ax5 = Axis(
    fig[1, 3],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)

ax6 = Axis(
    fig[2, 3],
    xlabel = L"\tau",
    ylabel = L"1 - (\sigma/\sigma_0)^6",
    yscale = log10,
    xscale = log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
)
ts_fit = exp.(range(-1, 6, length = 20))

colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)

# Mass dependence
data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ0.5_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ1.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ5.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\mu = 1/2", L"\mu = 1", L"\mu = 2", L"\mu = 5"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax1,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax1,
        τs[idx],
        C_att * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )

end

data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ0.5_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ1.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ5.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\mu = 1/2", L"\mu = 1", L"\mu = 2", L"\mu = 5"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax2,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax2,
        τs[idx],
        C_rep * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )

end
# Width Scaling

data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ1_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ8_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\lambda = 1", L"\lambda = 2", L"\lambda = 4", L"\lambda = 8"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax3,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax3,
        τs[idx],
        C_att * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )


end

data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ1_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ4_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ8_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\lambda = 1", L"\lambda = 2", L"\lambda = 4", L"\lambda = 8"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax4,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax4,
        τs[idx],
        C_rep * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )

end
## Depth Scaling
data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-250_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-1000_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-2000_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\Phi_0 = -250", L"\Phi_0 = -500", L"\Phi_0 = -1000", L"\Phi_0 = -2000"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax5,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax5,
        τs[idx],
        C_att * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )


end

data = [
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0250_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ01000_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
    load_object(
        "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ02000_μ2.0_d60_ωTnothing_τ250.jld2",
    ),
]
labs = [L"\Phi_0 = 250", L"\Phi_0 = 500", L"\Phi_0 = 1000", L"\Phi_0 = 2000"]
colors = [my_vermillion, my_orange, my_green, my_blue]
for n = 1:length(data)
    d = data[n]
    τs = d.τs
    σs = abs.(d.σs |> vec)
    idx = amp_idx(σs)
    scatter!(
        ax6,
        τs[idx],
        abs.(σs[idx] .^ 6 .- σs[1] .^ 6) ./ σs[1] .^ 6,
        color = colors[n],
        markersize = 5,
        label = labs[n],
    )
    lines!(
        ax6,
        τs[idx],
        C_rep * (d.Φ * d.λ)^2 / (d.μ) * τs[idx],
        color = colors[n],
        linewidth = 2,
    )

end

txt = ["(a)", "(d)", "(b)", "(e)", "(c)", "(f)"]
axs = [ax1, ax2, ax3, ax4, ax5, ax6]
for n = 1:length(axs)
    ax = axs[n]
    ax.xticks = [1, 10, 100, 1000]
    xlims!(ax, (4e-1, 250))
    ylims!(ax, (1e-4, 1.25))
    text!(ax, txt[n], position = (1, 1 / 2), textsize = 16, font = "CMU Serif")
    axislegend(ax, position = :rb)

end

save("Amplitude_Decay.pdf", fig)


# Attractive vs repulsive
repulsive = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0500_μ2.0_d60_ωTnothing_τ250.jld2",
)
attractive = load_object(
    "data/Single_Non_Thermal/Single_σ0[100]_τ0Inf_λ2_Φ0-500_μ2.0_d60_ωTnothing_τ250.jld2",
)

σs_attractive = attractive.σs |> vec
σs_repulsive = repulsive.σs |> vec

ρs_attractive = attractive.ρs |> vec
ρs_repulsive = repulsive.ρs |> vec
Φ0_attractive = attractive.Φ
Φ0_repulsive = repulsive.Φ
λ = attractive.λ

force_attractive =
    -Φ0_attractive .* exp.(-(ρs_attractive - σs_attractive) .^ 2 ./ (2 * λ^2)) .*
    (ρs_attractive - σs_attractive) ./ λ^2
force_repulsive =
    -Φ0_repulsive .* exp.(-(ρs_repulsive - σs_repulsive) .^ 2 ./ (2 * λ^2)) .*
    (ρs_repulsive - σs_repulsive) ./ λ^2

n_pts = length(σs_attractive)
τs = attractive.τs
idx = findall(x -> (x >= 0) && (x <= 1), τs)

fig = Figure(resolution = (400, 300), font = "CMU Serif", fontsize = 14)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"-dU/d\sigma")
ax2 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\rho", yaxisposition = :right)


lines!(ax1, τs[idx], force_attractive[idx], color = my_vermillion, linewidth = 2, label = L"\Phi_0 < 0")
lines!(ax1, τs[idx], force_repulsive[idx], color = my_blue, linewidth = 2, label = L"\Phi_0 > 0")
lines!(ax2, τs[idx], ρs_attractive[idx], color = my_vermillion, linewidth = 2, linestyle = "-")
lines!(ax2, τs[idx], ρs_repulsive[idx], color = my_blue, linewidth = 2, linestyle = "-")
xlims!(ax1, (0.23, 0.27))
xlims!(ax2, (0.23, 0.27))
ylims!(ax1, (-200, 200))
ylims!(ax2, (-0.2, 0.2))
l = axislegend(ax1, position = :lb)
save("Attractive_Repulsive.pdf", fig)