using CairoMakie
using Colors
using Distributions
using JLD2
using LaTeXStrings
using LinearAlgebra
using ProgressMeter
using QuadGK
using Statistics
using StatsBase
## Parameters
η = 1e-12                       # Small number for integration
## Friendly colors
my_red = colorant"rgba(204, 121, 167, 0.75)"
my_vermillion = colorant"rgba(213, 94, 0, 0.75)"
my_orange = colorant"rgba(230, 159, 0, 0.75)"
my_yellow = colorant"rgba(240, 228, 66, 0.75)"
my_green = colorant"rgba(0, 158, 115, 0.75)"
my_sky = colorant"rgba(86, 180, 233, 0.75)"
my_blue = colorant"rgba(0, 114, 178, 0.75)"
my_black = colorant"rgba(0, 0, 0, 0.35)"

## Types
struct ChainSystem
    ωmin::Float64               # Smallest mode frequency
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Time step
    Γ::Vector{Float64}          # Response array
end

struct ThermalTrajectory
    ωmin::Float64               # Smallest mode frequency
    ωmax::Float64               # Largest mode frequency
    μ::Float64                  # Mass of the chain atoms
    δ::Float64                  # Time step
    ρHs::Vector{Float64}        # Thermal trajectory
    ωT::Union{Nothing,Float64}  # Temperature
end

struct SystemSolution
    ωmin::Float64               # Smallest mode frequency
    ωmax::Float64               # Largest mode frequency
    μ::Float64                  # Mass of the chain atoms
    τs::Vector{Float64}         # Time steps
    τ0::Float64                 # Memory
    Φ::Float64                  # Magnitude of the Gaussian potential
    λ::Float64                  # Standard deviation of the potential
    σs::Matrix{Float64}         # Positions of mobile atoms. 
    ρs::Vector{Float64}         # Positions of the chain atom
    ωT::Union{Nothing,Float64}  # Temperature
end

## Functions
# Frequency as a function of momentum
@inline function Ω(ωmin, ωmax, x)
    return sqrt(ωmin^2 * cos(x)^2 + ωmax^2 * sin(x)^2)
end

# Chain recoil function
function Γ(τ, ωmin, ωmax)
    int_fun(x) = sin(2 * π * τ * x) / (√((x^2 - ωmin^2) * (ωmax^2 - x^2)))
    res = quadgk(int_fun, ωmin + η, ωmax - η)
    return (res[1] * 2 / π)
end

# Function for assembiling a ChainSystem
function mkChainSystem(ωmin, ωmax, τ_max, d)
    δ = (1 / ωmax) / d                  # Time step in units of t_M
    n_pts = floor(τ_max / δ) |> Int     # Number of time steps given τ_max and δ
    # Precomputing the memory term
    Γ_list = @showprogress pmap(jj -> Γ(δ * jj, ωmin, ωmax), 1:n_pts)
    return ChainSystem(ωmin, ωmax, δ, Γ_list)
end

# Mode amplitude
function ζq(ωq, ωT, μ)
    # Subtract a small number from p. The reason is that for low ωT, p ≈ 1,
    # causing issues with the rand() generator
    n = rand(Geometric(1 - exp(-ωq / ωT) - η))
    res = √(n + 1 / 2) * √(2 / ωq / μ)
    return res
end

# Homogeneous displacement of the active chain atom at time step n given a set of ωs
# and the corresponding ζs and ϕs as a sum of normal coordinates. 
function ζH(n, δ, ζs, ϕs, ωs)
    n_ζ = length(ζs)
    res = [ζs[x] * cos(2 * π * δ * n * ωs[x] + ϕs[x]) / √(n_ζ) for x = 1:n_ζ] |> sum
    return res
end

# Function that calculates the trajectories of the mobile atoms and the chain particle.
# τ0 determines the memory length and τ is the simulation length, both in units of
# the trap period.
function motion_solver(
    system::ChainSystem,
    Φ0::T where {T<:Real},
    λ::T where {T<:Real},
    σ0::Vector{T} where {T<:Real},
    tTraj::ThermalTrajectory,
    τ0::T where {T<:Real},
    τ::T where {T<:Real},
)
    ωmin = system.ωmin      # Spring force constant
    ωmax = system.ωmax      # Confining potential force constant
    δ = system.δ            # Time step
    Γ_list = system.Γ       # Memory term
    μ = tTraj.μ             # Mass of the mobile particle
    # Check that the thermal trajectory is for the correct system
    if (ωmin != tTraj.ωmin || ωmax != tTraj.ωmax || δ != tTraj.δ)
        error("Thermal trajectory describes a different system. Check your input.")
    else
        ρHs = tTraj.ρHs
    end

    n_pts = floor(τ / δ) |> Int         # Number of time steps
    τs = δ .* (1:n_pts) |> collect      # Times
    if length(ρHs) < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
    end

    τ0_pts = max(floor(τ0 / δ), 1)    # Memory time points.
    # Even if τ0 == 0, we have to keep a single time point to make sure arrays work.
    # For zero memory, the recoil contribution is dropped (see line 170)

    # If the precomputed memory is shorter than the simulation time AND shorter
    # than the desired memory, terminate the calculation.
    if (length(Γ_list) < n_pts && length(Γ_list) < τ0_pts)
        error("Chosen memory and the simulation time exceed the precomputed range.")
    else
        # The number of memory pts can be limited by the total simulation time.
        τ0_pts = min(τ0_pts, n_pts) |> Int
        Γ_list = (2 * π * δ / μ) .* Γ_list[1:τ0_pts] |> reverse
    end

    # Interaction terms
    function dU_dρ(r)
        return (-Φ0 * exp(-r^2 / (2 * λ^2)) * r / λ^2)
    end

    ## Arrays
    σs = zeros(n_pts, length(σ0))         # Mobile mass position
    ρs = ρHs[1:n_pts]                     # Chain mass position
    U_pr_mob = zeros(n_pts, length(σ0))   # Force on the mobile particles
    U_pr_chain = zeros(n_pts)             # Force on the chain atom

    ## Initial values
    σs[1, :] = σ0              # Initialize the starting R
    σs[2, :] = σ0              # Starting at rest
    # Initialize the interaction terms
    U_pr1 = dU_dρ.(ρs[1] .- σs[1, :])
    U_pr2 = dU_dρ.(ρs[2] .- σs[2, :])
    U_pr_mob[1, :] = -U_pr1
    U_pr_mob[2, :] = -U_pr2
    U_pr_chain[1] = sum(U_pr1)
    U_pr_chain[2] = sum(U_pr2)

    for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index
        # Get the memory term elements. If the number of elapsed steps is smaller
        # than the number of elements in the memory term, take the corresponding
        # number of the most recent elements
        Γs = Γ_list[(τ0_pts-min(τ0_pts, curr)+1):end]
        # Get the interaction force elements. If the number of elapsed steps is
        # smaller than the number of elements in the memory term, take the available
        # values.
        Us = U_pr_chain[(curr-min(τ0_pts, curr)+1):curr]
        # If mem == 0, drop the recoil contribution
        ρs[nxt] -= dot(Γs, Us) * (τ0 != 0)
        σs[nxt, :] =
            (2 * π * δ)^2 .* (-U_pr_mob[curr, :] - σs[curr, :]) + 2 .* σs[curr, :] -
            σs[curr-1, :]
        U_pr = dU_dρ.(ρs[nxt] .- σs[nxt, :])
        U_pr_chain[nxt] = sum(U_pr)
        U_pr_mob[nxt, :] = -U_pr
    end
    return SystemSolution(ωmin, ωmax, μ, τs, τ0, Φ0, λ, σs, ρs, tTraj.ωT)
end

function auto_corr(l, lag)
    res = cor(l[1:(end-lag), :], l[(1+lag):end, :])
    return res
end
