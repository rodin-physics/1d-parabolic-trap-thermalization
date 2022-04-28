using Distributed
using Random

proc_num = 9
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")
system = load_object("precomputed/systems/System_ωmin2_ωmax20_d60.jld2")
## Width and Depth Dependence
λ = [1, 2, 4, 8]                                       # Interaction width
Φ = [-2000, -1000, -500, -250, 250, 500, 1000, 2000]  # Interaction height
parameters = [(x, y) for x in λ, y in Φ] |> vec

σ0 = [100]      # Starting point
μ = 2.0         # Mass of the chain particles

τ0 = Inf        # Memory
τ = 250         # Simulation time
d = 60          # Number of time steps for the fastest chain mode

@showprogress pmap(parameters) do param
    λ = param[1]
    Φ0 = param[2]

    ωmin = system.ωmin |> Float64
    ωmax = system.ωmax |> Float64
    δ = system.δ |> Float64

    n_pts = floor(τ / δ) |> Int  # Number of time steps given τ and δ
    ρHs = zeros(n_pts)

    tTraj = ThermalTrajectory(ωmin, ωmax, μ, δ, ρHs, nothing)
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, σ0, tTraj, τ0, τ)
        save_object(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end


## Mass Dependence
λ = 2                                       # Interaction width
Φ = [-500, 500]                             # Interaction height

σ0 = [100]                      # Starting point
μ = [1 / 2, 1, 2, 5]            # Mass of the chain particles
parameters = [(x, y) for x in μ, y in Φ] |> vec

τ0 = Inf        # Memory
τ = 250         # Simulation time
d = 60          # Number of time steps for the fastest chain mode

@showprogress pmap(parameters) do param
    μ = param[1]
    Φ0 = param[2]

    ωmin = system.ωmin |> Float64
    ωmax = system.ωmax |> Float64
    δ = system.δ |> Float64

    n_pts = floor(τ / δ) |> Int  # Number of time steps given τ and δ
    ρHs = zeros(n_pts)

    tTraj = ThermalTrajectory(ωmin, ωmax, μ, δ, ρHs, nothing)
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, σ0, tTraj, τ0, τ)
        save_object(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Memory dependence
λ = 4
Φ0 = -500
μ = 2.0
τ = 250

ωmin = system.ωmin |> Float64
ωmax = system.ωmax |> Float64
δ = system.δ |> Float64

n_pts = floor(τ / δ) |> Int  # Number of time steps given τ and δ
ρHs = zeros(n_pts)

tTraj = ThermalTrajectory(ωmin, ωmax, μ, δ, ρHs, nothing)

σ0 = [73.4561560502402]
if (
    !isfile(
        "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
    )
)
    res = motion_solver(system, Φ0, λ, σ0, tTraj, τ0, τ)
    save_object(
        "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
        res,
    )
end

σ0 = [100]
mems = [0, 0.05, 0.1, 0.25, 0.5, 0.51, 0.55, 1, 1.01, 1.05, 1.55, 2, 10, 50]

@showprogress pmap(mems) do τ0
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, Φ0, λ, σ0, tTraj, τ0, τ)
        save_object(
            "data/Single_Non_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end


## Multiparticle
Φs = [-500, 500]
λ = 4
μ = 2.0

# Starting positions of the mobile particles
Random.seed!(50)
σ0 = 20 .* randn(25) .+ 100
τ0 = Inf        # Memory
d = 60          # Number of time steps for the fastest chain mode
τ = 1000

@showprogress pmap(Φs) do Φ0
    ωmin = system.ωmin |> Float64
    ωmax = system.ωmax |> Float64
    δ = system.δ |> Float64

    n_pts = floor(τ / δ) |> Int  # Number of time steps given τ and δ
    ρHs = zeros(n_pts)

    tTraj = ThermalTrajectory(ωmin, ωmax, μ, δ, ρHs, nothing)
    if (
        !isfile(
            "data/Multi_Non_Thermal/Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, σ0, tTraj, τ0, τ)
        save_object(
            "data/Multi_Non_Thermal/Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end