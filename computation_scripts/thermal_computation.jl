using Distributed
using Random
proc_num = 7
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

## Single particle

system = load_object("precomputed/systems/System_ωmin2_ωmax20_d60.jld2")
d = 60
tTraj = [
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT1.0e-5_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT2.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT5.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT10.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT50.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT100.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT250.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT500.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT1000.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT2500.0_τ1000.jld2"),
]
τ0s = [Inf, 50, 1, 0.05, 0]
Φs = [-500, 500]
λ = 4
μ = 2.0

# Starting positions of the mobile particles
σ0 = [100]
τ = 1000
parameters = reshape([(t, Φ, τ0) for t in tTraj, Φ in Φs, τ0 in τ0s], (:, 1)) |> vec

## Calculation
@showprogress pmap(parameters) do param
    traj = param[1]
    Φ0 = param[2]
    τ0 = param[3]
    if (
        !isfile(
            "data/Single_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
        )
    )

        @time res = motion_solver(system, Φ0, λ, σ0, traj, τ0, τ)
        save_object(
            "data/Single_Thermal/Single_σ0$(σ0)_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
            res,
        )
    end
end

## Multiple particles


system = load_object("precomputed/systems/System_ωmin2_ωmax20_d60.jld2")
d = 60
tTraj = [
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT1.0e-5_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT1.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT2.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT5.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT10.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT50.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT100.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT250.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT500.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT1000.0_τ1000.jld2"),
    load_object("precomputed/rH/rH_ωmin2_ωmax20_μ2.0_d60_ωT2500.0_τ1000.jld2"),
]

τ0s = [Inf, 0.05, 0]
Φs = [-500, 500]
λ = 4
μ = 2.0

# Starting positions of the mobile particles
Random.seed!(50)
σ0 = 20 .* randn(25) .+ 100
τ = 1000

parameters = reshape([(t, Φ, τ0) for t in tTraj, Φ in Φs, τ0 in τ0s], (:, 1)) |> vec

## Calculation
@showprogress pmap(parameters) do param
    traj = param[1]
    Φ0 = param[2]
    τ0 = param[3]
    if (
        !isfile(
            "data/Multi_Thermal/Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, Φ0, λ, σ0, traj, τ0, τ)

        save_object(
            "data/Multi_Thermal/Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
            res,
        )
    end
end

# Starting positions of the mobile particles
τ0s = [Inf]
Φs = [-500, 500]
λ = 4
μ = 2.0

# Starting positions of the mobile particles
Random.seed!(10)
σ0 = 1 .* randn(25)
τ = 1000

parameters = reshape([(t, Φ, τ0) for t in tTraj, Φ in Φs, τ0 in τ0s], (:, 1)) |> vec

## Calculation
@showprogress pmap(parameters) do param
    traj = param[1]
    Φ0 = param[2]
    τ0 = param[3]
    if (
        !isfile(
            "data/Multi_Thermal/WarmingUp_Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, Φ0, λ, σ0, traj, τ0, τ)

        save_object(
            "data/Multi_Thermal/WarmingUp_Multi_$(length(σ0))_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_μ$(μ)_d$(d)_ωT$(traj.ωT)_τ$(τ).jld2",
            res,
        )
    end
end