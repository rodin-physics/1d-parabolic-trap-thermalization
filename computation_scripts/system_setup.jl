using Distributed
using Random
proc_num = 8
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term

d = 60                              # Number of time steps in the fastest chain mode
M = 1 / 2                           # Mass of the mobile particles
K_M = 1 / 10                        # Trap force constant
ΩM = √(K_M / M)                     # Trap frequency
t_M = 2 * π / ΩM                    # Period of the trapped mass

k = 20                              # Spring force constant
K = 1                               # Confining potential force constant

ms = [1, 1 / 2, 1 / 4]              # Mass of the chain atoms
t_max = 1000 * t_M                  # Simulation time

for m in ms
    res = mkChainSystem(K, k, m, t_max, d)
    save_object("precomputed/systems/System_K$(K)_k$(k)_m$(m)_d$(d).jld2", res)
end

## Precompute the thermal trajectories
m = 1                           # Mass of the chain atoms
τ = 1000;
t_max = τ * t_M                 # Simulation time

Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency
δ = (2 * π / Ωmax) / d          # Time step
n_pts = floor(t_max / δ) |> Int # Number of time steps given t_max and δ

## Thermal Trajectory
n_masses = 1000000              # Number of chain masses for simulating r0
qs = range(0, π / 2, length = round(n_masses / 2) |> Integer)
Ωs = Ω.(K, k, m, qs)
ΩTs = [1e-5, 1, 2, 5, 10, 20, 50, 100, 150, 250, 500, 1000]
for ii = 1:length(ΩTs)
    println(ii)
    ΩT = ΩTs[ii]
    # Seeding the RNG
    Random.seed!(150)
    ϕs = 2 * π * rand(length(qs))
    ζs = ζq.(Ωs, ΩT, ħ)
    rHs = @showprogress pmap(n -> ζH(n, δ, ζs, ϕs, Ωs) / √(m), 1:n_pts)
    res = ThermalTrajectory(k, K, m, δ, rHs, ΩT, ħ)
    save_object(
        "precomputed/rH/rH_K$(K)_k$(k)_m$(m)_d$(d)_ΩT$(ΩT)_τ$(τ)_hbar$(ħ).jld2",
        res,
    )
end