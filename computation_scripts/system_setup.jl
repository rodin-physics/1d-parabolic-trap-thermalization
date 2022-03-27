using Distributed
using Random
proc_num = 9
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term

d = 60          # Number of time steps in the fastest chain mode
ωmin = 2        # Smallest mode frequency
ωmax = 20       # Largest mode frequency
τ_max = 1000    # Simulation time
if (!isfile("precomputed/systems/System_ωmin$(ωmin)_ωmax$(ωmax)_d$(d).jld2"))
    res = mkChainSystem(ωmin, ωmax, τ_max, d)
    save_object("precomputed/systems/System_ωmin$(ωmin)_ωmax$(ωmax)_d$(d).jld2", res)
end

## Precompute the thermal trajectories
μ = 2.0                         # Mass of the chain atoms
τ = 1000                        # Time length
δ = (1 / ωmax) / d              # Time step in units of t_M
n_pts = floor(τ / δ) |> Int      # Number of time steps given t_max and δ

## Thermal Trajectory
n_masses = 1000000              # Number of chain masses for simulating ρ0
qs = range(0, π / 2, length = round(n_masses / 2) |> Integer)
ωs = Ω.(ωmin, ωmax, qs)
ωTs = [1e-5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 2500]
for ii = 1:length(ωTs)
    println(ii)
    ωT = ωTs[ii]

    if (!isfile("precomputed/rH/rH_ωmin$(ωmin)_ωmax$(ωmax)_μ$(μ)_d$(d)_ωT$(ωT)_τ$(τ).jld2"))
        # Seeding the RNG
        Random.seed!(150)
        ϕs = 2 * π * rand(length(ωs))
        ζs = ζq.(ωs, ωT, μ)
        ρHs = @showprogress pmap(n -> ζH(n, δ, ζs, ϕs, ωs), 1:n_pts)
        res = ThermalTrajectory(ωmin, ωmax, μ, δ, ρHs, ωT)
        save_object(
            "precomputed/rH/rH_ωmin$(ωmin)_ωmax$(ωmax)_μ$(μ)_d$(d)_ωT$(ωT)_τ$(τ).jld2",
            res,
        )
    end
end
