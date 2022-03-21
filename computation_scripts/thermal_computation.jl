using Distributed
using Random
proc_num = 8
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

## Single particle

system = load_object("precomputed/systems/System_K1_k20_m1.0_d60.jld2")
d = 60
tTraj = [
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT1.0e-5_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT1.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT2.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT5.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT10.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT50.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT100.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT250.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT500.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT1000.0_τ1000_hbar0.005.jld2"),
]

mems = [Inf, 50, 1, 0.05, 0]
Fs = [-1, 1]
s = 1 / 4
M = 1 / 2           # Mass of the mobile particles
K_M = 1 / 10        # Trap force constant

# Starting positions of the mobile particles
x0 = [10]
τ = 1000

parameters = reshape([(t, F, mem) for t in tTraj, F in Fs, mem in mems], (:, 1)) |> vec

## Calculation
@showprogress pmap(parameters) do param
    traj = param[1]
    F = param[2]
    mem = param[3]

    if (
        !isfile(
            "data/Single_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(traj.m)_d$(d)_ΩT$(traj.ΩT)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, traj, mem, τ)

        save_object(
            "data/Single_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(res.m)_d$(d)_ΩT$(res.ΩT)_τ$(τ).jld2",
            res,
        )
    end

end

## Multiple particles

system = load_object("precomputed/systems/System_K1_k20_m1.0_d60.jld2")
d = 60
tTraj = [
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT5.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT10.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT50.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT100.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT250.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT500.0_τ1000_hbar0.005.jld2"),
    load_object("precomputed/rH/rH_K1_k20_m1_d60_ΩT1000.0_τ1000_hbar0.005.jld2"),
]

mems = [0, 0.05, Inf]
Fs = [-1, 1]
s = 1 / 4
M = 1 / 2           # Mass of the mobile particles
K_M = 1 / 10        # Trap force constant

# Starting positions of the mobile particles
Random.seed!(50)
x0 = 2 .* randn(25) .+ 10
τ = 1000

parameters = reshape([(t, F, mem) for t in tTraj, F in Fs, mem in mems], (:, 1)) |> vec

## Calculation
@showprogress pmap(parameters) do param
    traj = param[1]
    F = param[2]
    mem = param[3]
    if (
        !isfile(
            "data/Multi_Thermal/Multi_$(length(x0))_Mem$(mem)TM_s$(s)_F$(F)_m$(traj.m)_d$(d)_ΩT$(traj.ΩT)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, traj, mem, τ)

        save_object(
            "data/Multi_Thermal/Multi_$(length(x0))_Mem$(mem)TM_s$(s)_F$(F)_m$(res.m)_d$(d)_ΩT$(res.ΩT)_τ$(τ).jld2",
            res,
        )
    end
end

# Starting positions of the mobile particles
Random.seed!(10)
x0 = 0.1 .* randn(25)
τ = 1000

mems = [Inf]
Fs = [-1, 1]
s = 1 / 4
M = 1 / 2           # Mass of the mobile particles
K_M = 1 / 10        # Trap force constant

parameters = reshape([(t, F, mem) for t in tTraj, F in Fs, mem in mems], (:, 1)) |> vec

## Calculation
map(parameters) do param
    # @showprogress pmap(parameters) do param
    traj = param[1]
    F = param[2]
    mem = param[3]
    if (
        !isfile(
            "data/Multi_Thermal/WarmingUp_Multi_$(length(x0))_Mem$(mem)TM_s$(s)_F$(F)_m$(traj.m)_d$(d)_ΩT$(traj.ΩT)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, traj, mem, τ)

        save_object(
            "data/Multi_Thermal/WarmingUp_Multi_$(length(x0))_Mem$(mem)TM_s$(s)_F$(F)_m$(res.m)_d$(d)_ΩT$(res.ΩT)_τ$(τ).jld2",
            res,
        )
    end
end