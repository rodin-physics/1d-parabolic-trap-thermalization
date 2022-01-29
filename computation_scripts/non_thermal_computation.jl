using Distributed

proc_num = 9
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

d = 60

M = 1 / 2                           # Mass of the mobile particles
K_M = 1 / 10                        # Trap force constant
ΩM = √(K_M / M)                     # Trap frequency
t_M = 2 * π / ΩM                    # Period of the trapped mass

## Width and Depth Dependence
system = load_object("precomputed/systems/System_K1_k20_m1.0_d60.jld2")

s = [1 / 2, 1 / 4, 1 / 8, 1 / 16]
F = [-1, 1]
parameters = [(x, y) for x in s, y in F] |> vec

x0 = [10]

mem = Inf
τ = 250
t_max = τ * t_M                         # Simulation time

@showprogress pmap(parameters) do param
    s = param[1]
    F = param[2]

    m = system.m |> Float64
    k = system.k |> Float64
    K = system.K |> Float64
    δ = system.δ |> Float64

    n_pts = floor(t_max / δ) |> Int  # Number of time steps given t_max and δ
    rHs = zeros(n_pts)

    tTraj = ThermalTrajectory(k, K, m, δ, rHs, nothing, ħ)
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, tTraj, mem, τ)
        save_object(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

s = [1 / 4]
F = [-4, -2, -1, -1/2, 1/2, 1,2,4]
parameters = [(x, y) for x in s, y in F] |> vec

x0 = [10]

mem = Inf
τ = 250
t_max = τ * t_M                         # Simulation time

@showprogress pmap(parameters) do param
    s = param[1]
    F = param[2]

    m = system.m |> Float64
    k = system.k |> Float64
    K = system.K |> Float64
    δ = system.δ |> Float64

    n_pts = floor(t_max / δ) |> Int  # Number of time steps given t_max and δ
    rHs = zeros(n_pts)

    tTraj = ThermalTrajectory(k, K, m, δ, rHs, nothing, ħ)
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, tTraj, mem, τ)
        save_object(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Mass Dependence

systems = [
    load_object("precomputed/systems/System_K1_k20_m0.25_d60.jld2"),
    load_object("precomputed/systems/System_K1_k20_m0.5_d60.jld2"),
    load_object("precomputed/systems/System_K1_k20_m1.0_d60.jld2"),
]

s = 1 / 16
F = -1

x0 = [10]

mem = Inf
τ = 250
t_max = τ * t_M                         # Simulation time

@showprogress pmap(systems) do system
    m = system.m |> Float64
    k = system.k |> Float64
    K = system.K |> Float64
    δ = system.δ |> Float64

    n_pts = floor(t_max / δ) |> Int  # Number of time steps given t_max and δ
    rHs = zeros(n_pts)

    tTraj = ThermalTrajectory(k, K, m, δ, rHs, nothing, ħ)

    if (
        !isfile(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, tTraj, mem, τ)
        save_object(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Memory dependence

system = load_object("precomputed/systems/System_K1_k20_m1.0_d60.jld2")
d = 60

s = 1 / 4
F = -1

τ = 250
t_max = τ * t_M                         # Simulation time

m = system.m |> Float64
k = system.k |> Float64
K = system.K |> Float64
δ = system.δ |> Float64

n_pts = floor(t_max / δ) |> Int  # Number of time steps given t_max and δ
rHs = zeros(n_pts)

tTraj = ThermalTrajectory(k, K, m, δ, rHs, nothing, ħ)

x0 = [7.422582732835284]
res = motion_solver(system, F, s, M, K_M, x0, tTraj, mem, τ)
save_object(
    "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
    res,
)

x0 = [10]
mems = [0, 0.05, 0.1, 0.25, 0.5, 0.51, 0.55, 1, 1.01, 1.05, 1.55, 2, 10, 50]

@showprogress pmap(mems) do mem
    if (
        !isfile(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, K_M, x0, tTraj, mem, τ)
        save_object(
            "data/Single_Non_Thermal/Single_R0$(x0)_Mem$(mem)TM_s$(s)_F$(F)_m$(m)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end
