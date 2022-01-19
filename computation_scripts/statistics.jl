using Distributed
proc_num = 5
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

files = [
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F-1_m1.0_d60_ΩT10.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT50.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_MemInfTM_s0.25_F1_m1.0_d60_ΩT10.0_τ1000.jld2",
]

for f in files
    data = load_object(f)

    ΩT = data.ΩT
    F = data.F
    s = data.s
    Rs = data.Rs
    rs = data.rs
    k = data.k
    K = data.K
    ts = data.ts
    δ = ts[2] - ts[1]
    K_M = data.K_M
    M = data.M
    m = data.m
    ħ = data.ħ

    Ωmin = sqrt(K / m)                  # Smallest chain frequency
    Ωmax = sqrt(4 * k / m + K / m)      # Largest chain frequency

    ΩM = √(K_M / M)                 # Trap frequency
    t_M = 2 * π / ΩM                # Period of the trapped mass

    times = data.ts ./ t_M
    Rs = data.Rs
    rs = data.rs
    nPts = 5000
    num_τs = length(times) / 10 |> floor |> Int
    τs = 0:Int(floor(num_τs / nPts)):num_τs |> collect
    ts_corr = τs * (times[2] - times[1])

    prt = size(Rs)[2]

    kin_en = vcat(
        zeros(1, prt),
        ((Rs[2:end, :] - Rs[1:(end-1), :]) ./ δ) .^ 2 ./ 2 .* M ./ (ħ * ΩT),
    )
    pot_en = (Rs .^ 2 / 2 * K_M + F .* exp.(-(rs .- Rs) .^ 2 ./ (2 * s^2))) ./ (ħ * ΩT)
    tot_en = (pot_en + kin_en)

    if (
        !isfile(
            "data/Correlations/Energy_Corr_F$(data.F)_K$(data.K)_k$(data.k)_m$(data.m)_ΩT$(data.ΩT)_hbar$(data.ħ).jld2",
        )
    )
        res = @showprogress pmap(x -> auto_corr(tot_en, x), τs)
        save_object(
            "data/Correlations/Energy_Corr_F$(data.F)_K$(data.K)_k$(data.k)_m$(data.m)_ΩT$(data.ΩT)_hbar$(data.ħ).jld2",
            (res, ts_corr),
        )
    end

    if (
        !isfile(
            "data/Correlations/Corr_F$(data.F)_K$(data.K)_k$(data.k)_m$(data.m)_ΩT$(data.ΩT)_hbar$(data.ħ).jld2",
        )
    )
        res = @showprogress pmap(x -> auto_corr(Rs, x), τs)
        save_object(
            "data/Correlations/Corr_F$(data.F)_K$(data.K)_k$(data.k)_m$(data.m)_ΩT$(data.ΩT)_hbar$(data.ħ).jld2",
            (res, ts_corr),
        )
    end
end
