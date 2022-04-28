using Distributed
proc_num = 5
addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

files = [
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0500_μ2.0_d60_ωT2500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT100.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT250.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT500.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT1000.0_τ1000.jld2",
    "data/Multi_Thermal/Multi_25_τ0Inf_λ4_Φ0-500_μ2.0_d60_ωT2500.0_τ1000.jld2",
]

for f in files
    data = load_object(f)

    σs = data.σs
    ρs = data.ρs
    τs = data.τs
    δ = τs[2] - τs[1]

    max_τ_index = length(τs) / 10 |> floor |> Int
    
    nPts = 5000
    τs_idx = 0:Int(floor(max_τ_index / nPts)):max_τ_index |> collect
    τs_corr = τs_idx .* δ
    
    if (
        !isfile(
            "data/Correlations/Corr_Φ0$(data.Φ)_λ$(data.λ)_μ$(data.μ)_ωT$(data.ωT).jld2",
        )
    )
        res = @showprogress pmap(x -> auto_corr(σs, x), τs_idx)
        save_object(
            "data/Correlations/Corr_Φ0$(data.Φ)_λ$(data.λ)_μ$(data.μ)_ωT$(data.ωT).jld2",
            (res, τs_corr),
        )
    end
end
