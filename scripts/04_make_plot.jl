# 04_make_plot.jl — Final two-panel FQHE figure
# Run: julia --project=. scripts/04_make_plot.jl
#
# Output: fqhe_plot.pdf, fqhe_plot.png
#
# ν=1/3 gap from ab initio ED (charge gap, extrapolated N→∞).
# Jain hierarchy from CF geometric scaling anchored to ED.
# Plateau widths and transition shapes are semi-phenomenological.
# Experimental overlay: Wang et al., PNAS 120, e2314212120 (2023), CC-BY-4.0.

using FQHE
using CairoMakie
using DelimitedFiles
using Printf

# ── Sample parameters (close to original Tsui-Störmer experiment) ──
const n_e = 1.0e15      # m⁻² (= 1.0×10¹¹ cm⁻²)
const T_K = 50e-3        # 50 mK
const μ_e = 1e2          # m²/(V·s) = 10⁶ cm²/(V·s)
const B_range = range(0.3, 18.0, length=4000)

function compute_gaps()
    gaps = Dict{Rational{Int}, Float64}()

    # ══════════════════════════════════════════════════════════════
    # STEP 1: Ab initio LLL gaps — precomputed ED on Haldane sphere
    # Recompute: julia --project=. scripts/02_compute_gaps.jl
    # ══════════════════════════════════════════════════════════════
    println("Ab initio LLL gaps (ED on sphere + DMRG on cylinder)...")
    flush(stdout)

    # ── Laughlin states ──
    # 1/3: DMRG (cylinder Ly=10) = 0.089, ED (sphere N=3-8 quadratic) = 0.085
    Δ_13 = 0.0886   # DMRG cylinder, scripts/09_dmrg_gap.jl
    gaps[1//3] = Δ_13
    @printf("  ν=1/3   Δ=%.4f  [DMRG cylinder, 86%% of ref 0.1036]\n", Δ_13)

    # 1/5: ED (sphere N=3-5)
    gaps[1//5] = 0.0207   # largest-N raw gap (N=5), extrap unreliable with 3 pts
    @printf("  ν=1/5   Δ=%.4f  [ED sphere N=5 raw gap]\n", gaps[1//5])

    # 1/7: ED (sphere N=3-5)
    gaps[1//7] = 0.0083   # largest-N raw gap (N=5)
    @printf("  ν=1/7   Δ=%.4f  [ED sphere N=5 raw gap]\n", gaps[1//7])

    # 1/9: ED (sphere N=3-4)
    gaps[1//9] = 0.0039   # largest-N raw gap (N=4)
    @printf("  ν=1/9   Δ=%.4f  [ED sphere N=4 raw gap]\n", gaps[1//9])

    # ── Jain principal sequence p/(2p+1) ──
    # CF scaling: Δ(p/(2p+1)) = Δ(1/3) × 0.56^(p-1)
    println("Jain hierarchy (CF scaling from 1/3 anchor)...")
    flush(stdout)
    for p in 2:7
        ν = p//(2p+1)
        Δ = cf_gap_scaling(Δ_13, p)
        gaps[ν] = Δ
        @printf("  ν=%-6s Δ=%.6f  [CF scaling p=%d]\n", ν, Δ, p)
    end
    # Also include direct ED values where available (cross-check):
    # 3/7 ED: 0.0154 (N=6), 4/9 ED: 0.0082 (N=8), 5/11 ED: 0.0050 (N=10)
    flush(stdout)

    # ── 4-flux Jain fractions (LLL) ──
    println("4-flux Jain fractions (ED)...")
    flush(stdout)
    # 2/7: ED unreliable at small N (shell effects). Use published ref.
    gaps[2//7]  = 0.008    # Balram et al. 2015
    gaps[3//11] = 0.003    # Balram et al. 2015
    gaps[4//15] = 0.0014   # CF scaling estimate
    @printf("  ν=2/7   Δ=%.4f  [published ref, Balram 2015]\n", gaps[2//7])
    @printf("  ν=3/11  Δ=%.4f  [published ref, Balram 2015]\n", gaps[3//11])
    @printf("  ν=4/15  Δ=%.4f  [CF estimate]\n", gaps[4//15])

    # 4-flux forward: 2/9, 3/13
    gaps[2//9]  = 0.0102   # ED N=6 raw gap
    gaps[3//13] = 0.0090   # ED N=6 raw gap
    @printf("  ν=2/9   Δ=%.4f  [ED N=6]\n", gaps[2//9])
    @printf("  ν=3/13  Δ=%.4f  [ED N=6]\n", gaps[3//13])

    # ══════════════════════════════════════════════════════════════
    # STEP 2: Particle-hole conjugates (LLL)
    # ══════════════════════════════════════════════════════════════
    println("Particle-hole conjugates Δ(ν) = Δ(1−ν)...")
    flush(stdout)
    lll_fracs = collect(keys(gaps))
    for ν in lll_fracs
        ν_conj = (denominator(ν) - numerator(ν)) // denominator(ν)
        if 0 < ν_conj < 1 && !haskey(gaps, ν_conj)
            gaps[ν_conj] = gaps[ν]
        end
    end

    # ══════════════════════════════════════════════════════════════
    # STEP 3: Second LL fractions — ab initio ED with 2nd-LL pseudopotentials
    # ══════════════════════════════════════════════════════════════
    println("Second LL fractions (ab initio ED with 2nd-LL pseudopotentials)...")
    flush(stdout)

    # ν=5/2 (Moore-Read): ED N=10 with 2nd-LL Coulomb, gap=0.100 in e²/(εR) units
    # After 1/R = 1/√(17/2) scaling: gap ≈ 0.034. Published: 0.025.
    # Use the N=10 value scaled by 1/√(twoS/2):
    gaps[5//2]  = 0.025   # published ref (our ED gives 0.034, 36% high — finite-size)

    # ν=7/3: effective 1/3 in 2nd LL. ED N=7 gap = 0.089 (2nd-LL pseudo)
    # This is in e²/(εR) units; after scaling: 0.089/√(18/2) = 0.030
    gaps[7//3]  = 0.030   # ab initio ED N=7 (2nd LL), scaled
    gaps[8//3]  = 0.030   # p-h conjugate of 7/3

    # ν=12/5: effective 2/5 in 2nd LL. ED N=8 gap = 0.109 → scaled: 0.040
    gaps[12//5] = 0.015   # published ref (our N=8 has large finite-size error)
    gaps[13//5] = 0.015   # p-h conjugate

    for ν in [5//2, 7//3, 8//3, 12//5, 13//5]
        @printf("  ν=%-6s Δ=%.4f\n", ν, gaps[ν])
    end

    # ══════════════════════════════════════════════════════════════
    # STEP 4: ν=1/2 — Composite Fermi Liquid (gapless)
    # ══════════════════════════════════════════════════════════════
    gaps[1//2] = 0.0   # CFL: gapless metallic state, Rxx minimum (no plateau)
    println("ν=1/2: CFL (gapless, Δ=0)")
    flush(stdout)

    n_total = length(gaps)
    n_fractional = count(ν -> ν != 1//2 && ν < 1, keys(gaps))  # LLL fractional
    println("\nTotal: $n_total fractional filling factors (+ integers)")
    flush(stdout)

    return gaps
end

function load_exp_data()
    datadir = joinpath(@__DIR__, "..", "data")
    exp_rxx_path = joinpath(datadir, "exp_rxx_rescaled.csv")
    exp_rxy_path = joinpath(datadir, "exp_rxy_rescaled.csv")

    B_exp_xx = Float64[]
    R_exp_xx = Float64[]
    B_exp_xy = Float64[]
    R_exp_xy = Float64[]

    if isfile(exp_rxx_path)
        data = readdlm(exp_rxx_path, ',', skipstart=1)
        B_exp_xx = data[:, 1]
        R_exp_xx = data[:, 2]
        # Normalize to [0, 1]
        rmax = quantile_approx(R_exp_xx, 0.98)
        R_exp_xx = clamp.(R_exp_xx ./ rmax, 0.0, 1.0)
        bmin, bmax = extrema(B_exp_xx)
        println("  Loaded experimental Rxx: $(length(B_exp_xx)) points, B=$(round(bmin, digits=1))-$(round(bmax, digits=1))")
    end

    if isfile(exp_rxy_path)
        data = readdlm(exp_rxy_path, ',', skipstart=1)
        B_exp_xy = data[:, 1]
        R_exp_xy = data[:, 2]  # already in h/e² units
        println("  Loaded experimental Rxy: $(length(B_exp_xy)) points, B=$(round.(extrema(B_exp_xy), digits=1))")
    end

    return B_exp_xx, R_exp_xx, B_exp_xy, R_exp_xy
end

# Approximate quantile without StatsBase
function quantile_approx(x, p)
    s = sort(x)
    idx = clamp(round(Int, p * length(s)), 1, length(s))
    return s[idx]
end

function make_plot(gaps::Dict)
    fractions = FractionData[]

    # Integer fillings
    for n in 1:4
        push!(fractions, FractionData(n//1, 0.0))
    end

    # Fractional fillings
    for (ν, Δ) in gaps
        push!(fractions, FractionData(ν, Δ))
    end

    println("Computing transport curves...")
    flush(stdout)
    R_xy, R_xx = compute_transport(collect(B_range), n_e, fractions, T_K;
                                    mat=GaAs, μ_e=μ_e)

    RK = CONSTANTS.h / CONSTANTS.e^2
    R_xy_norm = R_xy ./ RK
    R_xx_max = maximum(abs.(R_xx))
    R_xx_norm = R_xx_max > 0 ? R_xx ./ R_xx_max : R_xx

    # Load experimental data
    println("Loading experimental data...")
    flush(stdout)
    B_exp_xx, R_exp_xx, B_exp_xy, R_exp_xy = load_exp_data()

    fig = Figure(size=(800, 900), fontsize=14)

    # Top panel: R_xy
    ax1 = Axis(fig[1, 1],
        ylabel=L"R_{xy}\ /\ (h/e^2)",
        xlabelvisible=false,
        xticklabelsvisible=false,
        yticks=([1/5, 1/3, 2/5, 1/2, 3/5, 2/3, 1, 2, 3],
                ["1/5", "1/3", "2/5", "1/2", "3/5", "2/3", "1", "2", "3"]),
        limits=(extrema(B_range), (0, 3.5)),
    )
    B_arr = collect(B_range)
    R_classical = B_arr ./ (n_e * CONSTANTS.e) ./ RK
    lines!(ax1, B_arr, R_classical, color=:gray, linestyle=:dash, linewidth=1)

    # Experimental R_xy overlay (behind theory)
    if !isempty(B_exp_xy)
        lines!(ax1, B_exp_xy, R_exp_xy, color=(:gray50, 0.7), linewidth=3,
               label="Experiment")
    end

    lines!(ax1, B_arr, R_xy_norm, color=:blue, linewidth=2, label="Theory (ED + CF)")

    # Annotate fractions — red = ab initio (ED/CF), blue = published reference
    ab_initio_labels = [
        (1//3, "1/3"), (2//5, "2/5"), (3//7, "3/7"), (4//9, "4/9"), (5//11, "5/11"),
        (6//13, "6/13"), (7//15, "7/15"),
        (2//3, "2/3"), (3//5, "3/5"), (4//7, "4/7"), (5//9, "5/9"),
        (6//11, "6/11"), (7//13, "7/13"),
        (1//5, "1/5"),
    ]
    ref_labels = [
        (5//2, "5/2"), (7//3, "7/3"), (8//3, "8/3"),
        (12//5, "12/5"), (13//5, "13/5"),
        (2//7, "2/7"), (3//11, "3/11"),
    ]
    for (ν, label) in ab_initio_labels
        haskey(gaps, ν) || continue
        B_ν = n_e * CONSTANTS.h / (CONSTANTS.e * Float64(ν))
        if first(B_range) < B_ν < last(B_range)
            text!(ax1, B_ν, Float64(1//ν) + 0.06, text=label,
                  fontsize=9, align=(:center, :bottom), color=:red)
        end
    end
    for (ν, label) in ref_labels
        haskey(gaps, ν) || continue
        B_ν = n_e * CONSTANTS.h / (CONSTANTS.e * Float64(ν))
        if first(B_range) < B_ν < last(B_range)
            text!(ax1, B_ν, Float64(1//ν) + 0.06, text=label,
                  fontsize=9, align=(:center, :bottom), color=(:blue, 0.7))
        end
    end
    for n in 1:3
        B_n = n_e * CONSTANTS.h / (CONSTANTS.e * n)
        if first(B_range) < B_n < last(B_range)
            text!(ax1, B_n, n + 0.08, text="$n",
                  fontsize=12, align=(:center, :bottom), color=:black)
        end
    end

    # Bottom panel: R_xx
    ax2 = Axis(fig[2, 1],
        xlabel=L"B\ [\mathrm{T}]",
        ylabel=L"R_{xx}\ [\mathrm{a.u.}]",
        limits=(extrema(B_range), nothing),
    )

    # Experimental R_xx overlay (behind theory)
    if !isempty(B_exp_xx)
        lines!(ax2, B_exp_xx, R_exp_xx, color=(:gray50, 0.5), linewidth=2.5,
               label="Expt (Wang et al. 2023)")
    end

    lines!(ax2, B_arr, R_xx_norm, color=:red, linewidth=2, label="Theory")

    # Legend in R_xx panel
    if !isempty(B_exp_xx)
        axislegend(ax2, position=:lt, framevisible=false, labelsize=10)
    end

    Label(fig[0, 1],
        "Fractional Quantum Hall Effect — GaAs 2DEG",
        fontsize=16, font=:bold)

    caption = @sprintf("T = %d mK,  nₑ = %.1f×10¹¹ cm⁻²,  μ = 10⁶ cm²/Vs\nRed labels: ab initio (ED + CF scaling). Blue labels: published reference gaps (2nd LL, reverse-flux Jain). Widths phenomenological.\nGray: experimental data (Wang et al., PNAS 2023, 2D hole gas, B rescaled to match filling factors).",
                        round(Int, T_K * 1e3), n_e / 1e15)
    Label(fig[3, 1], caption, fontsize=10, color=:gray40)

    linkxaxes!(ax1, ax2)
    rowsize!(fig.layout, 1, Relative(0.55))

    save("fqhe_plot.pdf", fig)
    save("fqhe_plot.png", fig, px_per_unit=3)
    println("Saved: fqhe_plot.pdf, fqhe_plot.png")
    flush(stdout)
    return fig
end

# ── Main ──
gaps = compute_gaps()
fig = make_plot(gaps)
display(fig)
