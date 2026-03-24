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

    # ── Step 1: ν=1/3 gap from precomputed ED (N=3..8) ──
    # Precomputed via numerical CG + Lanczos (N=8 takes ~18 min).
    # To recompute: julia --project=. scripts/02_compute_gaps.jl
    Ns_13 = [3, 4, 5, 6, 7, 8]
    ng_13 = [0.11899158, 0.09353321, 0.09311677, 0.08203499, 0.08082519, 0.08159550]
    Δ_13 = extrapolate_gap(Ns_13, ng_13; order=2)
    @printf("ν=1/3: Δ(∞) = %.6f  (quadratic, N=3–8, ref: 0.1036, ratio: %.1f%%)\n", Δ_13, Δ_13/0.1036*100)
    flush(stdout)
    gaps[1//3] = Δ_13

    # ── Step 2: ν=1/5 gap from precomputed ED (N=3..5) ──
    Ns_15 = [3, 4, 5]
    ng_15 = [0.02770232, 0.01904592, 0.02062640]
    Δ_15 = extrapolate_gap(Ns_15, ng_15; order=1)
    @printf("ν=1/5: Δ(∞) = %.6f  (linear, N=3–5, ref: 0.0244)\n", Δ_15)
    flush(stdout)
    gaps[1//5] = Δ_15

    # ── Step 3: Jain principal sequence via CF scaling ──
    # Geometric decay: Δ(p/(2p+1)) = Δ(1/3) × 0.56^(p-1)
    # Calibrated from published Δ(2/5)/Δ(1/3), validated at p=3,4.
    println("Jain hierarchy from CF scaling (anchored to ED 1/3 gap)...")
    flush(stdout)
    for p in 2:7
        ν = p//(2p+1)
        Δ = cf_gap_scaling(Δ_13, p)
        gaps[ν] = Δ
        @printf("  ν=%-5s (p=%d)  Δ=%.6f\n", ν, p, Δ)
        flush(stdout)
    end

    # ── Step 4: Particle-hole conjugates (LLL) ──
    ph_map = [(2//3, 1//3), (3//5, 2//5), (4//7, 3//7),
              (5//9, 4//9), (6//11, 5//11), (7//13, 6//13)]
    for (ν_conj, ν_orig) in ph_map
        haskey(gaps, ν_orig) || continue
        gaps[ν_conj] = gaps[ν_orig]
    end
    println("  LLL conjugates: Δ(ν) = Δ(1−ν)")
    flush(stdout)

    # ── Step 5: Second Landau level fractions (PUBLISHED REFERENCE GAPS) ──
    # These are NOT computed by our pipeline. They use effective 2nd-LL
    # pseudopotentials which differ from LLL Coulomb. Values from:
    #   ν=5/2:  Morf PRB 58 R1465 (1998); DMRG Zaletel et al. PRB 91 (2015)
    #   ν=7/3:  Morf & d'Ambrumenil PRL 74 5116 (1995)
    #   ν=12/5: Rezayi & Read PRB 79 075306 (2009)
    # TODO: replace with own DMRG calculations via ITensors.jl
    println("Second LL fractions (published reference gaps)...")
    flush(stdout)
    sll_gaps = [
        (5//2,  0.025),   # Moore-Read / Pfaffian state
        (7//3,  0.045),   # Laughlin-like in 2nd LL
        (8//3,  0.045),   # p-h conjugate of 7/3
        (12//5, 0.015),   # Jain 2/5 in 2nd LL
        (13//5, 0.015),   # p-h conjugate of 12/5
    ]
    for (ν, Δ) in sll_gaps
        gaps[ν] = Δ
        @printf("  ν=%-5s  Δ=%.4f  [published ref]\n", ν, Δ)
        flush(stdout)
    end

    # ── Step 6: Reverse-flux Jain fractions (LLL, published refs) ──
    # ν = p/(2p−1) sequence: weaker states, very small gaps.
    #   Balram et al. PRL 115 186805 (2015)
    # TODO: compute via ED once second-LL machinery is in place
    println("Reverse-flux Jain fractions (published reference gaps)...")
    flush(stdout)
    rev_jain = [
        (2//7,  0.008),   # reverse-flux p=2
        (3//11, 0.003),   # reverse-flux p=3
    ]
    for (ν, Δ) in rev_jain
        gaps[ν] = Δ
        @printf("  ν=%-5s  Δ=%.4f  [published ref]\n", ν, Δ)
        flush(stdout)
    end

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

    # Digitized Störmer Fig. 18 R_xy overlay (behind theory)
    stormer_path = joinpath(@__DIR__, "..", "data", "stormer_rxy_rescaled.csv")
    if isfile(stormer_path)
        sd = readdlm(stormer_path, ',', skipstart=1)
        lines!(ax1, sd[:, 1], sd[:, 2], color=(:gray50, 0.6), linewidth=3,
               label="Störmer (RMP 1999, rescaled)")
    end

    # Wang et al. R_xy overlay (1/3 plateau region only)
    if !isempty(B_exp_xy)
        lines!(ax1, B_exp_xy, R_exp_xy, color=(:orange, 0.7), linewidth=2.5,
               label="Wang et al. (PNAS 2023)")
    end

    lines!(ax1, B_arr, R_xy_norm, color=:blue, linewidth=2, label="Theory (ED + CF)")
    axislegend(ax1, position=:lt, framevisible=false, labelsize=9)

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

    caption = @sprintf("T = %d mK,  nₑ = %.1f×10¹¹ cm⁻²,  μ = 10⁶ cm²/Vs\nRed labels: ab initio (ED + CF). Blue labels: published reference gaps (2nd LL). Gray: Störmer RMP 1999 Fig. 18 (digitized, B rescaled).\nOrange: Wang et al. PNAS 2023 (2D hole gas). R_xx: red = theory, gray = experiment.",
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
