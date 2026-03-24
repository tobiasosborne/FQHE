# 04_make_plot.jl — Final two-panel FQHE figure
# Run: julia --project=. scripts/04_make_plot.jl
#
# Output: fqhe_plot.pdf, fqhe_plot.png
#
# ν=1/3 gap from ab initio ED (charge gap, extrapolated N→∞).
# Jain hierarchy from CF geometric scaling anchored to ED.
# Plateau widths and transition shapes are semi-phenomenological.

using FQHE
using CairoMakie
using Printf

# ── Sample parameters (close to original Tsui-Störmer experiment) ──
const n_e = 1.0e15      # m⁻² (= 1.0×10¹¹ cm⁻²)
const T_K = 50e-3        # 50 mK
const μ_e = 1e2          # m²/(V·s) = 10⁶ cm²/(V·s)
const B_range = range(0.3, 18.0, length=4000)

function compute_gaps()
    gaps = Dict{Rational{Int}, Float64}()

    # ── Step 1: Charge gap for ν=1/3 from ED (Laughlin state) ──
    println("Computing ν=1/3 charge gap from ED...")
    flush(stdout)
    Ns_13 = [3, 4, 5, 6]
    cg_13 = Float64[]
    for N in Ns_13
        _, Δ = charge_gap(N, 1//3)
        push!(cg_13, Δ)
        @printf("  N=%d  Δ_charge=%.6f\n", N, Δ)
        flush(stdout)
    end
    Δ_13 = extrapolate_gap(Ns_13, cg_13)
    @printf("  Δ(1/3, ∞) = %.6f  (ref: 0.1036)\n", Δ_13)
    flush(stdout)
    gaps[1//3] = Δ_13

    # ── Step 2: Charge gap for ν=1/5 from ED (Laughlin state) ──
    println("Computing ν=1/5 charge gap from ED...")
    flush(stdout)
    Ns_15 = [3, 4]
    cg_15 = Float64[]
    for N in Ns_15
        _, Δ = charge_gap(N, 1//5)
        push!(cg_15, Δ)
        @printf("  N=%d  Δ_charge=%.6f\n", N, Δ)
        flush(stdout)
    end
    Δ_15 = extrapolate_gap(Ns_15, cg_15)
    @printf("  Δ(1/5, ∞) = %.6f  (ref: 0.0244)\n", Δ_15)
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

    # ── Step 4: Particle-hole conjugates ──
    ph_map = [(2//3, 1//3), (3//5, 2//5), (4//7, 3//7),
              (5//9, 4//9), (6//11, 5//11), (7//13, 6//13)]
    for (ν_conj, ν_orig) in ph_map
        haskey(gaps, ν_orig) || continue
        gaps[ν_conj] = gaps[ν_orig]
    end
    println("  Conjugates: Δ(ν) = Δ(1−ν)")
    flush(stdout)

    return gaps
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
    lines!(ax1, B_arr, R_xy_norm, color=:blue, linewidth=2)

    # Annotate fractions
    frac_labels = [
        (1//3, "1/3"), (2//5, "2/5"), (3//7, "3/7"), (4//9, "4/9"), (5//11, "5/11"),
        (6//13, "6/13"), (7//15, "7/15"),
        (2//3, "2/3"), (3//5, "3/5"), (4//7, "4/7"), (5//9, "5/9"),
        (6//11, "6/11"), (7//13, "7/13"),
        (1//5, "1/5"),
    ]
    for (ν, label) in frac_labels
        haskey(gaps, ν) || continue
        B_ν = n_e * CONSTANTS.h / (CONSTANTS.e * Float64(ν))
        if first(B_range) < B_ν < last(B_range)
            text!(ax1, B_ν, Float64(1//ν) + 0.06, text=label,
                  fontsize=9, align=(:center, :bottom), color=:red)
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
    lines!(ax2, B_arr, R_xx_norm, color=:red, linewidth=2)

    Label(fig[0, 1],
        "Fractional Quantum Hall Effect — GaAs 2DEG",
        fontsize=16, font=:bold)

    caption = @sprintf("T = %d mK,  nₑ = %.1f×10¹¹ cm⁻²,  μ = 10⁶ cm²/Vs\n1/3 gap from ED (charge gap N→∞). Jain sequence (p=1–7) via CF scaling. Widths phenomenological.",
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
