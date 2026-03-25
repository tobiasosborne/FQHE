# 02_compute_gaps.jl — Ab initio ED for ALL FQHE fractions
# Run: julia --project=. scripts/02_compute_gaps.jl
#
# Computes neutral gaps for 13 LLL fractions via exact diagonalisation
# on the Haldane sphere. Fast mode: dim < 5000, total ~10 min.
# Results saved to data/ed_gaps.csv for reuse by plot script.
#
# Sources:
#   P04: Fano et al., PRB 34, 2670 (1986)
#   P14: Balram et al., PRB (various)

using FQHE
using Printf

# ── All LLL fractions for ab initio ED ──
# N values chosen so: (1) twoS = N/ν - q is a positive integer
#                      (2) Hilbert space dim(Lz=0) < 5000
const ED_FRACTIONS = [
    # Jain principal p/(2p+1)
    1//3, 2//5, 3//7, 4//9, 5//11,
    # Laughlin 1/m
    1//5, 1//7, 1//9,
    # 4-flux reverse p/(4p-1)
    2//7, 3//11, 4//15,
    # 4-flux forward p/(4p+1)
    2//9, 3//13,
]

const N_RANGE = Dict(
    # Jain principal
    1//3  => [3, 4, 5, 6, 7],     # dim: 5, 21, 99, 500, 2652
    2//5  => [4, 6, 8],            # dim: 2, 42, 804
    3//7  => [6],                  # dim: 3
    4//9  => [8],                  # dim: 4
    5//11 => [10],                 # dim: ~5
    # Laughlin
    1//5  => [3, 4, 5],            # dim: 15, 113, 969 (N=6→8855, skip)
    1//7  => [3, 4, 5],            # dim: 30, 332, 4095
    1//9  => [3, 4],               # dim: 51, 731 (N=5→11781, skip)
    # 4-flux reverse
    2//7  => [4, 6],               # dim: 8, 333 (N=8→14535, skip)
    3//11 => [6],                  # dim: 77
    4//15 => [8],                  # dim: 804
    # 4-flux forward
    2//9  => [4, 6],               # dim: 21, 1428
    3//13 => [6],                  # dim: 214
)

# Published reference gaps for comparison (neutral gaps, e²/(ε ℓ_B))
const REFERENCE_GAPS = Dict(
    1//3  => 0.1036,   # Fano 1986
    2//5  => 0.058,    # Balram 2018
    3//7  => 0.030,    # Balram 2018
    1//5  => 0.0244,   # Fano 1986
    1//7  => 0.012,    # estimated
    2//7  => 0.008,    # Balram 2015
    3//11 => 0.003,    # Balram 2015
)

function run_all_gaps()
    results = Dict{Rational{Int}, NamedTuple}()

    for ν in ED_FRACTIONS
        Ns = N_RANGE[ν]
        println("\n" * "="^60)
        println("ν = $ν")
        println("="^60)

        n_gaps = Float64[]

        for N in Ns
            twoS = sphere_flux(N, ν)
            basis = enumerate_fock_states(N, twoS; twoLz=0)
            dim = length(basis.states)

            # Skip if too large
            if dim > 5000
                @printf("  N=%2d  2S=%2d  dim=%6d  SKIPPED (dim > 5000)\n", N, twoS, dim)
                continue
            end

            t = @elapsed begin
                E0_pp, Δ_n = neutral_gap(basis, coulomb_pseudopotentials(twoS))
            end
            push!(n_gaps, Δ_n)
            @printf("  N=%2d  2S=%2d  dim=%6d  Δ_neutral=%8.6f  (%.1fs)\n",
                    N, twoS, dim, Δ_n, t)
        end

        # Filter out N values that were skipped
        Ns_used = Ns[1:length(n_gaps)]

        # Extrapolate
        Δn_inf = length(n_gaps) >= 2 ? extrapolate_gap(Ns_used, n_gaps) : n_gaps[end]
        ref = get(REFERENCE_GAPS, ν, NaN)
        ratio_str = isnan(ref) ? "—" : @sprintf("%.0f%%", Δn_inf / ref * 100)
        @printf("  → Δ(∞) = %.6f  (ref: %s, ratio: %s)\n",
                Δn_inf, isnan(ref) ? "—" : @sprintf("%.4f", ref), ratio_str)

        results[ν] = (Ns=Ns_used, gaps=n_gaps, gap_inf=Δn_inf)
    end

    return results
end

function save_results(results)
    outdir = joinpath(@__DIR__, "..", "data")
    mkpath(outdir)
    outfile = joinpath(outdir, "ed_gaps.csv")

    open(outfile, "w") do io
        println(io, "nu_p,nu_q,N_values,gaps,gap_extrapolated")
        for ν in sort(collect(keys(results)); by=x->(denominator(x), numerator(x)))
            r = results[ν]
            p, q = numerator(ν), denominator(ν)
            ns_str = join(r.Ns, ";")
            gs_str = join([@sprintf("%.8f", g) for g in r.gaps], ";")
            @printf(io, "%d,%d,%s,%s,%.8f\n", p, q, ns_str, gs_str, r.gap_inf)
        end
    end
    println("\nSaved to $outfile")
end

function print_summary(results)
    println("\n" * "="^70)
    println("SUMMARY: Ab initio neutral gaps (e²/(ε ℓ_B))")
    println("="^70)
    @printf("%-7s  %-8s  %-12s  %-10s  %-8s  %s\n",
            "ν", "N range", "Δ_neutral(∞)", "Reference", "Ratio", "Method")
    println("-"^70)

    for ν in ED_FRACTIONS
        haskey(results, ν) || continue
        r = results[ν]
        ref = get(REFERENCE_GAPS, ν, NaN)
        ratio = isnan(ref) ? "—" : @sprintf("%.0f%%", r.gap_inf / ref * 100)
        ns_str = "$(r.Ns[1])-$(r.Ns[end])"
        @printf("%-7s  %-8s  %10.6f    %-10s  %-8s  ED\n",
                ν, ns_str, r.gap_inf,
                isnan(ref) ? "—" : @sprintf("%.4f", ref), ratio)
    end

    # P-H conjugates
    println("\nParticle-hole conjugates (Δ(ν) = Δ(1−ν)):")
    for ν in sort(collect(keys(results)); by=x->Float64(x))
        ν_conj = (denominator(ν) - numerator(ν)) // denominator(ν)
        ν_conj != ν && 0 < ν_conj < 1 || continue
        @printf("  Δ(%s) = Δ(%s) = %.6f\n", ν_conj, ν, results[ν].gap_inf)
    end
end

# ── Main ──
results = run_all_gaps()
print_summary(results)
save_results(results)
