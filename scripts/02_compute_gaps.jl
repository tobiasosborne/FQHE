# 02_compute_gaps.jl — Run ED for all target fractions
# Run: julia --project=. scripts/02_compute_gaps.jl
#
# Computes neutral gaps AND charge gaps for FQHE fractions via exact
# diagonalisation on the Haldane sphere with Coulomb pseudopotentials.
# Extrapolates to N→∞ via Δ(N) = Δ(∞) + a/N.
#
# Sources:
#   P04: Fano et al., PRB 34, 2670 (1986)
#   P14: Balram et al., PRB (various)

using FQHE
using Printf

# Fractions computed directly via ED
const ED_FRACTIONS = [1//3, 2//5, 3//7, 1//5]

# N values for each fraction (chosen so twoS is integer and Hilbert space is tractable)
# N values capped so charge gap Hilbert spaces stay ≤ ~13K
const N_RANGE = Dict(
    1//3 => [3, 4, 5, 6, 7, 8],
    2//5 => [4, 6, 8],
    3//7 => [6, 9],
    1//5 => [3, 4, 5, 6],
)

# Published reference values for comparison (charge gaps, e²/(ε ℓ_B))
const REFERENCE_GAPS = Dict(
    1//3 => 0.1036,   # Fano 1986 p.2678
    2//5 => 0.058,    # Balram 2018
    3//7 => 0.030,    # Balram 2018
    1//5 => 0.0244,   # Fano 1986 Table VIII
)

function run_all_gaps()
    neutral_results = Dict{Rational{Int}, NamedTuple}()
    charge_results  = Dict{Rational{Int}, NamedTuple}()

    for ν in ED_FRACTIONS
        Ns = N_RANGE[ν]
        println("\n" * "="^60)
        println("ν = $ν")
        println("="^60)

        n_gaps = Float64[]
        c_gaps = Float64[]
        c_E0s  = Float64[]

        for N in Ns
            twoS = sphere_flux(N, ν)

            # --- Neutral gap ---
            basis = enumerate_fock_states(N, twoS; twoLz=0)
            dim = length(basis.states)
            t_n = @elapsed begin
                E0_pp, Δ_n = neutral_gap(basis, coulomb_pseudopotentials(twoS))
            end
            push!(n_gaps, Δ_n)
            @printf("  N=%2d  2S=%2d  dim=%6d  Δ_neutral=%8.6f  (%.1fs)\n",
                    N, twoS, dim, Δ_n, t_n)

            # --- Charge gap ---
            t_c = @elapsed begin
                E0_c, Δ_c = charge_gap(N, ν)
            end
            push!(c_gaps, Δ_c)
            push!(c_E0s, E0_c)
            @printf("  %23s  Δ_charge=%8.6f  (%.1fs)\n", "", Δ_c, t_c)
        end

        # Extrapolate
        Δn_inf = extrapolate_gap(Ns, n_gaps)
        Δc_inf = extrapolate_gap(Ns, c_gaps)
        @printf("  Extrapolated:  Δ_neutral(∞)=%.6f  Δ_charge(∞)=%.6f  (ref=%.4f)\n",
                Δn_inf, Δc_inf, get(REFERENCE_GAPS, ν, NaN))

        neutral_results[ν] = (Ns=Ns, gaps=n_gaps, gap_inf=Δn_inf)
        charge_results[ν]  = (Ns=Ns, gaps=c_gaps, E0s=c_E0s, gap_inf=Δc_inf)
    end

    # Particle-hole conjugates
    println("\n" * "="^60)
    println("Particle-hole conjugates (Δ(ν) ≈ Δ(1−ν))")
    println("="^60)
    for (ν, ν_conj) in [(2//3, 1//3), (3//5, 2//5)]
        @printf("  Δ(%s) = Δ(%s) = %.6f\n", ν, ν_conj, charge_results[ν_conj].gap_inf)
    end

    # Summary
    println("\n" * "="^60)
    println("Summary: extrapolated gaps in units of e²/(ε ℓ_B)")
    println("="^60)
    @printf("%-6s  %-12s  %-12s  %-10s  %-10s\n",
            "ν", "Δ_charge(∞)", "Δ_neutral(∞)", "Reference", "Ratio")
    println("-"^60)
    for ν in ED_FRACTIONS
        ref = get(REFERENCE_GAPS, ν, NaN)
        Δc = charge_results[ν].gap_inf
        Δn = neutral_results[ν].gap_inf
        ratio = isnan(ref) ? NaN : Δc / ref
        @printf("%-6s  %10.6f    %10.6f    %10s  %8s\n",
                ν, Δc, Δn,
                isnan(ref) ? "—" : @sprintf("%.4f", ref),
                isnan(ratio) ? "—" : @sprintf("%.3f", ratio))
    end
    for (ν, ν_conj) in [(2//3, 1//3), (3//5, 2//5)]
        Δc = charge_results[ν_conj].gap_inf
        ref = get(REFERENCE_GAPS, ν_conj, NaN)
        @printf("%-6s  %10.6f    %10s    %10s  (= %s)\n",
                ν, Δc, "—", isnan(ref) ? "—" : @sprintf("%.4f", ref), ν_conj)
    end

    return neutral_results, charge_results
end

neutral_results, charge_results = run_all_gaps()
