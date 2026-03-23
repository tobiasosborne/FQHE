# 02_compute_gaps.jl — Run ED for all target fractions
# Run: julia --project=. scripts/02_compute_gaps.jl
#
# Computes charge gaps for FQHE fractions via exact diagonalisation
# on the Haldane sphere with Coulomb pseudopotentials (Fano eq. 25).

using FQHE
using Printf

const TARGET_FRACTIONS = [1//3, 2//5, 3//7, 1//5, 2//3, 3//5]
const N_RANGE = Dict(
    1//3 => [3, 4, 5, 6, 7, 8],
    2//5 => [4, 6, 8],
    3//7 => [6, 9],
    1//5 => [3, 4, 5, 6],
    2//3 => [4, 6, 8],
    3//5 => [4, 6, 8],
)

function run_all_gaps()
    results = Dict{Rational{Int}, GapTable}()

    for ν in TARGET_FRACTIONS
        Ns = N_RANGE[ν]
        E0s = Float64[]
        gaps = Float64[]

        println("\n=== ν = $ν ===")
        for N in Ns
            twoS = sphere_flux(N, ν)
            @printf("  N=%2d  2S=%2d  ", N, twoS)

            VJ = coulomb_pseudopotentials(twoS)
            basis = enumerate_fock_states(N, twoS; twoLz=0)
            dim = length(basis.states)
            @printf("dim=%6d  ", dim)

            t = @elapsed begin
                E0_pp, Δ = charge_gap(N, ν)
            end
            push!(E0s, E0_pp)
            push!(gaps, Δ)

            @printf("E₀/N = %10.6f  Δ̃ = %8.6f  (%.1fs)\n", E0_pp, Δ, t)
        end

        gap_inf = extrapolate_gap(Ns, gaps)
        @printf("  Δ̃(∞) = %.6f\n", gap_inf)

        results[ν] = GapTable(ν, Ns, E0s, gaps, gap_inf)
    end

    # Summary table
    println("\n" * "="^60)
    println("Summary: extrapolated gaps in units of e²/(ε ℓ_B)")
    println("="^60)
    @printf("%-8s  %-12s  %-12s\n", "ν", "Δ̃(∞)", "Fano ref")
    println("-"^40)
    refs = Dict(1//3 => 0.1036, 1//5 => 0.0244, 2//5 => NaN,
                3//7 => NaN, 2//3 => NaN, 3//5 => NaN)
    for ν in TARGET_FRACTIONS
        ref = get(refs, ν, NaN)
        @printf("%-8s  %10.6f    %s\n", ν, results[ν].gap_extrap,
                isnan(ref) ? "—" : @sprintf("%.4f", ref))
    end

    return results
end

results = run_all_gaps()
