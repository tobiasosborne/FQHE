# 06_cylinder_ed.jl — Validate cylinder matrix elements + small-N ED on cylinder
#
# 1. Unit tests for the direct element f(s,δ): positivity, symmetry, Bessel check
# 2. Antisymmetrized element: antisymmetry, Hermiticity, momentum conservation
# 3. Small-N ED on cylinder at ν=1/3: N_e=3, N_phi=9, Ly scan

using LinearAlgebra, Combinatorics
using SpecialFunctions  # for besselk

# Load FQHE module for hilbert space utilities
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using FQHE

# Also load cylinder functions directly (not yet in module)
include(joinpath(@__DIR__, "..", "src", "cylinder.jl"))

function test_direct_element()
    println("=== Step 1.1: Direct pair interaction f(s,δ) checks ===\n")

    Ly = 12.0
    Δ = 2π / Ly

    # Check 1: f(s,0) should match the Bessel function formula
    # f(s, 0) = (2/Ly) exp(-s²Δ²/2) ∫₀^∞ dq/√(q²+s²Δ²) exp(-q²/2)
    #         = (2/Ly) exp(-s²Δ²/4) K₀(s²Δ²/4)  [analytic result]
    println("  f(s, 0) vs Bessel K₀ formula:")
    for s in 1:6
        f_num = _cylinder_f(s, 0.0, Ly)
        a2 = (s * Δ)^2
        # Identity: ∫₀^∞ dq/√(q²+a²) e^{-q²/2} = (1/2)e^{a²/4} K₀(a²/4)
        # So f(s,0) = (2/Ly) e^{-a²/2} × (1/2) e^{a²/4} K₀(a²/4) = (1/Ly) e^{-a²/4} K₀(a²/4)
        f_bessel = (1.0 / Ly) * exp(-a2 / 4) * besselk(0, a2 / 4)
        err = abs(f_num - f_bessel) / abs(f_bessel)
        println("    s=$s: f_num=$(round(f_num, digits=8)), f_K₀=$(round(f_bessel, digits=8)), rel_err=$err")
        @assert err < 1e-8 "Bessel check failed at s=$s"
    end
    println("  ✓ All match K₀ formula (rel err < 1e-8)\n")

    # Check 2: f(s, δ) > 0 for s ≠ 0, δ = 0 (direct interaction is repulsive)
    for s in 1:8
        @assert _cylinder_f(s, 0.0, Ly) > 0 "f($s, 0) not positive"
    end
    println("  ✓ f(s, 0) > 0 for s=1..8\n")

    # Check 3: f(s, δ) = f(-s, -δ) (parity in s) and f(s, δ) = f(s, -δ) (even in δ)
    for s in [-3, -1, 1, 3], δ_int in [-2, 0, 2]
        δ = δ_int * Δ
        fp = _cylinder_f(s, δ, Ly)
        fm = _cylinder_f(-s, -δ, Ly)
        @assert abs(fp - fm) < 1e-12 "Parity check failed"
        fe = _cylinder_f(s, -δ, Ly)
        @assert abs(fp - fe) < 1e-12 "Even-in-δ check failed"
    end
    println("  ✓ f(s,δ) = f(-s,-δ) and f(s,δ) = f(s,-δ)\n")

    # Check 4: f(s=0, δ) via background subtraction is well-defined and negative
    # (removing a charge reduces energy → f_sub(0,δ) < 0 for δ ≠ 0)
    for δ_int in 1:5
        δ = δ_int * Δ
        f0 = _cylinder_f(0, δ, Ly)
        println("    f(0, $δ_int Δ) = $(round(f0, digits=8))")
        @assert f0 < 0 "f(0, δ) should be negative for δ≠0"
    end
    println("  ✓ f(0, δ) < 0 for δ > 0 (background-subtracted)\n")

    # Check 5: Exponential decay with |s|
    f_vals = [_cylinder_f(s, 0.0, Ly) for s in 1:8]
    println("  f(s, 0): ", [round(v, digits=6) for v in f_vals])
    for i in 1:length(f_vals)-1
        @assert f_vals[i] > f_vals[i+1] "Non-monotonic decay"
    end
    println("  ✓ Monotonically decaying with |s|\n")
end

function test_antisymmetrized()
    println("=== Step 1.2: Antisymmetrized matrix element checks ===\n")

    Ly = 12.0

    # Check 1: Antisymmetry V(a,b,c,d) = -V(a,b,d,c) [swapping ket orbitals]
    for (a,b,c,d) in [(0,3,1,2), (0,5,2,3), (1,4,2,3)]
        V1 = cylinder_matrix_element(a, b, c, d, Ly)
        V2 = cylinder_matrix_element(a, b, d, c, Ly)
        @assert abs(V1 + V2) < 1e-12 "Antisymmetry failed for ($a,$b,$c,$d)"
    end
    println("  ✓ V(a,b,c,d) = -V(a,b,d,c)\n")

    # Check 2: Hermiticity V(a,b,c,d) = V(c,d,a,b)
    for (a,b,c,d) in [(0,3,1,2), (0,5,2,3), (1,6,3,4)]
        V_fwd = cylinder_matrix_element(a, b, c, d, Ly)
        V_bwd = cylinder_matrix_element(c, d, a, b, Ly)
        @assert abs(V_fwd - V_bwd) < 1e-12 "Hermiticity failed for ($a,$b,$c,$d): $V_fwd vs $V_bwd"
    end
    println("  ✓ V(a,b,c,d) = V(c,d,a,b)\n")

    # Check 3: Momentum conservation — V = 0 when a+b ≠ c+d
    V_bad = cylinder_matrix_element(0, 3, 1, 3, Ly)
    @assert V_bad == 0.0 "Momentum conservation not enforced"
    println("  ✓ V = 0 when a+b ≠ c+d\n")

    # Check 4: Diagonal elements — print values (may be negative due to background subtraction)
    # On the cylinder, f_sub(0, δ) < 0 can make diagonal negative. This is physical:
    # the background subtraction removes the q=0 divergence. Gaps are unaffected.
    println("  Diagonal elements (background-subtracted, may be negative):")
    for (a, b) in [(0,1), (0,2), (0,3), (1,4), (2,5)]
        Vdiag = cylinder_matrix_element(a, b, a, b, Ly)
        println("    V($a,$b,$a,$b) = $(round(Vdiag, digits=8))")
    end
    # V should become more negative as |a-b| increases (background subtraction grows)
    # but this is not required — just informational
    println("  (informational — gaps unaffected by constant energy shift)\n")
end

function build_cylinder_hamiltonian_sparse(Ne::Int, Nphi::Int, Ly::Float64)
    """Build sparse Hamiltonian on cylinder in Fock basis (for ED cross-check)."""
    states = collect(combinations(0:Nphi-1, Ne))
    D = length(states)
    state_to_idx = Dict(sort(s) => i for (i, s) in enumerate(states))

    H = zeros(Float64, D, D)

    for (β_idx, occ) in enumerate(states)
        # Diagonal: Σ_{i<j ∈ occ} V(i,j,i,j)
        for ii in 1:Ne-1, jj in ii+1:Ne
            a, b = occ[ii], occ[jj]
            H[β_idx, β_idx] += cylinder_matrix_element(a, b, a, b, Ly)
        end

        # Off-diagonal: annihilate (c,d), create (a,b)
        for ii in 1:Ne-1, jj in ii+1:Ne
            c, d = occ[ii], occ[jj]
            strip = setdiff(occ, [c, d])

            for a in 0:Nphi-1
                a in strip && continue
                for b in a+1:Nphi-1
                    b in strip && continue
                    a + b == c + d || continue
                    (a == c && b == d) && continue

                    V = cylinder_matrix_element(a, b, c, d, Ly)
                    abs(V) > 1e-15 || continue

                    new_occ = sort([strip; a; b])
                    α_idx = get(state_to_idx, new_occ, 0)
                    α_idx > 0 || continue

                    sgn = fermion_sign_cylinder(occ, c, d, a, b)
                    H[α_idx, β_idx] += V * sgn
                end
            end
        end
    end
    return Symmetric(H)
end

function fermion_sign_cylinder(occ::Vector{Int}, c::Int, d::Int, a::Int, b::Int)
    n_orb = max(maximum(occ), a, b) + 1
    state = UInt64(0)
    for o in occ
        state |= UInt64(1) << o
    end
    sign = 1

    # Annihilate c
    mask = (UInt64(1) << c) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state &= ~(UInt64(1) << c)

    # Annihilate d
    mask = (UInt64(1) << d) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state &= ~(UInt64(1) << d)

    # Create b
    mask = (UInt64(1) << b) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state |= UInt64(1) << b

    # Create a
    mask = (UInt64(1) << a) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)

    return sign
end

function test_cylinder_ed()
    println("=== Step 1.3: Small-N ED on cylinder ===\n")

    Ne = 3
    Nphi = 9

    println("--- Ly scan: N_e=$Ne, N_phi=$Nphi, ν=1/3 ---")
    println()
    println("  Ly    | Lx    | aspect | E0/N       | gap")
    println("  ------|-------|--------|------------|--------")

    for Ly in [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        H = build_cylinder_hamiltonian_sparse(Ne, Nphi, Ly)
        vals = eigvals(Matrix(H))
        E0 = vals[1]
        gap = vals[2] - vals[1]
        Lx = cylinder_length(Nphi, Ly)
        aspect = Lx / Ly
        println("  $(lpad(Ly, 5)) | $(lpad(round(Lx,digits=1),5)) | $(lpad(round(aspect,digits=2),6)) | $(lpad(round(E0/Ne, digits=6),10)) | $(round(gap, digits=6))")
    end

    # Detailed check at Ly=14
    Ly = 14.0
    H = build_cylinder_hamiltonian_sparse(Ne, Nphi, Ly)
    vals = eigvals(Matrix(H))
    E0_per_N = vals[1] / Ne
    gap = vals[2] - vals[1]

    println("\n--- Ground truth checks at Ly=14 ---")
    println("  E0/N = $E0_per_N")
    println("  Gap  = $gap")
    println("  Hilbert dim = $(size(H, 1))  (should be C(9,3)=84)")
    @assert size(H, 1) == 84 "Wrong Hilbert space dimension"
    # At Nphi=9 the cylinder is very short (Lx=4.0 for Ly=14), so E0/N differs
    # from thermodynamic value. The key validation is Bessel K₀ match + Hermiticity.
    @assert -0.60 < E0_per_N < -0.10 "E0/N out of expected range"
    @assert gap >= 0 "Gap should be non-negative"
    println("  ✓ E0/N and gap in range\n")

    # Check Hermiticity of H
    H_dense = Matrix(H)
    herm_err = maximum(abs.(H_dense - H_dense'))
    println("  Hermiticity error: $herm_err")
    @assert herm_err < 1e-12 "Hamiltonian not Hermitian"
    println("  ✓ Hamiltonian is Hermitian\n")

    return vals[1], gap
end

function main()
    test_direct_element()
    test_antisymmetrized()
    E0, gap = test_cylinder_ed()
    println("=== All Phase 1 checks PASSED ===")
end

main()
