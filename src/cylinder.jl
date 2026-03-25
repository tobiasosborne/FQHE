# cylinder.jl — Cylinder geometry and Coulomb matrix elements in LLL
#
# Uses the Haldane PSEUDOPOTENTIAL decomposition: H = Σ_m V_m P_m
# where V_m are geometry-independent Haldane pseudopotentials and P_m
# are pair projectors at relative angular momentum m.
#
# On the cylinder, the pair projector involves Hermite functions of the
# relative orbital index:
#   ⟨a,b|P_m|c,d⟩ = δ_{K,K'} g_m(a-b) g_m(c-d)
#   g_m(r) = N_m H_m(rΔ/√2) exp(-r²Δ²/4)
#
# This approach gives positive-definite matrix elements (V_m > 0) and
# avoids all background/Madelung issues.
#
# Sources:
#   P01: Haldane, PRL 51, 605 (1983) — pseudopotentials
#   P28: Rezayi & Haldane, PRL 84, 4685 (2000) — cylinder geometry
#   Seidel & Lee, PRL 97, 056804 (2006) — thin cylinder limit

export cylinder_matrix_element, cylinder_orbital_spacing, cylinder_length
export plane_pseudopotentials

using Combinatorics: binomial
using SpecialFunctions: loggamma

"""
    cylinder_orbital_spacing(Ly) → Δ

Spacing between adjacent LLL orbital guiding centers on a cylinder
of circumference `Ly` (in units of ℓ_B).
"""
cylinder_orbital_spacing(Ly::Real) = 2π / Ly

"""
    cylinder_length(Nphi, Ly) → L_x

Physical length of a cylinder hosting `Nphi` LLL orbitals
at circumference `Ly` (both in units of ℓ_B).
"""
cylinder_length(Nphi::Integer, Ly::Real) = Nphi * 2π / Ly

"""
    plane_pseudopotentials(mmax) → Vector{Float64}

Haldane pseudopotentials V_m for the Coulomb interaction 1/|r| on the
infinite plane, projected into the LLL.

  V_m = C(2m, m) / 4^m

Returns V[m+1] for m = 0, 1, ..., mmax. Units: e²/(ε ℓ_B).
For fermions, only odd m contribute (even m are zero by Pauli exclusion).
"""
function plane_pseudopotentials(mmax::Int)
    Vm = zeros(mmax + 1)
    for m in 0:mmax
        Vm[m+1] = binomial(2m, m) / 4.0^m
    end
    return Vm
end

"""
    _hermite(n, x) → Float64

Physicists' Hermite polynomial H_n(x) via recurrence.
"""
function _hermite(n::Int, x::Float64)
    n == 0 && return 1.0
    n == 1 && return 2x
    H0, H1 = 1.0, 2x
    for k in 2:n
        H0, H1 = H1, 2x * H1 - 2(k - 1) * H0
    end
    return H1
end

"""
    _pair_amplitude(m, r, Ly) → Float64

Normalized pair amplitude g_m(r) on the cylinder:
  g_m(r) = √(Δ / (√(2π) 2^m m!)) × H_m(rΔ/√2) × exp(-r²Δ²/4)

For odd m: g_m(-r) = -g_m(r) (antisymmetric, as required for fermions).
For even m: g_m(-r) = g_m(r) (symmetric, excluded by Pauli for identical fermions).
"""
function _pair_amplitude(m::Int, r::Int, Ly::Float64)
    Δ = 2π / Ly
    x = r * Δ / sqrt(2.0)
    # Normalization: N_m = √(Δ / (√(2π) × 2^m × m!))
    # Use log to avoid overflow for large m
    log_norm = 0.5 * (log(Δ) - 0.5 * log(2π) - m * log(2.0) - loggamma(m + 1))
    N_m = exp(log_norm)
    H_m = _hermite(m, x)
    gauss = exp(-r^2 * Δ^2 / 4)
    return N_m * H_m * gauss
end

"""
    cylinder_matrix_element(a, b, c, d, Ly; mmax=30) → Float64

Two-body Coulomb matrix element ⟨ab|V|cd⟩ on the cylinder via pseudopotential
decomposition.

  V_{abcd} = δ_{a+b,c+d} × Σ_{m odd} V_m × g_m(a-b) × g_m(c-d)

0-indexed orbitals. Units: e²/(ε ℓ_B).
"""
function cylinder_matrix_element(a::Integer, b::Integer, c::Integer, d::Integer,
                                  Ly::Real; mmax::Int=30)
    a + b == c + d || return 0.0

    Vm = plane_pseudopotentials(mmax)
    r_bra = a - b
    r_ket = c - d

    val = 0.0
    for m in 1:2:mmax   # odd m only (fermions)
        gm_bra = _pair_amplitude(m, r_bra, Float64(Ly))
        gm_ket = _pair_amplitude(m, r_ket, Float64(Ly))
        val += Vm[m+1] * gm_bra * gm_ket
    end
    return val
end
