# composite_fermion.jl — Jain CF theory: gap predictions, CF filling map
#
# Sources:
#   P3: Jain, PRL 63, 199 (1989)
#   P8: Jain, Composite Fermions, Ch. 7

export cf_effective_field, cf_gap_estimate

"""
    cf_effective_field(B, n_e) → B* [T]

Effective magnetic field seen by composite fermions: B* = B - 2n_e φ₀.
"""
function cf_effective_field(B::Real, n_e::Real)
    φ₀ = CONSTANTS.h / CONSTANTS.e  # flux quantum
    return B - 2 * n_e * φ₀
end

"""
    cf_gap_estimate(ν, B, mat) → Δ [J]

Composite fermion gap estimate for Jain fraction ν = p/(2p+1).
Δ_CF ≈ ℏeB*/(m_CF) where m_CF ≈ C × √B × m_e.
"""
function cf_gap_estimate(ν::Rational, B::Real, mat::MaterialParams=GaAs)
    # TODO: implement with proper CF effective mass from P8 Ch. 7
    error("Not yet implemented")
end
