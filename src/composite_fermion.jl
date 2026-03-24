# composite_fermion.jl — Jain CF theory: gap predictions, CF filling map
#
# Sources:
#   P3: Jain, PRL 63, 199 (1989)
#   P8: Jain, Composite Fermions, Ch. 7

export cf_effective_field, cf_gap_scaling

"""
    cf_effective_field(B, n_e) → B* [T]

Effective magnetic field seen by composite fermions: B* = B - 2n_e φ₀.
"""
function cf_effective_field(B::Real, n_e::Real)
    φ₀ = CONSTANTS.h / CONSTANTS.e  # flux quantum
    return B - 2 * n_e * φ₀
end

"""
    cf_gap_scaling(Δ_anchor, p; r=0.56) → Δ [e²/(ε ℓ_B)]

Gap for Jain fraction ν = p/(2p+1) via geometric CF scaling:
    Δ(p) = Δ_anchor × r^(p-1)

where r ≈ 0.56 is calibrated from Δ(2/5)/Δ(1/3) in published ED data
(Fano 1986, Balram 2018).  Predictions vs published:
    p=2 (2/5): r   = 0.56 → 0.058/0.1036 = 0.56  (calibrated)
    p=3 (3/7): r²  = 0.31 → 0.030/0.1036 = 0.29  (within 7%)
    p=4 (4/9): r³  = 0.18 → 0.017/0.1036 = 0.16  (within 10%)

Anchor: Δ_anchor = charge gap for ν=1/3 (p=1) from ED.
"""
function cf_gap_scaling(Δ_anchor::Real, p::Int; r::Real=0.56)
    p >= 1 || error("p must be ≥ 1, got $p")
    return Δ_anchor * r^(p - 1)
end
