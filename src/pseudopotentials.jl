# pseudopotentials.jl — Haldane pseudopotentials for Coulomb interaction in LLL
#
# Closed-form formula from Fano, Ortolani & Colombo, PRB 34, 2670 (1986).
# Eq. (25) = Eq. (A15):
#   V_J^(S) = 2 * C(4S-2J, 2S-J) * C(4S+2J+2, 2S+J+1) / C(4S+2, 2S+1)^2
#
# V_J is the pseudopotential for pair angular momentum J on a sphere
# with monopole strength S = twoS/2.  Units: e²/(ε R) where R = √S ℓ_B.
# For fermions, only J with odd relative angular momentum (twoS - J odd) contribute.
#
# Sources:
#   P04: Fano et al., PRB 34, 2670 (1986), Eq. (25), Appendix
#   P01: Haldane, PRL 51, 605 (1983)

export coulomb_pseudopotentials, two_body_element

"""
    coulomb_pseudopotentials(twoS) → Vector{Float64}

Haldane pseudopotentials V_J for Coulomb 1/r on the sphere (Fano eq. 25).
Returns VJ[J+1] for J = 0, 1, ..., twoS.  Non-fermionic channels are zero.
Units: dimensionless (multiply by 1/R = 1/√S for energy in e²/(ε ℓ_B)).
"""
function coulomb_pseudopotentials(twoS::Integer)
    twoS >= 1 || error("Need twoS ≥ 1, got $twoS")

    # Use BigInt to avoid overflow in binomial coefficients
    t = BigInt(twoS)
    denom = binomial(2t + 2, t + 1)^2

    VJ = zeros(Float64, twoS + 1)
    for J in 0:twoS
        isodd(twoS - J) || continue          # fermions only: odd relative ang. mom.
        num = 2 * binomial(2t - 2J, t - J) * binomial(2t + 2J + 2, t + J + 1)
        VJ[J + 1] = Float64(num) / Float64(denom)
    end
    return VJ
end

"""
    two_body_element(twoS, i_a, i_b, i_c, i_d, VJ) → Float64

Two-body matrix element ⟨a,b|V|c,d⟩ on the sphere (Fano eq. 24).
Orbital indices i ∈ 0:twoS, with m_i = (2i - twoS)/2.
Includes the 1/R = 1/√(twoS/2) factor.  Units: e²/(ε ℓ_B).
"""
function two_body_element(twoS::Integer, i_a::Int, i_b::Int,
                          i_c::Int, i_d::Int, VJ::Vector{Float64})
    # m-values (may be half-integer)
    S = twoS / 2.0
    m_a = (2i_a - twoS) / 2.0
    m_b = (2i_b - twoS) / 2.0
    m_c = (2i_c - twoS) / 2.0
    m_d = (2i_d - twoS) / 2.0

    M = m_a + m_b
    abs(M - (m_c + m_d)) < 1e-10 || return 0.0

    val = 0.0
    Jmin = ceil(Int, abs(M))
    for J in Jmin:twoS
        isodd(twoS - J) || continue
        cg1 = clebschgordan(S, m_a, S, m_b, Float64(J), M)
        cg2 = clebschgordan(S, m_c, S, m_d, Float64(J), M)
        val += VJ[J + 1] * cg1 * cg2
    end

    return val / sqrt(S)   # 1/R factor: energy in e²/(ε ℓ_B)
end
