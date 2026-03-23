# sphere.jl — Haldane sphere geometry
#
# Sources:
#   P01: Haldane, PRL 51, 605 (1983)
#   P04: Fano et al., PRB 34, 2670 (1986), Sec. II

export sphere_flux, sphere_radius, shift

"""
    sphere_flux(N, ν) → 2S::Int

Monopole strength 2S for N electrons at filling ν on the Haldane sphere.
2S = N/ν − shift(ν).  Uses exact rational arithmetic; errors if non-integer.
"""
function sphere_flux(N::Integer, ν::Rational)
    twoS_rat = N // ν - shift(ν)
    isinteger(twoS_rat) || error("N=$N incompatible with ν=$ν: 2S = $twoS_rat is not integer")
    twoS = Int(twoS_rat)
    twoS >= 0 || error("N=$N too small for ν=$ν (2S=$twoS < 0)")
    return twoS
end

"""
    shift(ν) → S_shift::Int

Topological shift for filling fraction ν on the Haldane sphere.
For all Jain fractions ν = p/(2p±1) and their particle-hole conjugates,
and for Laughlin states ν = 1/m, the shift equals the denominator.

Ref: Jain, Composite Fermions, Table 4.2.
"""
function shift(ν::Rational)
    p, q = numerator(ν), denominator(ν)
    0 < p < q || error("Need 0 < p < q for ν = $p/$q")

    # Laughlin: ν = 1/q (q odd)
    p == 1 && isodd(q) && return q

    # Principal Jain sequence: ν = p/(2p+1)
    q == 2p + 1 && return q

    # Second Jain sequence: ν = p/(2p-1)  [includes particle-hole conjugates]
    q == 2p - 1 && return q

    # Particle-hole conjugate: ν = (q-p)/q of a known fraction
    ν_conj = (q - p) // q
    return shift(ν_conj)
end

"""
    sphere_radius(twoS) → R/ℓ_B

Radius of Haldane sphere in units of magnetic length: R = √(twoS/2).
"""
sphere_radius(twoS::Integer) = sqrt(twoS / 2.0)
