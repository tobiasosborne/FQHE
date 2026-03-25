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
The shift equals the denominator q for ALL odd-denominator FQH fractions
(Laughlin, Jain principal/reverse at any flux attachment, and their
particle-hole conjugates). This is the universal rule from Wen's
topological classification.

Special cases:
- ν = 1/2 (CFL): gapless composite Fermi liquid, no shift defined.
- Even-denominator fractions (e.g., 5/2): handled via effective LL mapping.

Ref: Jain, Composite Fermions, Table 4.2; Wen, Adv. Phys. 44, 405 (1995).
"""
function shift(ν::Rational)
    p, q = numerator(ν), denominator(ν)
    0 < p < q || error("Need 0 < p < q for ν = $p/$q")

    # ν = 1/2: Composite Fermi Liquid — gapless, no shift
    p == 1 && q == 2 && error("ν = 1/2 is a gapless CFL state; shift is undefined")

    # Universal rule: shift = denominator for all odd-denominator FQH fractions
    isodd(q) && return q

    # Even denominators (e.g., 5/2 mapped to effective ν*) need special handling
    error("ν = $ν has even denominator; use effective LL mapping")
end

"""
    sphere_radius(twoS) → R/ℓ_B

Radius of Haldane sphere in units of magnetic length: R = √(twoS/2).
"""
sphere_radius(twoS::Integer) = sqrt(twoS / 2.0)
