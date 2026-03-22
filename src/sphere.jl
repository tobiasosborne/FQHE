# sphere.jl — Haldane sphere geometry
#
# Sources:
#   P1: Haldane, PRL 51, 605 (1983)
#   P8: Jain, Composite Fermions, Ch. 4

export sphere_flux, sphere_radius, shift

"""
    sphere_flux(N, ν) → 2S

Monopole strength (twice the flux quanta) for N electrons at filling ν on the sphere.
2S = N/ν - S_shift, where S_shift is the topological shift.
"""
function sphere_flux(N::Integer, ν::Rational)
    S_shift = shift(ν)
    twoS = N ÷ numerator(ν) * denominator(ν) - S_shift
    return twoS
end

"""
    shift(ν) → S_shift

Topological shift for filling fraction ν.
- Laughlin 1/q: S_shift = q
- Jain p/(2p+1): S_shift = 2p+1  (Jain's convention)
- Jain p/(2p-1): S_shift = 2p-1
"""
function shift(ν::Rational)
    p, q = numerator(ν), denominator(ν)
    if p == 1
        return q  # Laughlin state
    elseif q == 2p + 1
        return q  # principal Jain sequence
    elseif q == 2p - 1
        return q  # second Jain sequence
    else
        error("Shift not implemented for ν = $ν")
    end
end

"""
    sphere_radius(S) → R/ℓ_B

Radius of Haldane sphere in units of magnetic length: R = √S ℓ_B.
"""
sphere_radius(S::Real) = sqrt(S)
