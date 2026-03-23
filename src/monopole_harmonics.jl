# monopole_harmonics.jl — Single-particle basis on Haldane sphere
#
# For LLL exact diagonalisation we don't need explicit harmonics —
# the two-body matrix elements reduce to Haldane pseudopotentials.
# This file provides angular momentum algebra utilities.
#
# Sources:
#   P01: Haldane, PRL 51, 605 (1983)
#   P04: Fano et al., PRB 34, 2670 (1986), Eq. (21)

export lll_orbitals, orbital_m

"""
    lll_orbitals(twoS) → 0:twoS

Orbital indices for LLL on sphere with monopole strength S = twoS/2.
There are twoS + 1 orbitals.  Use `orbital_m(i, twoS)` for the m-value.
"""
lll_orbitals(twoS::Integer) = 0:twoS

"""
    orbital_m(i, twoS) → Float64

Angular momentum m-value of orbital i: m = (2i − twoS)/2.
"""
orbital_m(i::Integer, twoS::Integer) = (2i - twoS) / 2.0
