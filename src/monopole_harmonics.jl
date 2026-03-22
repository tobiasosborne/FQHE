# monopole_harmonics.jl — Single-particle basis on Haldane sphere
#
# For LLL exact diagonalisation we don't need explicit harmonics —
# the two-body matrix elements reduce to Haldane pseudopotentials.
# This file provides angular momentum algebra utilities.
#
# Sources:
#   P1: Haldane, PRL 51, 605 (1983)
#   Wu & Yang, Nucl. Phys. B 107, 365 (1976)

export lll_orbitals

"""
    lll_orbitals(twoS) → Vector{Int}

Return the m-values of LLL orbitals on a sphere with monopole strength S = twoS/2.
m ranges from -S to S in integer steps, giving 2S+1 orbitals.
"""
function lll_orbitals(twoS::Integer)
    S = twoS
    return collect(-S:2:S)  # m values in units of 1/2, stored as integers 2m
end
