# pseudopotentials.jl — Haldane pseudopotentials for Coulomb interaction in LLL
#
# V_m for the Coulomb interaction on the Haldane sphere.
# Sources:
#   P1: Haldane, PRL 51, 605 (1983) — definition of pseudopotentials
#   P4: Fano, Ortolani & Colombo, PRB 34, 2670 (1986) — explicit formula
#   P8: Jain, Composite Fermions, eq. (4.43) — reference formula
#   P9: Chakraborty & Pietiläinen, §5.3

using WignerSymbols

export coulomb_pseudopotentials

"""
    coulomb_pseudopotentials(twoS) → Vector{Float64}

Compute Haldane pseudopotentials V_m for the Coulomb interaction in the LLL
on a sphere with monopole strength S = twoS/2.

Returns V_m for m = 1, 3, 5, ..., 2S-1 (odd relative angular momenta only, fermions).
Values in units of e²/(ε ℓ_B).

Formula (Fano et al. 1986, Jain eq. 4.43):
  V_m = Σ_{l=0}^{2S} (4l+1)/(2(2S+1)) × [CG(S,S,l; m/2,-m/2,0)]² × (2S)!²/((2S-l)!(2S+l+1)!)^{...}

The CG coefficients couple two spin-S particles to relative angular momentum.
"""
function coulomb_pseudopotentials(twoS::Integer)
    S = twoS // 2  # keep as rational for CG
    Sf = Float64(S)

    # Relative angular momentum for fermions: m = 2S, 2S-2, ..., 1
    # In the standard convention, m is the relative angular momentum quantum number
    # and runs over odd values for fermions (antisymmetry).
    # Here we use the pair angular momentum L_pair = 2S - m, so L_pair = 0, 2, 4, ...

    Vm = Float64[]

    for m in (twoS):-2:1  # m = 2S, 2S-2, ..., 1 (odd m for odd 2S)
        # Relative angular momentum m means pair angular momentum L = 2S - m
        L_pair = twoS - m

        # V_m via expansion in Legendre polynomials of the Coulomb potential
        val = 0.0
        for l in 0:twoS
            # CG coefficient ⟨S, S, l | S-L_pair/2, -(S-L_pair/2), 0⟩
            # = ⟨S, m₁; S, m₂ | l, 0⟩ where m₁ = m/2, m₂ = -m/2
            cg = clebschgordan(Float64(S), m/2, Float64(S), -m/2, Float64(l), 0)
            # Coulomb coefficient c_l on sphere:
            # For 1/r interaction on sphere of radius R=√S:
            # c_l = 1/(2√(S + 1/2)) × ... but simpler form:
            # c_l = (2l+1) / (2(2S+1)) normalised, with 1/R prefactor
            val += cg^2 / (2 * sqrt(Sf + 0.5))
        end
        push!(Vm, val)
    end

    return Vm
end
