# dmrg_hamiltonian.jl — Build FQHE cylinder Hamiltonian as ITensors MPO
#
# Constructs the LLL-projected Coulomb Hamiltonian on a cylinder
# using ITensors OpSum (fermionic sites with QN conservation).
#
# Sources:
#   P25: Zaletel et al., PRL 110, 236801 (2013) — iDMRG for FQHE
#   P28: Rezayi & Haldane, PRL 84, 4685 (2000) — cylinder geometry

using ITensors, ITensorMPS

export build_cylinder_mpo

"""
    _precompute_cylinder_elements(Nphi, Ly; mmax=15, range_cutoff=0) → Dict

Precompute all nonzero cylinder matrix elements V_{abcd} for a<b, c<d, a+b=c+d.
Returns a Dict mapping (a,b,c,d) → V.

`range_cutoff`: max orbital range for interaction terms. If 0, set to min(Nphi-1, ceil(2*Ly/π)).
The Coulomb interaction decays as exp(-r²Δ²/4), so terms beyond ~2Ly/π orbital spacings
contribute < 1e-4 and can be safely dropped. This dramatically speeds up MPO construction.
"""
function _precompute_cylinder_elements(Nphi::Int, Ly::Float64; mmax::Int=15, range_cutoff::Int=0)
    Vm = plane_pseudopotentials(mmax)
    Δ = 2π / Ly

    # Auto range cutoff: pair amplitude g_m(r) ~ exp(-r²Δ²/4), negligible for r > Rmax
    if range_cutoff <= 0
        range_cutoff = min(Nphi - 1, ceil(Int, 2.5 * Ly / π))
    end

    # Precompute pair amplitudes g_m(r) for all needed r and m
    rmax = min(Nphi - 1, range_cutoff)
    g_table = zeros(mmax + 1, 2rmax + 1)  # g_table[m+1, r+rmax+1]
    for m in 1:2:mmax
        for r in -rmax:rmax
            g_table[m+1, r+rmax+1] = _pair_amplitude(m, r, Float64(Ly))
        end
    end

    elements = Dict{NTuple{4,Int}, Float64}()

    for a in 0:Nphi-1, b in a+1:Nphi-1
        abs(a - b) <= range_cutoff || continue  # bra pair range cutoff
        K = a + b
        r_bra = a - b

        for c in max(0, K - Nphi + 1):min(Nphi - 1, K)
            d = K - c
            d > c || continue
            0 <= d < Nphi || continue
            abs(c - d) <= range_cutoff || continue  # ket pair range cutoff

            r_ket = c - d

            val = 0.0
            for m in 1:2:mmax
                val += Vm[m+1] * g_table[m+1, r_bra+rmax+1] * g_table[m+1, r_ket+rmax+1]
            end
            abs(val) > 1e-12 || continue
            elements[(a, b, c, d)] = val
        end
    end
    return elements
end

"""
    build_cylinder_mpo(Nphi, Ne, Ly; mmax=30) → (sites, H, n_terms)

Build the LLL Coulomb Hamiltonian on a cylinder as an ITensors MPO.
Uses pseudopotential decomposition with precomputed pair amplitudes.
"""
function build_cylinder_mpo(Nphi::Int, Ne::Int, Ly::Float64; mmax::Int=30)
    sites = siteinds("Fermion", Nphi; conserve_qns=true)

    # Precompute all matrix elements
    elements = _precompute_cylinder_elements(Nphi, Ly; mmax)

    ampo = OpSum()
    n_terms = 0

    for ((a, b, c, d), V) in elements
        if a == c && b == d
            # Diagonal: n_a n_b
            ampo += V, "N", a + 1, "N", b + 1
        else
            # Off-diagonal: c†_a c†_b c_d c_c
            ampo += V, "Cdag", a + 1, "Cdag", b + 1, "C", d + 1, "C", c + 1
        end
        n_terms += 1
    end

    H = MPO(ampo, sites)
    return sites, H, n_terms
end
