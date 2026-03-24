# hamiltonian.jl — Sparse Hamiltonian assembly from pseudopotentials
#
# Two-body Coulomb Hamiltonian in the LLL on Haldane sphere.
# H = (1/2) Σ_{abcd} ⟨ab|V|cd⟩ a†_a a†_b a_d a_c
# Fano et al. (1986) Eq. (23)–(24).
#
# Sources:
#   P01: Haldane, PRL 51, 605 (1983)
#   P04: Fano et al., PRB 34, 2670 (1986)

using SparseArrays, LinearAlgebra

export build_hamiltonian

"""
    build_hamiltonian(basis, VJ) → Hermitian sparse matrix

Build the Coulomb Hamiltonian in the Fock basis.
VJ: pseudopotentials from `coulomb_pseudopotentials(twoS)`.
Returns H in units of e²/(ε ℓ_B).
"""
function build_hamiltonian(basis::FockBasis, VJ::Vector{Float64})
    twoS = basis.twoS
    n_orb = basis.n_orb
    D = length(basis.states)

    # --- Step 1: Precompute two-body matrix elements ---
    Vmat = _precompute_matrix_elements(twoS, VJ)

    # --- Step 2: Assemble sparse Hamiltonian ---
    I_idx = Int[];  J_idx = Int[];  V_val = Float64[]
    sizehint!(I_idx, D * 20)   # rough estimate of nonzeros
    sizehint!(J_idx, D * 20)
    sizehint!(V_val, D * 20)

    for (β_idx, β) in enumerate(basis.states)
        occ = occupied_orbitals(β, n_orb)
        N_e = length(occ)

        # Diagonal: Σ_{a<b ∈ occ} 2 V_{ab,ab}
        diag = 0.0
        for ii in 1:N_e-1, jj in ii+1:N_e
            diag += 2.0 * get(Vmat, (occ[ii], occ[jj], occ[ii], occ[jj]), 0.0)
        end
        if abs(diag) > 1e-15
            push!(I_idx, β_idx); push!(J_idx, β_idx); push!(V_val, diag)
        end

        # Off-diagonal: annihilate (c,d), create (a,b)
        for ii in 1:N_e-1, jj in ii+1:N_e
            c, d = occ[ii], occ[jj]
            β_strip = β & ~(UInt64(1) << c) & ~(UInt64(1) << d)
            twoM_cd = (2c - twoS) + (2d - twoS)

            for a in 0:n_orb-1
                (β_strip >> a) & 1 == 0 || continue
                for b in a+1:n_orb-1
                    (β_strip >> b) & 1 == 0 || continue
                    (2a - twoS) + (2b - twoS) == twoM_cd || continue
                    (a == c && b == d) && continue   # diagonal handled above

                    v = get(Vmat, (a, b, c, d), 0.0)
                    abs(v) > 1e-15 || continue

                    α = β_strip | (UInt64(1) << a) | (UInt64(1) << b)
                    α_idx = state_index(basis, α)
                    α_idx > 0 || continue

                    sgn = _fermion_sign(β, c, d, a, b)
                    push!(I_idx, α_idx); push!(J_idx, β_idx); push!(V_val, 2.0 * v * sgn)
                end
            end
        end
    end

    H = sparse(I_idx, J_idx, V_val, D, D)
    return Hermitian(H, :U)   # uses upper triangle → exact symmetry
end

# ---------- internal helpers ----------

"""Precompute V_{ab,cd} for all ordered pairs (a<b), (c<d) with same 2M."""
function _precompute_matrix_elements(twoS::Integer, VJ::Vector{Float64})
    n_orb = twoS + 1
    S = twoS / 2.0

    # Group orbital pairs by 2M
    pairs_by_twoM = Dict{Int, Vector{Tuple{Int,Int}}}()
    for a in 0:n_orb-1, b in a+1:n_orb-1
        twoM = (2a - twoS) + (2b - twoS)
        push!(get!(pairs_by_twoM, twoM, Tuple{Int,Int}[]), (a, b))
    end

    Vmat = Dict{NTuple{4,Int}, Float64}()
    inv_sqrtS = 1.0 / sqrt(S)

    for (twoM, pairs) in pairs_by_twoM
        M = twoM / 2.0
        Jmin = ceil(Int, abs(M))

        # Precompute CG coefficients for this M
        n_pairs = length(pairs)
        # CG_table[pair_index, J_index]
        fermionic_Js = [J for J in Jmin:twoS if isodd(twoS - J)]
        isempty(fermionic_Js) && continue

        cg_table = zeros(Float64, n_pairs, length(fermionic_Js))
        for (pi, (a, b)) in enumerate(pairs)
            m_a = (2a - twoS) / 2.0
            m_b = (2b - twoS) / 2.0
            for (ji, J) in enumerate(fermionic_Js)
                cg_table[pi, ji] = clebschgordan(S, m_a, S, m_b, Float64(J), M)
            end
        end

        # V_{ab,cd} = Σ_J VJ[J+1] * cg(ab,J) * cg(cd,J) / √S
        for pi in 1:n_pairs, qi in 1:n_pairs
            val = 0.0
            for (ji, J) in enumerate(fermionic_Js)
                val += VJ[J + 1] * cg_table[pi, ji] * cg_table[qi, ji]
            end
            val *= inv_sqrtS
            if abs(val) > 1e-15
                Vmat[(pairs[pi]..., pairs[qi]...)] = val
            end
        end
    end
    return Vmat
end

"""
Fermion sign of a†_a a†_b a_d a_c |state⟩ where c < d are annihilated, a < b created.
Operator order: first a_c, then a_d, then a†_b, then a†_a.
"""
function _fermion_sign(state::UInt64, c::Int, d::Int, a::Int, b::Int)
    s = state
    sign = 1

    # Annihilate c (rightmost operator acts first)
    mask = (UInt64(1) << c) - 1
    isodd(count_ones(s & mask)) && (sign = -sign)
    s &= ~(UInt64(1) << c)

    # Annihilate d
    mask = (UInt64(1) << d) - 1
    isodd(count_ones(s & mask)) && (sign = -sign)
    s &= ~(UInt64(1) << d)

    # Create b
    mask = (UInt64(1) << b) - 1
    isodd(count_ones(s & mask)) && (sign = -sign)
    s |= UInt64(1) << b

    # Create a
    mask = (UInt64(1) << a) - 1
    isodd(count_ones(s & mask)) && (sign = -sign)

    return sign
end
