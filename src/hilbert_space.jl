# hilbert_space.jl — Many-body Fock space construction
#
# Bit-string representation of Fock states in the LLL on Haldane sphere.
# Orbital i ∈ 0:twoS has angular momentum m_i = (2i − twoS)/2.
# We work with 2×Lz (integer) to avoid half-integer issues.
#
# Sources:
#   P01: Haldane, PRL 51, 605 (1983)
#   P04: Fano et al., PRB 34, 2670 (1986), Sec. IV

using Combinatorics

export FockBasis, enumerate_fock_states, occupied_orbitals, state_index

"""
Fock space basis for N fermions in n_orb = twoS+1 LLL orbitals.
States stored as UInt64 bitstrings (bit i = orbital i occupied).
"""
struct FockBasis
    N::Int                          # number of electrons
    n_orb::Int                      # number of orbitals = twoS + 1
    twoS::Int                       # twice the monopole strength
    states::Vector{UInt64}          # Fock states in this Lz sector
    target_twoLz::Int               # twice the total Lz
    state_to_index::Dict{UInt64,Int}  # fast reverse lookup
end

"""
    enumerate_fock_states(N, twoS; twoLz=0) → FockBasis

Enumerate all N-fermion states in twoS+1 orbitals with total 2Lz = twoLz.
"""
function enumerate_fock_states(N::Integer, twoS::Integer; twoLz::Integer=0)
    n_orb = twoS + 1
    n_orb <= 64 || error("n_orb=$n_orb > 64: need UInt128 (not implemented)")
    N <= n_orb || error("N=$N > n_orb=$n_orb: impossible")

    states = UInt64[]
    for combo in combinations(0:n_orb-1, N)
        total_twoLz = sum(2*idx - twoS for idx in combo)
        if total_twoLz == twoLz
            bits = UInt64(0)
            for idx in combo
                bits |= UInt64(1) << idx
            end
            push!(states, bits)
        end
    end

    lookup = Dict{UInt64,Int}(s => i for (i, s) in enumerate(states))
    return FockBasis(N, n_orb, twoS, states, twoLz, lookup)
end

"""
    occupied_orbitals(state, n_orb) → Vector{Int}

Sorted list of occupied orbital indices from a bitstring state.
"""
function occupied_orbitals(state::UInt64, n_orb::Integer)
    orbs = Int[]
    for i in 0:n_orb-1
        if (state >> i) & 1 == 1
            push!(orbs, i)
        end
    end
    return orbs
end

"""
    state_index(basis, state) → Int or 0

Look up basis index for a bitstring state.  Returns 0 if not found.
"""
state_index(basis::FockBasis, state::UInt64) = get(basis.state_to_index, state, 0)
