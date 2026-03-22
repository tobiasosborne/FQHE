# hilbert_space.jl — Many-body Fock space construction
#
# Bit-string representation of Fock states in the LLL.
# Sources:
#   P1: Haldane, PRL 51, 605 (1983)
#   P9: Chakraborty & Pietiläinen, The Quantum Hall Effects, Ch. 5

using Combinatorics

export FockBasis, enumerate_fock_states, lz_sector

"""
Fock space basis for N fermions in n_orb orbitals.
States stored as UInt64 bitstrings (sufficient for n_orb ≤ 64).
"""
struct FockBasis
    N::Int           # number of electrons
    n_orb::Int       # number of orbitals (= 2S+1)
    twoS::Int        # twice the monopole strength
    states::Vector{UInt64}   # bit-string Fock states in this sector
    total_Lz::Int    # total Lz of this sector (in units of 1, for integer S)
end

"""
    enumerate_fock_states(N, twoS; Lz=0) → FockBasis

Enumerate all N-fermion Fock states in 2S+1 LLL orbitals with total Lz = Lz.
Orbitals labeled m = -S, -S+1, ..., S (integer spacing for integer 2S).
"""
function enumerate_fock_states(N::Integer, twoS::Integer; Lz::Integer=0)
    n_orb = twoS + 1
    S = twoS ÷ 2  # works for even twoS; need half-integer handling otherwise

    states = UInt64[]

    # Enumerate all C(n_orb, N) configurations
    for combo in combinations(0:n_orb-1, N)
        # combo contains orbital indices 0, 1, ..., n_orb-1
        # map to m-values: m_i = index_i - S
        total_lz = sum(idx - S for idx in combo)
        if total_lz == Lz
            # encode as bitstring
            bits = UInt64(0)
            for idx in combo
                bits |= UInt64(1) << idx
            end
            push!(states, bits)
        end
    end

    return FockBasis(N, n_orb, twoS, states, Lz)
end

"""
    occupied_orbitals(state, n_orb) → Vector{Int}

Return sorted list of occupied orbital indices from a bitstring state.
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
