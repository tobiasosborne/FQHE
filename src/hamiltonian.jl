# hamiltonian.jl — Sparse Hamiltonian assembly from pseudopotentials
#
# Two-body Hamiltonian in the LLL on Haldane sphere.
# Sources:
#   P1: Haldane, PRL 51, 605 (1983)
#   P8: Jain, Composite Fermions, Ch. 4

using SparseArrays, LinearAlgebra, WignerSymbols

export build_hamiltonian

# TODO: Implement after pseudopotentials and hilbert_space are validated
function build_hamiltonian(basis::FockBasis, Vm::Vector{Float64})
    error("Not yet implemented — awaiting validation of pseudopotentials (Step 3)")
end
