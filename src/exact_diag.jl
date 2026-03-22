# exact_diag.jl — Lanczos/Arnoldi for ground state and gap extraction
#
# Sources:
#   KrylovKit.jl documentation
#   P9: Chakraborty & Pietiläinen, §5.4

using KrylovKit

export compute_spectrum, compute_gap

# TODO: Implement after hamiltonian.jl is validated
function compute_spectrum(H; nev=5)
    error("Not yet implemented — awaiting Hamiltonian assembly (Step 5)")
end

function compute_gap(H)
    error("Not yet implemented")
end
