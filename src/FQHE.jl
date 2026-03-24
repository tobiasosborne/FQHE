module FQHE

# Physical constants and GaAs material parameters
include("materials.jl")

# Landau levels, filling factor, energy scales
include("landau.jl")

# Haldane sphere geometry
include("sphere.jl")

# Monopole harmonics (single-particle basis on sphere)
include("monopole_harmonics.jl")

# Many-body Fock space construction
include("hilbert_space.jl")

# Numerical Clebsch-Gordan coefficients (replaces WignerSymbols.jl)
include("clebsch_gordan.jl")

# Haldane pseudopotentials for Coulomb interaction
include("pseudopotentials.jl")

# Sparse Hamiltonian assembly
include("hamiltonian.jl")

# Lanczos/Arnoldi exact diagonalisation
include("exact_diag.jl")

# Composite fermion theory (cross-check)
include("composite_fermion.jl")

# Laughlin wavefunction overlaps (validation)
include("laughlin.jl")

# Integer QHE (analytic)
include("integer_qhe.jl")

# Transport model: gaps → R_xx, R_xy
include("transport.jl")

end # module
