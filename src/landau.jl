# landau.jl — Landau levels, filling factor, energy scales
#
# Sources:
#   Standard QHE theory, e.g. Girvin lectures (P11) §2

export filling_factor, B_at_filling, landau_level_energy

"""
    filling_factor(B, n_e) → ν

Filling factor ν = n_e h / (eB) = n_e × 2π ℓ_B².
n_e in m⁻² (SI), B in Tesla.
"""
filling_factor(B::Real, n_e::Real) = n_e * CONSTANTS.h / (CONSTANTS.e * B)

"""
    B_at_filling(ν, n_e) → B [T]

Magnetic field at which filling factor equals ν.
"""
B_at_filling(ν::Real, n_e::Real) = n_e * CONSTANTS.h / (CONSTANTS.e * ν)

"""
    landau_level_energy(n, B, mat) → E_n [J]

Energy of n-th Landau level: E_n = (n + 1/2) ħω_c.
"""
landau_level_energy(n::Integer, B::Real, mat::MaterialParams=GaAs) =
    (n + 0.5) * cyclotron_energy(B, mat)
