# materials.jl — Physical constants and GaAs/AlGaAs parameters
#
# Sources:
#   D3: m*=0.067 m_e, ε=12.9, g*=-0.44
#   Ando, Fowler & Stern, Rev. Mod. Phys. 54 (1982)
#   CODATA 2018 for fundamental constants

export PhysicalConstants, CONSTANTS, GaAs, MaterialParams, magnetic_length, cyclotron_energy, coulomb_energy

"""
Fundamental physical constants (SI units, CODATA 2018).
"""
struct PhysicalConstants
    ħ::Float64    # reduced Planck constant [J·s]
    e::Float64    # elementary charge [C]
    m_e::Float64  # electron mass [kg]
    k_B::Float64  # Boltzmann constant [J/K]
    ε_0::Float64  # vacuum permittivity [F/m]
    h::Float64    # Planck constant [J·s]
end

const CONSTANTS = PhysicalConstants(
    1.054571817e-34,   # ħ
    1.602176634e-19,   # e
    9.1093837015e-31,  # m_e
    1.380649e-23,      # k_B
    8.8541878128e-12,  # ε_0
    6.62607015e-34,    # h
)

"""
Material parameters for a 2DEG heterostructure.
"""
struct MaterialParams
    m_star::Float64    # effective mass [kg]
    ε_r::Float64       # relative dielectric constant
    g_star::Float64    # effective g-factor
    name::String
end

"""GaAs/AlGaAs heterostructure — the standard FQHE system."""
const GaAs = MaterialParams(
    0.067 * CONSTANTS.m_e,   # m* = 0.067 m_e
    12.9,                     # ε_r = 12.9
    -0.44,                    # g* = -0.44
    "GaAs/AlGaAs",
)

"""
    magnetic_length(B) → ℓ_B [m]

Magnetic length ℓ_B = √(ħ/eB).
At B = 1 T: ℓ_B ≈ 25.66 nm.
"""
magnetic_length(B::Real) = sqrt(CONSTANTS.ħ / (CONSTANTS.e * B))

"""
    cyclotron_energy(B, mat::MaterialParams) → ħω_c [J]

Cyclotron energy ħω_c = ħeB/m*.
At B = 1 T for GaAs: ħω_c/k_B ≈ 20 K.
"""
cyclotron_energy(B::Real, mat::MaterialParams=GaAs) = CONSTANTS.ħ * CONSTANTS.e * B / mat.m_star

"""
    coulomb_energy(B, mat::MaterialParams) → E_C [J]

Coulomb energy scale e²/(4πε₀ε_r ℓ_B).
This is the natural energy unit for FQHE gaps.
"""
function coulomb_energy(B::Real, mat::MaterialParams=GaAs)
    ℓ_B = magnetic_length(B)
    return CONSTANTS.e^2 / (4π * CONSTANTS.ε_0 * mat.ε_r * ℓ_B)
end
