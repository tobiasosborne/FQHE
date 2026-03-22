# integer_qhe.jl — Integer QHE: Landau level gaps (analytic)
#
# No ED needed — the gap is simply ℏω_c.
# Sources:
#   Standard, e.g. Girvin lectures (P11) §2

export iqhe_gap, iqhe_Rxy

"""
    iqhe_gap(n, B, mat) → Δ [J]

Gap for integer filling ν = n: Δ = ℏω_c = ℏeB/m*.
"""
iqhe_gap(n::Integer, B::Real, mat::MaterialParams=GaAs) = cyclotron_energy(B, mat)

"""
    iqhe_Rxy(n) → R_xy [Ω]

Hall resistance at integer filling ν = n: R_xy = h/(ne²).
"""
iqhe_Rxy(n::Integer) = CONSTANTS.h / (n * CONSTANTS.e^2)
