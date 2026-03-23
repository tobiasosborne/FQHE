# transport.jl — Gaps → R_xx, R_xy vs B
#
# Semi-phenomenological transport model.
# Ab initio: plateau positions (h/νe²), gap magnitudes, hierarchy.
# Phenomenological: plateau widths (via mobility), transition shapes (semicircle law).
#
# Sources:
#   P10: Dykhne & Ruzin, PRB 50, 2369 (1994) — semicircle law
#   P11: Girvin lectures, §2-3

export compute_transport, FractionData

"""Data for one QH fraction used by the transport model."""
struct FractionData
    ν::Rational{Int}
    gap_dimless::Float64     # gap in units of e²/(ε ℓ_B)
    is_integer::Bool
end

"""
    compute_transport(B_range, n_e, fractions, T; mat=GaAs, μ_e=1e2)
    → (R_xy, R_xx)   [arrays same length as B_range]

Compute R_xy and R_xx vs B using gaps + semicircle transport law.
μ_e in m²/(V·s)  (1e6 cm²/Vs = 1e2 m²/Vs).
"""
function compute_transport(B_range::AbstractVector, n_e::Real,
                           fractions::Vector{FractionData}, T::Real;
                           mat::MaterialParams=GaAs, μ_e::Real=1e2)
    h = CONSTANTS.h
    e = CONSTANTS.e
    kB = CONSTANTS.k_B
    RK = h / e^2                    # von Klitzing constant ≈ 25812.807 Ω

    # Sort fractions by ν (descending) → increasing B
    sorted = sort(fractions, by=f -> Float64(f.ν), rev=true)
    ν_vals = [Float64(f.ν) for f in sorted]
    B_centers = [n_e * h / (e * ν) for ν in ν_vals]
    gap_phys = [f.gap_dimless * coulomb_energy(B, mat) / kB   # gap in Kelvin
                for (f, B) in zip(sorted, B_centers)]

    # Plateau half-widths in ν-space: proportional to gap/disorder
    Γ = CONSTANTS.ħ * e / (mat.m_star * μ_e)   # disorder broadening [J]
    w_ν = [min(Float64(f.ν) * 0.4,    # cap at 40% of ν
               f.gap_dimless > 0 ? 0.05 * f.gap_dimless * coulomb_energy(B, mat) / Γ : 0.0)
           for (f, B) in zip(sorted, B_centers)]

    R_xy = similar(B_range, Float64)
    R_xx = similar(B_range, Float64)

    for (bi, B) in enumerate(B_range)
        B <= 0 && (R_xy[bi] = 0.0; R_xx[bi] = 0.0; continue)

        ν_B = n_e * h / (e * B)

        # Find which plateau we're on (if any)
        on_plateau = false
        for (fi, f) in enumerate(sorted)
            νf = ν_vals[fi]
            if abs(ν_B - νf) < w_ν[fi] && w_ν[fi] > 0
                on_plateau = true
                R_xy[bi] = RK / νf
                # Activated R_xx at plateau center
                Δ_K = gap_phys[fi]
                R_xx[bi] = Δ_K > 0 ? 500.0 * exp(-Δ_K / (2 * T)) : 0.0
                break
            end
        end
        on_plateau && continue

        # Between plateaus: find bounding fractions
        ν_above = nothing; ν_below = nothing
        for fi in 1:length(ν_vals)
            if ν_vals[fi] > ν_B
                ν_above = ν_vals[fi]
            end
            if ν_vals[fi] < ν_B && ν_below === nothing
                ν_below = ν_vals[fi]
            end
        end

        # Classical Hall line as fallback
        if ν_above === nothing || ν_below === nothing
            R_xy[bi] = B / (n_e * e)
            # SdH oscillations for high filling
            if ν_B > 1
                τ = mat.m_star * μ_e / e
                ωc = e * B / mat.m_star
                dingle = exp(-π / (ωc * τ))
                thermal = ν_B > 0.5 ? 2π^2 * kB * T / (CONSTANTS.ħ * ωc) : 1e10
                amp = thermal < 50 ? dingle * thermal / sinh(thermal) : 0.0
                R_xx[bi] = RK * (1.0 - 4.0 * amp * cos(2π * ν_B)) / ν_B^2 * 0.02
            else
                R_xx[bi] = 0.0
            end
            continue
        end

        # Semicircle interpolation (Dykhne-Ruzin)
        # Parameter t ∈ [0,1]: 0 at ν_above plateau edge, 1 at ν_below edge
        ν_mid = (ν_above + ν_below) / 2
        ν_span = ν_above - ν_below
        t = clamp((ν_above - ν_B) / ν_span, 0.0, 1.0)
        θ = π * t

        σ_above = ν_above * e^2 / h
        σ_below = ν_below * e^2 / h
        σ_mid = (σ_above + σ_below) / 2
        δσ = (σ_above - σ_below) / 2

        σ_xy = σ_mid + δσ * cos(θ)
        σ_xx = abs(δσ * sin(θ))

        # Convert conductivity → resistivity
        denom = σ_xx^2 + σ_xy^2
        if denom > 1e-30
            R_xy[bi] = σ_xy / denom
            R_xx[bi] = σ_xx / denom
        else
            R_xy[bi] = B / (n_e * e)
            R_xx[bi] = 0.0
        end
    end

    return R_xy, R_xx
end
