# transport.jl — Gaps → R_xx, R_xy vs B
#
# Semi-phenomenological transport model.
# Ab initio: plateau positions (h/νe²), gap magnitudes, hierarchy.
# Phenomenological: plateau widths (via mobility), transition shapes.
#
# Sources:
#   P10: Dykhne & Ruzin, PRB 50, 2369 (1994) — semicircle law
#   P11: Girvin lectures, §2-3

export compute_transport, FractionData

"""Data for one QH fraction used by the transport model."""
struct FractionData
    ν::Rational{Int}
    gap_dimless::Float64     # gap in units of e²/(ε ℓ_B), 0 for integer QHE
end

"""
    compute_transport(B_range, n_e, fractions, T; mat, μ_e) → (R_xy, R_xx)

Compute R_xy and R_xx vs B.
fractions: list of FractionData for all plateaus to include.
μ_e in m²/(V·s).  T in Kelvin.
"""
function compute_transport(B_range::AbstractVector, n_e::Real,
                           fractions::Vector{FractionData}, T::Real;
                           mat::MaterialParams=GaAs, μ_e::Real=1e2)
    h = CONSTANTS.h;  e = CONSTANTS.e;  kB = CONSTANTS.k_B
    RK = h / e^2

    # Sort fractions by ν (descending) → increasing B
    sorted = sort(fractions, by=f -> Float64(f.ν), rev=true)

    # Precompute: for each fraction, its B-center, physical gap, and plateau half-width in ν
    ν_list  = Float64[];  B_list = Float64[];  gap_K = Float64[];  hw = Float64[]
    for f in sorted
        ν_f = Float64(f.ν)
        B_f = n_e * h / (e * ν_f)
        # Physical gap in Kelvin
        if f.gap_dimless > 0
            Δ_K = f.gap_dimless * coulomb_energy(B_f, mat) / kB
        else
            # Integer QHE: gap = ℏω_c
            Δ_K = cyclotron_energy(B_f, mat) / kB
        end
        # Plateau half-width in ν: scales with gap/disorder
        Γ_K = CONSTANTS.ħ * e / (mat.m_star * μ_e) / kB   # disorder broadening in K
        w = clamp(0.03 * Δ_K / max(Γ_K, 1e-10), 0.003, ν_f * 0.12)
        push!(ν_list, ν_f);  push!(B_list, B_f)
        push!(gap_K, Δ_K);   push!(hw, w)
    end

    # Cap half-widths so adjacent plateaus don't swallow transition regions.
    # Without this, closely-spaced Jain fractions near ν=1/2 overlap and
    # produce no R_xx peaks between them.
    for fi in eachindex(ν_list)
        min_dist = Inf
        for fj in eachindex(ν_list)
            fi == fj && continue
            min_dist = min(min_dist, abs(ν_list[fi] - ν_list[fj]))
        end
        hw[fi] = min(hw[fi], 0.4 * min_dist)
    end

    R_xy = similar(B_range, Float64)
    R_xx = similar(B_range, Float64)

    for (bi, B) in enumerate(B_range)
        if B <= 0.01
            R_xy[bi] = 0.0;  R_xx[bi] = 0.0;  continue
        end

        ν_B = n_e * h / (e * B)

        # ── Check each plateau ──
        on_plateau = false
        for fi in eachindex(sorted)
            ν_f = ν_list[fi]
            δν = abs(ν_B - ν_f)
            w_f = hw[fi]
            if δν < w_f
                on_plateau = true
                R_xy[bi] = RK / ν_f
                # Activated R_xx: R_xx ∝ exp(-Δ/2kT) at plateau center
                R_xx[bi] = gap_K[fi] > 0 ? 200.0 * exp(-gap_K[fi] / (2T)) : 0.0
                break
            end
        end
        on_plateau && continue

        # ── Transition region: find bounding plateaus ──
        ν_hi = nothing;  ν_lo = nothing;  fi_hi = 0;  fi_lo = 0
        for fi in eachindex(ν_list)
            if ν_list[fi] > ν_B + hw[fi]
                ν_hi = ν_list[fi];  fi_hi = fi
            end
            if ν_list[fi] < ν_B - hw[fi] && ν_lo === nothing
                ν_lo = ν_list[fi];  fi_lo = fi
            end
        end

        if ν_hi !== nothing && ν_lo !== nothing
            # Semicircle transition between ν_hi and ν_lo
            σ_hi = ν_hi / RK;  σ_lo = ν_lo / RK
            σ_mid = (σ_hi + σ_lo) / 2
            δσ = (σ_hi - σ_lo) / 2

            # Map ν_B to parameter θ ∈ [0, π]
            edge_hi = ν_hi - hw[fi_hi]
            edge_lo = ν_lo + hw[fi_lo]
            t = clamp((edge_hi - ν_B) / max(edge_hi - edge_lo, 1e-10), 0.0, 1.0)
            θ = π * t

            σ_xy = σ_mid + δσ * cos(θ)
            σ_xx = abs(δσ) * sin(θ)

            denom = σ_xx^2 + σ_xy^2
            R_xy[bi] = denom > 1e-30 ? σ_xy / denom : B / (n_e * e)
            R_xx[bi] = denom > 1e-30 ? σ_xx / denom : 0.0
        elseif ν_B > ν_list[1]
            # Above highest fraction: classical Hall + SdH oscillations
            R_xy[bi] = B / (n_e * e)
            ωc_τ = e * B * μ_e / (CONSTANTS.ħ * e / (mat.m_star * μ_e) > 0 ? 1.0 : 1.0)
            τ = mat.m_star * μ_e / e
            ωc = e * B / mat.m_star
            X = 2π^2 * kB * T / (CONSTANTS.ħ * ωc)
            dingle = exp(-π / (ωc * τ))
            thermal = X > 0.01 ? X / sinh(X) : 1.0
            R_xx[bi] = abs((B / (n_e * e)) * 0.05 *
                       (1.0 - 4.0 * dingle * thermal * cos(2π * ν_B)))
        else
            # Below lowest fraction: classical fallback
            R_xy[bi] = B / (n_e * e)
            R_xx[bi] = 0.0
        end
    end

    return R_xy, R_xx
end
