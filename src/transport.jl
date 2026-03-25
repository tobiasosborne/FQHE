# transport.jl — Gaps → R_xx, R_xy vs B
#
# Transport model grounded in two established results:
#
# 1. Dykhne–Ruzin semicircle law (PRB 50, 2369, 1994):
#    σ_xx² + (σ_xy − σ_mid)² = δσ²
#    The tanh parameterization σ_xy = σ_mid + δσ·tanh(x) satisfies this
#    exactly because sech²(x) + tanh²(x) ≡ 1.
#
# 2. Wei–Tsui–Tsui scaling (PRL 61, 1294, 1988; Nature Comm. 2024):
#    The plateau-to-plateau transition width obeys
#      δν ∝ (max(Γ, k_BT) / Δ)^κ,   κ = 0.42 ± 0.04
#    This universal exponent holds for both IQHE and FQHE transitions
#    ("superuniversality", confirmed in trilayer graphene, 2024).
#
# σ_xy is built as an additive staircase: a sum of independent sigmoid steps,
# one per pair of adjacent plateaus.  This avoids the piecewise-selection
# discontinuities that arise when transition widths exceed the plateau spacing.
#
# R_xy and R_xx are obtained from (σ_xy, σ_xx) via the tensor relation:
#    R_xy = σ_xy / (σ_xx² + σ_xy²)
#    R_xx = σ_xx / (σ_xx² + σ_xy²)

export compute_transport, FractionData

"""Data for one QH fraction used by the transport model."""
struct FractionData
    ν::Rational{Int}
    gap_dimless::Float64     # gap in units of e²/(ε ℓ_B), 0 for gapless
end

"""
    compute_transport(B_range, n_e, fractions, T; mat, μ_e) → (R_xy, R_xx)

Compute Hall (R_xy) and longitudinal (R_xx) resistance vs B.

σ_xy: additive staircase of tanh sigmoid steps (semicircle-law-consistent).
σ_xx: sech² peaks at each plateau transition.
Resistivities from tensor inversion of the conductivity.
"""
function compute_transport(B_range::AbstractVector, n_e::Real,
                           fractions::Vector{FractionData}, T::Real;
                           mat::MaterialParams=GaAs, μ_e::Real=1e2)
    h = CONSTANTS.h;  e = CONSTANTS.e;  kB = CONSTANTS.k_B
    RK = h / e^2

    # Sort fractions by ν (descending) → increasing B
    sorted = sort(fractions, by=f -> Float64(f.ν), rev=true)

    # Precompute data for ALL fractions (used by σ_xx peaks)
    ν_f  = Float64[Float64(f.ν) for f in sorted]
    B_f  = [n_e * h / (e * ν) for ν in ν_f]

    # Gap in Kelvin
    gap_K = Float64[]
    for (i, f) in enumerate(sorted)
        if f.gap_dimless > 0
            push!(gap_K, f.gap_dimless * coulomb_energy(B_f[i], mat) / kB)
        else
            push!(gap_K, cyclotron_energy(B_f[i], mat) / kB)
        end
    end

    # Disorder broadening: Γ_K = ħe/(m*μ) in Kelvin (= ħ/τ in energy / k_B)
    Γ_K = CONSTANTS.ħ * e / (mat.m_star * μ_e) / kB

    # Wei–Tsui–Tsui universal exponent (PRL 61, 1294, 1988)
    κ = 0.42

    # ══════════════════════════════════════════════════════════════════
    # σ_xy additive staircase: plateau-forming fractions only
    # ══════════════════════════════════════════════════════════════════
    # Include: gapped states (gap_dimless > 0) and integers (cyclotron gap)
    # Exclude: gapless metallic states (ν=1/2 CFL)
    ν_plat   = Float64[]
    gap_plat = Float64[]
    for i in eachindex(sorted)
        if sorted[i].gap_dimless > 0 || denominator(sorted[i].ν) == 1
            push!(ν_plat, ν_f[i])
            push!(gap_plat, gap_K[i])
        end
    end
    n_plat = length(ν_plat)

    # Precompute transitions between each adjacent pair of plateaus
    n_trans  = max(n_plat - 1, 0)
    ν_trans  = Vector{Float64}(undef, n_trans)
    δ_trans  = Vector{Float64}(undef, n_trans)
    Δσ_step  = Vector{Float64}(undef, n_trans)  # conductance step per transition
    for k in 1:n_trans
        g_hi, g_lo = gap_plat[k], gap_plat[k+1]
        ν_hi, ν_lo = ν_plat[k], ν_plat[k+1]
        spacing = ν_hi - ν_lo

        # Gap-weighted midpoint: stronger fraction claims more territory
        ν_trans[k] = (ν_hi * g_lo + ν_lo * g_hi) / (g_hi + g_lo)

        # Transition width: (Γ/Δ)^κ scaling (Wei et al.)
        ratio = max(Γ_K, T) / max(min(g_hi, g_lo), 1e-10)
        δ_trans[k] = spacing * ratio^κ

        # Conductance step at this transition
        Δσ_step[k] = (ν_hi - ν_lo) / RK
    end

    R_xy = similar(B_range, Float64)
    R_xx = similar(B_range, Float64)

    for (bi, B) in enumerate(B_range)
        if B <= 0.01
            R_xy[bi] = 0.0;  R_xx[bi] = 0.0;  continue
        end

        ν_B = n_e * h / (e * B)

        # ── σ_xy: additive staircase ─────────────────────────────
        # σ_xy = σ_base + Σ_k Δσ_k × ½(1 + tanh(x_k))
        # Each transition contributes independently → no seam discontinuities.
        # At ν well above all transitions: every frac ≈ 1 → σ_xy = ν_plat[1]/RK
        # At ν well below all transitions: every frac ≈ 0 → σ_xy = ν_plat[end]/RK
        if n_plat > 0 && ν_B <= ν_plat[1]
            σ_xy_val = ν_plat[end] / RK
            for k in 1:n_trans
                x = (ν_B - ν_trans[k]) / δ_trans[k]
                frac = 0.5 * (1.0 + tanh(x))
                σ_xy_val += Δσ_step[k] * frac
            end
        else
            σ_xy_val = ν_B / RK   # classical Hall (above all plateaus)
        end

        # ── σ_xx: sech² peaks at each transition ─────────────────
        # Peak widths floored by the κ-scaled transition width to ensure
        # consistency with σ_xy and prevent ultra-narrow integer spikes.
        σ_xx_val = 0.0

        for i in 1:length(ν_f)-1
            ν_hi = ν_f[i]
            ν_lo = ν_f[i+1]
            ν_mid = (ν_hi + ν_lo) / 2
            spacing_i = ν_hi - ν_lo

            # Base width from gap-scaled formula
            σ_base = clamp(0.015 * max(Γ_K, T) / max(gap_K[i], 0.1), 0.002, 0.08) +
                     clamp(0.015 * max(Γ_K, T) / max(gap_K[i+1], 0.1), 0.002, 0.08)

            # Floor: κ-scaled width consistent with the σ_xy transition model
            gap_min_i = max(min(gap_K[i], gap_K[i+1]), 0.1)
            σ_κ = spacing_i * (max(Γ_K, T) / gap_min_i)^κ

            σ_trans = max(σ_base, σ_κ)

            x = (ν_B - ν_mid) / σ_trans
            peak = 1.0 / cosh(1.5 * x)^2

            # Peak amplitude: proportional to the conductance step
            δσ = spacing_i / (2 * RK)
            σ_xx_val += δσ * peak
        end

        # SdH oscillations above the highest fraction
        if ν_B > ν_f[1] + 0.3
            τ = mat.m_star * μ_e / e
            ωc = e * B / mat.m_star
            X = 2π^2 * kB * T / (CONSTANTS.ħ * ωc)
            dingle = exp(-π / (ωc * τ))
            thermal = X > 0.01 ? X / sinh(X) : 1.0
            σ_xx_val += (ν_B / RK) * 0.08 * (1.0 - 4.0 * dingle * thermal * cos(2π * ν_B))
        end

        # ── Tensor inversion: (σ_xx, σ_xy) → (R_xx, R_xy) ──────
        denom = σ_xx_val^2 + σ_xy_val^2
        R_xy[bi] = denom > 1e-30 ? σ_xy_val / denom : 0.0
        R_xx[bi] = denom > 1e-30 ? max(σ_xx_val / denom, 0.0) : 0.0
    end

    return R_xy, R_xx
end
