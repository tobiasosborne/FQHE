# transport.jl — Gaps → R_xx, R_xy vs B
#
# Smooth transport model: R_xy uses sigmoid transitions between plateaus,
# R_xx is derived from the dissipative conductivity σ_xx which peaks at
# transitions. Each plateau transition has width ∝ disorder/gap, giving
# naturally narrow peaks for strong fractions and broad ones for weak.
#
# Sources:
#   P10: Dykhne & Ruzin, PRB 50, 2369 (1994) — semicircle law
#   P11: Girvin lectures, §2-3

export compute_transport, FractionData

"""Data for one QH fraction used by the transport model."""
struct FractionData
    ν::Rational{Int}
    gap_dimless::Float64     # gap in units of e²/(ε ℓ_B), 0 for gapless
end

"""
    compute_transport(B_range, n_e, fractions, T; mat, μ_e) → (R_xy, R_xx)

Compute R_xy and R_xx vs B using a smooth transport model.

R_xy transitions smoothly between quantized plateaus via sigmoid functions.
R_xx is derived from the dissipative conductivity σ_xx, which peaks at
plateau-to-plateau transitions with width ∝ disorder/gap.
"""
function compute_transport(B_range::AbstractVector, n_e::Real,
                           fractions::Vector{FractionData}, T::Real;
                           mat::MaterialParams=GaAs, μ_e::Real=1e2)
    h = CONSTANTS.h;  e = CONSTANTS.e;  kB = CONSTANTS.k_B
    RK = h / e^2

    # Sort fractions by ν (descending) → increasing B
    sorted = sort(fractions, by=f -> Float64(f.ν), rev=true)

    # Precompute plateau data
    ν_f  = Float64[Float64(f.ν) for f in sorted]
    B_f  = [n_e * h / (e * ν) for ν in ν_f]

    # Gap in Kelvin for each fraction
    gap_K = Float64[]
    for (i, f) in enumerate(sorted)
        if f.gap_dimless > 0
            push!(gap_K, f.gap_dimless * coulomb_energy(B_f[i], mat) / kB)
        else
            push!(gap_K, cyclotron_energy(B_f[i], mat) / kB)
        end
    end

    # Disorder broadening in Kelvin
    Γ_K = CONSTANTS.ħ * e / (mat.m_star * μ_e) / kB

    # Transition width in filling factor: σ_ν ∝ Γ/Δ (wider for weaker fractions)
    σ_ν = Float64[]
    for i in eachindex(sorted)
        # Width controlled by ratio of disorder to gap, clamped
        σ = clamp(0.015 * max(Γ_K, T) / max(gap_K[i], 0.1), 0.002, 0.08)
        push!(σ_ν, σ)
    end

    R_xy = similar(B_range, Float64)
    R_xx = similar(B_range, Float64)

    for (bi, B) in enumerate(B_range)
        if B <= 0.01
            R_xy[bi] = 0.0;  R_xx[bi] = 0.0;  continue
        end

        ν_B = n_e * h / (e * B)

        # ── R_xy: smooth interpolation between plateaus ──
        # Start from classical Hall line, add plateau corrections
        R_xy_classical = B / (n_e * e)

        # Find the two nearest plateaus bracketing ν_B
        # and compute weighted plateau value
        R_xy_val = R_xy_classical
        total_weight = 0.0

        for i in eachindex(ν_f)
            # Sigmoid weight: how "locked" are we to this plateau?
            # Uses a steep tanh profile; steepness ∝ 1/σ_ν
            x = (ν_B - ν_f[i]) / σ_ν[i]
            # Plateau weight: Gaussian-like, strong near ν_f, falls off with σ_ν
            w = exp(-x^2 / 2)
            # Thermal suppression: at high T, plateaus weaken
            thermal = gap_K[i] > 0 ? min(1.0, gap_K[i] / (2T)) : 10.0
            w *= min(thermal, 1.0)
            total_weight += w
        end

        if total_weight > 0.01
            # Blend plateau values weighted by proximity
            R_blend = 0.0
            for i in eachindex(ν_f)
                x = (ν_B - ν_f[i]) / σ_ν[i]
                w = exp(-x^2 / 2)
                thermal = gap_K[i] > 0 ? min(1.0, gap_K[i] / (2T)) : 10.0
                w *= min(thermal, 1.0)
                R_plateau = RK / ν_f[i]
                R_blend += w * R_plateau
            end
            R_blend /= total_weight

            # Smooth crossover between classical and plateau
            lock_strength = min(total_weight, 1.0)
            R_xy_val = lock_strength * R_blend + (1 - lock_strength) * R_xy_classical
        end

        R_xy[bi] = R_xy_val

        # ── R_xx: dissipative conductivity peaks at transitions ──
        # σ_xx peaks where ν is between plateaus (transition regions)
        # Uses the semicircle law: σ_xx² + (σ_xy - σ_mid)² = (δσ)²

        σ_xx_val = 0.0

        # For each adjacent pair of plateaus, compute the transition contribution
        for i in 1:length(ν_f)-1
            ν_hi = ν_f[i]
            ν_lo = ν_f[i+1]
            ν_mid = (ν_hi + ν_lo) / 2

            # Transition region width in ν: sum of the two plateau σ values
            σ_trans = σ_ν[i] + σ_ν[i+1]

            # How far are we from this transition midpoint?
            x = (ν_B - ν_mid) / σ_trans

            # Semicircle peak shape: sech² gives smooth, naturally-shaped peaks
            peak = 1.0 / cosh(1.5 * x)^2

            # Peak amplitude: proportional to the conductance step
            δσ = (ν_hi - ν_lo) / (2 * RK)
            σ_xx_val += δσ * peak
        end

        # Also add SdH oscillation for ν > highest fraction
        if ν_B > ν_f[1] + 3 * σ_ν[1]
            τ = mat.m_star * μ_e / e
            ωc = e * B / mat.m_star
            X = 2π^2 * kB * T / (CONSTANTS.ħ * ωc)
            dingle = exp(-π / (ωc * τ))
            thermal = X > 0.01 ? X / sinh(X) : 1.0
            σ_xx_val += (ν_B / RK) * 0.08 * (1.0 - 4.0 * dingle * thermal * cos(2π * ν_B))
        end

        # Convert σ_xx to R_xx: R_xx = σ_xx / (σ_xx² + σ_xy²)
        σ_xy = 1.0 / R_xy_val
        denom = σ_xx_val^2 + σ_xy^2
        R_xx[bi] = denom > 1e-30 ? σ_xx_val / denom : 0.0

        # Ensure R_xx ≥ 0
        R_xx[bi] = max(R_xx[bi], 0.0)
    end

    return R_xy, R_xx
end
