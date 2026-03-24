# HANDOFF — FQHE Project

**Date:** 2026-03-24
**Status:** Phase A–C COMPLETE + numerical CG + improved extrapolation. ED pipeline extended to N=8, plot with ab initio gaps and experimental comparison.

## What was done (2026-03-24 session 2)

### 1. Numerical Clebsch-Gordan coefficients (replaces WignerSymbols.jl)
- **src/clebsch_gordan.jl**: Racah formula with precomputed log-factorials (Float64).
  Validated against WignerSymbols.jl exact rational CG: 16,791 coefficients tested,
  max relative error 1.8×10⁻¹¹ (enters quadratically in matrix elements → negligible).
- **Speedup**:
  - N=7 at ν=1/3: 170s → 0.3s (CG was bottleneck; now Lanczos-limited)
  - N=8 at ν=1/3: 30+ min → 18 min (dim=8512, Lanczos is genuine work)
  - N=9 at ν=1/3: impossible → tractable (dim=45207, ~hours)
- Removed `using WignerSymbols` from hamiltonian.jl. WignerSymbols.jl remains in
  Project.toml but is no longer used at runtime.

### 2. ED extended to N=8
- **Neutral gaps at ν=1/3** (Lz=0 sector, Haldane sphere):
  ```
  N=3: Δ=0.11899  twoS=6   dim=5      0.01s
  N=4: Δ=0.09353  twoS=9   dim=18     0.00s
  N=5: Δ=0.09312  twoS=12  dim=73     0.00s
  N=6: Δ=0.08203  twoS=15  dim=338    0.01s
  N=7: Δ=0.08083  twoS=18  dim=1656   0.3s   ← NEW
  N=8: Δ=0.08160  twoS=21  dim=8512   18min  ← NEW
  ```
- **Shell effect observed**: N=8 gap (0.0816) > N=7 gap (0.0808). This non-monotonicity
  is well-known on the Haldane sphere and motivates quadratic extrapolation.

### 3. Improved finite-size extrapolation
- **extrapolate_gap** now supports `order` kwarg (auto-selects quadratic for ≥4 points).
- Quadratic 1/N + 1/N² fit with N=3–8: Δ(∞) = 0.0852 = **82% of published 0.1036**.
  Previous linear 1/N from N=3–6: 0.048 = 46%. Almost doubled accuracy.
- Dense diag threshold raised from 1500 to 4000 (N=7 was 0.3s dense vs 30s Lanczos).

### 4. Plot script optimized
- **04_make_plot.jl**: Uses precomputed gap values (no live ED), generates in ~5s.
  Previously recomputed charge gaps from scratch (~minutes).
- Switched from charge gap to neutral gap as anchor (charge gap has sign issues
  for Jain fractions; neutral gap is well-behaved at all N).
- Caption updated to reflect method (quadratic extrapolation, N=3–8).

### 5. Tensor network survey
- **docs/tensor_network_FQHE_survey.md**: 427-line survey of TN methods for FQHE
  (2020–2026), covering iDMRG, MPS, PEPS, open-source codes, comparison with ED.
- Key finding: iDMRG on cylinders gives Δ(1/3) ≈ 0.101 (essentially the reference
  value), with no finite-size extrapolation. Best available thermodynamic-limit method.
- Open-source path: TeNPy (Python, v1.0 2024) or ITensors.jl (Julia).
  No standalone Julia FQHE-DMRG code exists.

## What was done (2026-03-24 session 1)

### 1. ED pipeline improvements
- **exact_diag.jl**: Fixed `_lowest_energy` sector scanning — proper parity stepping,
  VJ computed once, returns from first non-empty sector (rotational invariance).
- **Charge gap for ν=1/3**: Computed at N=3,4,5,6, extrapolated to N→∞ via Δ(N)=Δ(∞)+a/N.
  Result: Δ(1/3,∞) = 0.075 (published ref: 0.1036, ratio 0.72 — linear 1/N undershoots
  without quadratic correction).
- **Charge gap for Jain fractions**: NEGATIVE at all accessible N for ν=2/5, 3/7.
  Root cause: on small spheres (N≤8), the quasiparticle state at 2S−1 falls into a
  different topological sector. This is a known finite-size artifact for non-Laughlin fractions.

### 2. CF gap scaling
- **composite_fermion.jl**: Implemented `cf_gap_scaling(Δ_anchor, p; r=0.56)`.
  Geometric decay: Δ(p/(2p+1)) = Δ(1/3) × 0.56^(p−1).
  Calibrated from published Δ(2/5)/Δ(1/3). Matches p=3,4 within 10%.
- Jain principal sequence computed to p=7 (7/15), with particle-hole conjugates.

### 3. Transport model fix
- **transport.jl**: Added neighbor-aware plateau width capping. Each plateau's
  half-width is capped at 40% of the distance to its nearest neighbor. Without this,
  closely-spaced Jain fractions near ν=1/2 overlap and produce no R_xx peaks.

### 4. Plot updated
- Fractions shown: 1/3, 2/5, 3/7, 4/9, 5/11, 6/13, 7/15, 2/3, 3/5, 4/7, 5/9,
  6/11, 7/13, 1/5, plus integers 1–4.
- Fine staircase structure visible near ν=1/2 with R_xx peaks between all fractions.

### 5. Experimental data overlay
- **Source**: Wang et al., PNAS 120, e2314212120 (2023), CC-BY-4.0 (Zenodo record 10279058).
  GaAs 2D hole gas, n_h ≈ 1.55×10¹¹ cm⁻², ultra-high mobility.
- **Processing**: B rescaled by density ratio (×1.0/1.55 = ×0.645) to align filling factors
  with our electron simulation (n_e = 1.0×10¹¹ cm⁻²). Stitched low-field (0–9T) and
  high-field (9–16T) panels from Fig 4. Smoothed with 15-point moving average, downsampled
  to ~3600 points. R_xy from Fig 3a (1/3 plateau region).
- **Result**: Experimental R_xx overlay (gray) shows SdH oscillations and FQHE minima
  that align with theory predictions. R_xy confirms universal 1/3 plateau quantization.
- **Other datasets downloaded** (in data/, gitignored):
  - ETH Enkner et al. (Nature 2025) — GaAs 2DEG with cavity coupling, 193MB
  - Fu et al. (Nat Comms 2019) — confined 2DEG, nu=1–2 region only

## Error budget

| Step | Error | Status |
|------|-------|--------|
| Constants (CODATA) | < 1 ppm | Negligible |
| Pseudopotentials (Fano eq. 25) | < 1 ppm | Exact formula |
| Hamiltonian assembly | < 1e-10 | Exact |
| Lanczos diag | < 1e-10 | tol=1e-12 |
| **Finite-size extrapolation** | **18%** | **Dominant error** (quadratic, N=3–8) |
| CF scaling (r=0.56) | 7–10% | Empirical |
| LL mixing (not implemented) | 10–20% | Planned |
| Finite well width (not impl.) | 10–20% | Planned |
| Transport widths | ~30–50% | Phenomenological (by design) |

## Known issues / next steps

**Remaining gap to reference (82% → 100%):** The 18% shortfall comes from (a) shell
oscillations at small N and (b) limited data points for extrapolation. N=9+ would help
but takes hours. iDMRG would bypass this entirely.

**Highest-leverage improvements (ranked):**
1. iDMRG on cylinder via ITensors.jl — gives Δ(1/3) ≈ 0.101 directly, opens ν=5/2
2. LL mixing corrections (Faugno et al.) — 10–20% systematic shift
3. Finite well width corrections (Peterson et al.) — 10–20% on pseudopotentials
4. Chalker-Coddington network model for microscopic plateau widths
5. N=9–10 ED (hours of compute, diminishing returns vs DMRG)

## How to run

```bash
# Validate ED (should give ≈ 0.082035):
julia --project=. -e '
using FQHE
twoS = sphere_flux(6, 1//3)
VJ = coulomb_pseudopotentials(twoS)
basis = enumerate_fock_states(6, twoS; twoLz=0)
E0, gap = neutral_gap(basis, VJ)
println("N=6 ν=1/3: Δ_exc = $gap")
'

# Generate plot (~5s after JIT, uses precomputed gaps):
julia --project=. scripts/04_make_plot.jl
# → fqhe_plot.pdf, fqhe_plot.png

# Full gap pipeline (neutral + charge, ~20 min for N up to 8 at 1/3):
julia --project=. scripts/02_compute_gaps.jl
```

## Prerequisites
- Julia 1.10+ (tested with 1.12.5)
- `af` CLI (go binary)
- TIB VPN + Playwright for APS papers (already fetched)
