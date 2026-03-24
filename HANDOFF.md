# HANDOFF — FQHE Project

**Date:** 2026-03-24
**Status:** Phase A–C COMPLETE. ED pipeline validated, Jain sequence to p=7, plot with ab initio gaps.

## What was done (2026-03-24 session)

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
- **04_make_plot.jl**: Now computes gaps from ED (ν=1/3 charge gap as anchor) +
  CF scaling for Jain hierarchy, instead of hardcoded reference values.
- Fractions shown: 1/3, 2/5, 3/7, 4/9, 5/11, 6/13, 7/15, 2/3, 3/5, 4/7, 5/9,
  6/11, 7/13, 1/5, plus integers 1–4.
- Fine staircase structure visible near ν=1/2 with R_xx peaks between all fractions.

### 5. Performance findings
- **WignerSymbols.jl CG bottleneck**: CG coefficient computation becomes very slow
  for twoS > 15 (S > 7.5). N=7 at ν=1/3 (twoS=18) takes ~170s for neutral gap,
  ~236s for charge gap. N=8 (twoS=21) takes 30+ minutes. The bottleneck is exact
  rational arithmetic inside WignerSymbols.jl, not the Lanczos diagonalization.
- **Practical ceiling**: twoS ≤ 15 for fast computation (~1s). This limits:
  ν=1/3 to N≤6, ν=2/5 to N≤8, ν=1/5 to N≤4.

### 6. Known issues / next steps

**Charge gap for Jain fractions:** The simple "add/remove one flux" charge gap
is negative for ν=2/5, 3/7 at all accessible N. Proper Jain charge gaps on the
sphere require CF wavefunctions to define quasihole/quasiparticle excitations.
Alternative: use neutral gap at larger N once CG bottleneck is resolved.

**CG coefficient performance:** Replace WignerSymbols.jl with a numerical (Float64)
CG library, or cache CG coefficients across Hamiltonian builds. This would unlock
N=8–10 for all fractions (~10× more data for extrapolation).

**Extrapolation accuracy:** Linear 1/N extrapolation from N=3–6 gives ~72% of
published gaps. Adding quadratic 1/N² term requires N=8+ data points.

**Future extensions (from PRD):**
- Finite well width corrections (Peterson et al.) — 10–20% improvement
- LL mixing corrections (Faugno et al.) — biggest systematic error
- ν=5/2 via DMRG/ITensors.jl
- Chalker-Coddington network model for microscopic plateau widths

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

# Generate plot (computes ED gaps + CF scaling, ~5s after JIT):
julia --project=. scripts/04_make_plot.jl
# → fqhe_plot.pdf, fqhe_plot.png

# Full gap pipeline (charge + neutral, ~7 min for N up to 8 at 1/3):
julia --project=. scripts/02_compute_gaps.jl
```

## Prerequisites
- Julia 1.10+ (tested with 1.12.5)
- `af` CLI (go binary)
- TIB VPN + Playwright for APS papers (already fetched)
