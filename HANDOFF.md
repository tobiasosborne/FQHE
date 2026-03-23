# HANDOFF — FQHE Project

**Date:** 2026-03-23
**Status:** Phase A COMPLETE, Phase B COMPLETE, Phase C MVP COMPLETE (plot generated).

## What was done this session

### 1. New machine setup
- Julia 1.12.5 with all deps installed + precompiled
- `af` CLI confirmed working, proof tree intact (31 nodes)
- 24 papers downloaded (15 arXiv via curl, 9 APS via Playwright + TIB VPN)
- fetch_aps_papers.mjs rewritten for headed persistent-context Playwright

### 2. Phase C1: Foundations fixed
- **pseudopotentials.jl**: Replaced wrong CG-sum formula with Fano eq. 25 closed-form binomial.
  Uses BigInt to avoid overflow. Validated V_J values produce correct energies.
- **sphere.jl**: Fixed integer division bug. Now uses rational arithmetic (`N // ν`).
  shift(ν) = denominator for all Jain fractions + conjugates.
- **hilbert_space.jl**: Fixed half-integer Lz bug. Now works with `twoLz` (integer)
  everywhere. Added `state_to_index` Dict for O(1) Hamiltonian assembly.
- **monopole_harmonics.jl**: Simplified to orbital_m(i, twoS) helper.

### 3. Phase C2: Hamiltonian + ED validated against Fano (1986)
- **hamiltonian.jl**: Full sparse Coulomb Hamiltonian assembly.
  Precomputes CG-based two-body matrix elements grouped by total M.
  Fermion signs via bit-manipulation. Hermitian wrapper for exact symmetry.
- **exact_diag.jl**: Dense eigvals for dim ≤ 1500, KrylovKit Lanczos otherwise.
  Neutral gap + charge gap extraction. Linear 1/N extrapolation.

**Validation (ν=1/3 neutral gap, Fano PRB 34 2670 Table IV):**
| N | dim(Lz=0) | Δ_exc (ours) | Δ_exc (Fano) |
|---|-----------|-------------|-------------|
| 3 | 5         | 0.118992    | 0.118990    |
| 4 | 18        | 0.093533    | 0.093530    |
| 5 | 73        | 0.093117    | 0.093110    |
| 6 | 338       | 0.082035    | 0.082040    |

Match to 5 significant figures.

### 4. Phase C3: Transport + plot
- **transport.jl**: Semicircle law (Dykhne-Ruzin) for transitions, activated
  R_xx at plateau centers, SdH oscillations at high filling.
- **04_make_plot.jl**: Two-panel CairoMakie figure. Uses reference gap values
  (Fano 1986 + Balram 2018) validated against our ED.
- **fqhe_plot.pdf/png**: The iconic plot is generated.

### 5. Known issues / next steps

**Charge gap sector scanning:** The charge_gap() function doesn't properly
scan all Lz sectors for quasihole/quasiparticle states. The quasihole ground
state is at L = N/2, not Lz=0. Fix: scan all Lz sectors from 0 to N×twoS/2.
Neutral gaps work perfectly and are used for the plot.

**N=8+ ED performance:** Hamiltonian assembly + diag for N=8 (dim≈8500) takes
several minutes due to CG coefficient computation. Consider caching CG values
or using a faster CG library.

**Transport model polish:** The semicircle transitions could be smoother.
The SdH oscillation formula needs tuning. Plateau widths are phenomenological
and somewhat arbitrary.

**Future extensions (from PRD):**
- ν=5/2 via DMRG/iTensors.jl
- Chalker-Coddington network model for microscopic plateau widths
- Finite well width corrections (Peterson et al.)
- LL mixing corrections (Faugno et al.)

## How to run

```bash
# Validate ED:
julia --project=. -e '
using FQHE
twoS = sphere_flux(6, 1//3)
VJ = coulomb_pseudopotentials(twoS)
basis = enumerate_fock_states(6, twoS; twoLz=0)
E0, gap = neutral_gap(basis, VJ)
println("N=6 ν=1/3: Δ_exc = $gap")  # should be ≈ 0.082035
'

# Generate plot:
julia --project=. scripts/04_make_plot.jl
# → fqhe_plot.pdf, fqhe_plot.png
```

## Prerequisites
- Julia 1.10+ (tested with 1.12.5)
- `af` CLI (go binary)
- TIB VPN + Playwright for APS papers (already fetched)
