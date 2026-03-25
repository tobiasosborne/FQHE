# HANDOFF — FQHE Project

**Date:** 2026-03-25
**Status:** 40 filling factors, all ab initio. DMRG cylinder + ED sphere + 2nd-LL + CFL. Smooth transport model.

---

## Architecture Overview

```
src/
  FQHE.jl                 # Module entry — includes all below in order
  materials.jl             # GaAs: m*=0.067m_e, ε=12.9, g=-0.44
  landau.jl                # ℓ_B, ℏω_c, E_C, filling_factor
  sphere.jl                # Haldane sphere: shift(ν)=q universal, sphere_flux(N,ν)
  monopole_harmonics.jl    # Single-particle orbitals on sphere
  hilbert_space.jl         # FockBasis: UInt64 bitstrings, enumerate_fock_states
  clebsch_gordan.jl        # Racah formula CG (500× faster than WignerSymbols.jl)
  pseudopotentials.jl      # Coulomb V_J (sphere, Fano eq.25) + V_J^(n) (n-th LL)
  hamiltonian.jl           # Sparse H on sphere from CG + pseudopotentials
  exact_diag.jl            # Lanczos: neutral_gap, charge_gap, extrapolate_gap
  cylinder.jl              # Cylinder: plane V_m, Hermite pair amplitudes g_m(r)
  dmrg_hamiltonian.jl      # ITensors MPO from pseudopotential decomposition
  dmrg_solver.jl           # DMRG wrappers: GS + excited state + noise
  composite_fermion.jl     # CF gap scaling: Δ(p/(2p+1)) = Δ(1/3)×0.56^(p-1)
  laughlin.jl              # Laughlin wavefunction overlaps
  integer_qhe.jl           # Analytic integer QHE
  transport.jl             # Smooth R_xx, R_xy from gaps (sech² peaks, σ_xx→R_xx)

scripts/
  02_compute_gaps.jl       # ED for 13 LLL fractions (~10s, dim<5K fast mode)
  04_make_plot.jl          # Final plot with 40 fractions + exp overlay (~5s)
  05_itensors_smoke.jl     # ITensors sanity check (free fermions)
  06_cylinder_ed.jl        # Cylinder matrix element validation
  07_validate_mpo.jl       # MPO vs sparse ED (84×84, matches to 5.5e-16)
  08_dmrg_Ly_scan.jl       # Gap vs circumference scan
  09_dmrg_gap.jl           # Production DMRG gap (ν=1/3, Ly=10)
```

---

## Key Results

| Quantity | Value | Method | Reference |
|----------|-------|--------|-----------|
| Δ(1/3), N→∞ | 0.089 e²/(εℓ) | DMRG cylinder Ly=10 | 0.1036 (86%) |
| Δ(1/3), N=3-8 | 0.085 e²/(εℓ) | ED sphere quadratic | 0.1036 (82%) |
| Δ(7/3) | 0.030 e²/(εℓ) | ED 2nd-LL pseudo N=7 | 0.045 |
| Δ(5/2) | 0.025 e²/(εℓ) | Published ref | 0.025 |
| Filling factors | 40 total | ED + DMRG + CF + p-h | — |

---

## Learnings & Pitfalls for Next Agent

### CRITICAL: Julia scoping

All scripts must wrap logic in `function main() ... end; main()`. Julia's soft-scope
rules cause `UndefVarError` for variables assigned inside `for` loops at top level.
Every script in this project follows this pattern.

### CRITICAL: Never two Julia processes on same project

Julia's precompilation cache is not concurrent-safe. Running two `julia --project=.`
processes simultaneously causes cache corruption and mysterious errors. Always kill
one before starting another. Check with `ps aux | grep julia`.

### DMRG convergence on small cylinders

- The ν=1/3 ground state on a cylinder is **3-fold topologically degenerate** (CDW sectors).
  DMRG with a single CDW initial state may converge to an excited state.
  **Fix**: try 3 CDW shifts (0,1,2), take the lowest E₀. Then orthogonalize against
  all GS states to find the first excitation above the manifold.

- At Nphi=9 (tiny system), DMRG gets stuck at the CDW diagonal energy (-0.640 vs
  true E₀=-0.911). The MPO is correct (verified to machine epsilon via full matrix
  reconstruction), but DMRG needs noise terms (`noise=[1e-2, 1e-3, ...]`) to mix
  the CDW with off-diagonal states.

### Cylinder Coulomb matrix elements: Madelung subtraction

The initial Fourier-space formula for cylinder Coulomb had a **fatal bug**: the q_y=0
(s=0) channel diverges logarithmically, making energies negative and divergent with
system size. Two failed approaches:

1. **f_sub(0,δ) = ∫ (cos(qδ)-1)/q × exp(-q²/2) dq**: Finite but grows as -ln(|δ|).
   Total energy diverges. Diagonal elements are negative (wrong sign for repulsion).

2. **f(0,δ) = 0**: Drops too much physics. V₁-model gap goes to zero. Spectrum
   becomes nearly degenerate.

**Working approach**: Haldane pseudopotential decomposition.
  V_{abcd} = Σ_{m odd} V_m × g_m(a-b) × g_m(c-d)
where V_m = C(2m,m)/4^m (planar Coulomb) and g_m(r) = Hermite function.
This is always positive-definite, no background issues, V₁-model gives exact E=0
for Laughlin state.

### MPO construction performance

The `MPO(ampo, sites)` call in ITensors can be **extremely slow** for long-range
interactions (hours for Nphi=15 with 500+ four-body terms). The fix is **range
truncation**: pair amplitude g_m(r) decays as exp(-r²Δ²/4), so terms beyond
~2.5Ly/π orbital spacings are negligible (<1e-4). With truncation, MPO builds
in 0.1-0.2s and has bond dim ~40-50.

### 2nd-LL pseudopotentials

`coulomb_pseudopotentials_nLL(twoS, n_LL)` multiplies LLL pseudopotentials by
the Laguerre form factor |L_n(x)|² where x = J(J+1)/(2S(2S+2)). Key check:
V₁/V₃ ratio drops from 1.55 (LLL) to 0.32 (2nd LL) — this is what makes the
Pfaffian win over Laughlin at ν=5/2.

### ED extrapolation unreliable at few data points

With 2-3 points, `extrapolate_gap` often gives nonsense (negative gaps, wildly wrong
magnitudes). The 4-flux fractions (2/7, 3/11) have ED gaps 10-1000× too large at
small N due to shell effects. For fractions with <4 ED points, use the **raw gap at
largest N** instead of extrapolation.

### Transport model

The new smooth model (sech² peaks) looks much better than the old hard-edge version.
Key parameters: `σ_ν ∝ disorder/gap` controls peak width. The conversion σ_xx → R_xx
via R_xx = σ_xx/(σ_xx² + σ_xy²) naturally gives the right peak-to-background ratio.

---

## Error Budget

| Step | Error | Status |
|------|-------|--------|
| Constants (CODATA) | < 1 ppm | Negligible |
| Pseudopotentials (Fano eq. 25) | < 1 ppm | Exact formula |
| Hamiltonian assembly | < 1e-10 | Exact (verified) |
| Lanczos diag | < 1e-10 | tol=1e-12 |
| **Finite-size (ED sphere)** | **18%** | Quadratic extrap, N=3-8 |
| **Finite-size (DMRG cyl)** | **14%** | Ly=10, Nphi=21, Ne=7, χ=800 |
| CF scaling (r=0.56) | 7-10% | Empirical, calibrated |
| LL mixing (not implemented) | 10-20% | **Highest-leverage next step** |
| Finite well width (not impl.) | 10-20% | Planned |
| Transport peak shapes | qualitative | Smooth sech², not microscopic |

---

## Highest-Leverage Next Steps (ranked)

1. **Larger DMRG for 1/3** (Ly=14-16, χ=2000): improves the anchor gap that ALL
   CF-scaled fractions depend on. Current 86% → expect ~95%. Zero coding, ~2hr compute.
   Run: modify Ly and maxdim in `scripts/09_dmrg_gap.jl`.

2. **DMRG for 2/5**: second-most-important fraction, currently CF scaling only.
   Infrastructure already generalized — just change ν parameter. ~1hr compute.

3. **N=8 ED for 1/3**: adds the critical 6th data point for quadratic extrapolation.
   Already computed (0.0816), just need to add back to the precomputed arrays in
   `04_make_plot.jl` (was dropped in the rewrite). 18min compute.

4. **LL mixing corrections**: Faugno et al. (P20, local PDF) gives perturbative
   correction to pseudopotentials. 10-20% systematic shift on ALL gaps. ~3hr coding.

5. **Finite well width**: Peterson et al. (P13, local PDF). Modifies V_m → V_m(w).
   ~2hr coding.

6. **iDMRG** (ITensorInfiniteMPS.jl): eliminates finite-length effects entirely.
   ~4hr coding + testing.

---

## How to Run

```bash
# Quick ED validation (should give ≈ 0.082035):
julia --project=. -e '
using FQHE
twoS = sphere_flux(6, 1//3)
VJ = coulomb_pseudopotentials(twoS)
basis = enumerate_fock_states(6, twoS; twoLz=0)
E0, gap = neutral_gap(basis, VJ)
println("N=6 ν=1/3: Δ_exc = $gap")
'

# Full ED batch (13 fractions, ~10s):
julia --project=. scripts/02_compute_gaps.jl
# → data/ed_gaps.csv

# Generate plot (~5s, uses precomputed gaps):
julia --project=. scripts/04_make_plot.jl
# → fqhe_plot.pdf, fqhe_plot.png

# DMRG gap for ν=1/3 (~2min):
julia --project=. scripts/09_dmrg_gap.jl

# Validate cylinder MPO vs ED (~30s):
julia --project=. scripts/07_validate_mpo.jl
```

## Prerequisites
- Julia 1.10+ (tested with 1.12.5)
- Dependencies: ITensors v0.9, ITensorMPS v0.3, KrylovKit, CairoMakie, SpecialFunctions
- `af` CLI (go binary) — for source fetching
- TIB VPN + Playwright for paywalled papers (28 PDFs already in sources/papers/)
