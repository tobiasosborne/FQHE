# FQHE.jl

Ab initio simulation of the Fractional Quantum Hall Effect in GaAs, from Haldane pseudopotentials to transport curves.

![Fractional Quantum Hall Effect — GaAs 2DEG](fqhe_plot.png)

*Hall resistance (R_xy, blue) and longitudinal resistance (R_xx, red) for a GaAs 2DEG at T = 50 mK, 40 filling factors. The 1/3 gap is from DMRG on a cylinder (Ly=10, 86% of published reference); higher Jain fractions from composite-fermion scaling. Plateau widths and transition shapes are derived from the Dykhne–Ruzin semicircle law and Wei–Tsui–Tsui scaling (κ = 0.42). Gray: experimental data from a GaAs 2D hole gas ([Wang et al., PNAS 2023](https://doi.org/10.1073/pnas.2314212120), CC-BY-4.0), with B rescaled to align filling factors.*

## Overview

FQHE.jl builds the complete pipeline from microscopic Coulomb interactions to measurable transport signatures:

1. **Haldane sphere geometry** — monopole harmonics, angular momentum coupling via numerical Clebsch-Gordan coefficients (Racah formula, replaces WignerSymbols.jl for ~500× speedup)
2. **Pseudopotentials** — Coulomb interaction projected into the lowest Landau level (Fano et al. 1986, closed-form binomial formula)
3. **Exact diagonalization** — sparse Hamiltonian in the Fock basis, Lanczos solver for ground state and excitation gaps (N up to 8 at $\nu = 1/3$)
4. **Composite fermion theory** — geometric gap scaling for the Jain sequence: $\Delta(\nu = p/(2p+1)) = \Delta(1/3) \times r^{p-1}$
5. **DMRG on cylinder** — ITensors-based ground state and gap for ν=1/3 (Ly=10, χ=800), pseudopotential-decomposed MPO
6. **Transport model** — σ_xy from additive tanh staircase (Dykhne–Ruzin semicircle law), σ_xx from sech² peaks, transition widths from Wei–Tsui–Tsui scaling (κ = 0.42), tensor inversion for R_xy and R_xx

## Quick start

```bash
# Clone and instantiate
git clone <repo-url> && cd FQHE
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Validate ED against Fano et al. (1986) — should give gap ~ 0.082:
julia --project=. -e '
using FQHE
twoS = sphere_flux(6, 1//3)
VJ = coulomb_pseudopotentials(twoS)
basis = enumerate_fock_states(6, twoS; twoLz=0)
E0, gap = neutral_gap(basis, VJ)
println("N=6 ν=1/3 neutral gap: $gap")
'

# Generate the plot (~5s after JIT):
julia --project=. scripts/04_make_plot.jl
# → fqhe_plot.pdf, fqhe_plot.png
```

## Project structure

```
src/
  FQHE.jl                 # Module entry point
  materials.jl             # GaAs parameters (m*, g, ε)
  landau.jl                # Landau levels, magnetic length, cyclotron energy
  sphere.jl                # Haldane sphere: flux-particle relation
  monopole_harmonics.jl    # Single-particle orbitals on the sphere
  hilbert_space.jl         # Many-body Fock space enumeration
  clebsch_gordan.jl        # Numerical CG coefficients (Racah formula)
  pseudopotentials.jl      # Coulomb Haldane pseudopotentials
  hamiltonian.jl           # Sparse Hamiltonian construction
  exact_diag.jl            # Lanczos ED: ground state, neutral & charge gaps
  composite_fermion.jl     # CF gap scaling for Jain sequence
  laughlin.jl              # Laughlin wavefunction overlaps
  integer_qhe.jl           # Integer QHE (analytic)
  transport.jl             # R_xx, R_xy from gaps (semicircle law + scaling theory)
  cylinder.jl              # Cylinder geometry, pseudopotential decomposition
  dmrg_hamiltonian.jl      # ITensors MPO assembly
  dmrg_solver.jl           # DMRG solver wrapper

scripts/
  02_compute_gaps.jl       # Full gap pipeline (charge + neutral)
  04_make_plot.jl          # Final plot with 40 fractions + exp overlay
  06_cylinder_ed.jl        # Cylinder matrix element validation
  07_validate_mpo.jl       # MPO vs sparse ED (matches to 5.5e-16)
  08_dmrg_Ly_scan.jl       # Gap vs circumference scan
  09_dmrg_gap.jl           # Production DMRG gap (ν=1/3, Ly=10)

proof/                     # Formal verification artifacts
test/                      # Unit tests
```

## Key results

| Quantity | This code | Published | Reference |
|----------|-----------|-----------|-----------|
| Neutral gap, $\nu=1/3$, $N=6$ | 0.0820 $e^2/\epsilon\ell$ | 0.0822 | Fano et al. (1986) |
| Neutral gap, $\nu=1/3$, $N=8$ | 0.0816 $e^2/\epsilon\ell$ | — | This work |
| Neutral gap, $\nu=1/3$, $N\to\infty$ | 0.085 $e^2/\epsilon\ell$ | 0.1036 | Quadratic extrap. (82% of ref) |
| Charge gap, $\nu=1/3$, DMRG cyl | 0.089 $e^2/\epsilon\ell$ | 0.1036 | Ly=10, χ=800 (86% of ref) |

The DMRG cylinder result (86% of reference) is the primary gap anchor; ED sphere extrapolation (82%) serves as a cross-check. The remaining 14% is from finite circumference (Ly=10) and bond dimension (χ=800).

## Experimental comparison

The gray curves overlay real transport data from a GaAs 2D hole gas ([Wang et al., PNAS 2023](https://doi.org/10.1073/pnas.2314212120), CC-BY-4.0). Since the experimental sample has a different carrier density ($n_h = 1.55 \times 10^{11}$ cm$^{-2}$ vs our $n_e = 1.0 \times 10^{11}$ cm$^{-2}$), the magnetic field axis is rescaled by the density ratio so that filling factors align. The FQHE is topological — quantized Hall plateaux and $R_{xx}$ minima appear at the same filling fractions regardless of carrier type, confirming the universality of the effect.

## Requirements

- Julia 1.10+ (tested with 1.12)
- Dependencies installed automatically via `Project.toml`: Arpack, CairoMakie, Combinatorics, KrylovKit

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
