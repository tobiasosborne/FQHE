# CLAUDE.md — Fractional Quantum Hall Effect: The Iconic Plot from First Principles

## Mission

Compute the Hall resistance R_xy and longitudinal resistance R_xx as functions of
magnetic field B for a 2D electron gas, reproducing the iconic FQHE plot from
first principles. No fits to experimental data. The only inputs are fundamental
constants (ℏ, e, ε₀), material parameters (effective mass m*, dielectric constant ε),
and sample parameters (electron density n_e, temperature T, mobility μ_e).

The output is a two-panel plot:
- Top: R_xy / (h/e²) vs B showing integer and fractional plateaus
- Bottom: R_xx vs B showing Shubnikov-de Haas oscillations and FQHE minima

## Constraints

- 100% Julia. No Python, no Fortran, no wrappers.
- Julia 1.10+. Pin all deps in Project.toml.
- No single file > 250 lines.
- Every physical constant and formula must cite a source in `./sources/`.
- The exact diagonalisation must be genuinely ab initio in the lowest Landau level:
  Coulomb interaction only, no phenomenological parameters for the gaps.
- The transport model may be semi-phenomenological (this is stated honestly in the output).
- Target: GaAs/AlGaAs heterostructure (the standard system).

## Directory Structure

```
fqhe-iconic-plot/
├── CLAUDE.md
├── README.md
├── Project.toml
├── Manifest.toml
├── sources/
│   ├── papers/
│   ├── data/                  # Experimental data for comparison (Tsui, Stormer, etc.)
│   └── bibliography.bib
├── src/
│   ├── FQHE.jl               # Main module
│   ├── landau.jl             # Landau levels, filling factor, cyclotron energy
│   ├── sphere.jl             # Haldane sphere geometry
│   ├── monopole_harmonics.jl # Single-particle basis on sphere
│   ├── hilbert_space.jl      # Many-body Fock space construction at fixed N, 2S
│   ├── pseudopotentials.jl   # Haldane pseudopotentials V_m for Coulomb
│   ├── hamiltonian.jl        # Sparse Hamiltonian assembly from pseudopotentials
│   ├── exact_diag.jl         # Lanczos/Arnoldi for ground state + gap
│   ├── composite_fermion.jl  # Jain CF theory: gap predictions, CF filling map
│   ├── laughlin.jl           # Laughlin wavefunction overlaps (validation)
│   ├── integer_qhe.jl        # IQHE: Landau level gaps (analytic, no ED needed)
│   ├── transport.jl          # Gaps → R_xx, R_xy vs B
│   └── materials.jl          # GaAs parameters, physical constants
├── test/
│   ├── runtests.jl
│   ├── test_pseudopotentials.jl  # Compare V_m to published tables
│   ├── test_laughlin.jl          # Ground state overlap with Laughlin state
│   ├── test_gaps.jl              # Compare gaps to published ED results
│   ├── test_integer.jl           # R_xy = h/(ne²) exactly at ν = n
│   └── test_transport.jl         # Sanity checks on transport model
├── scripts/
│   ├── 01_fetch_sources.jl
│   ├── 02_compute_gaps.jl        # Run ED for all target fractions
│   ├── 03_build_transport.jl     # Assemble R_xx, R_xy curves
│   └── 04_make_plot.jl           # Final figure
├── notebooks/
│   └── fqhe_plot.jl              # Pluto notebook: interactive version
└── docs/
    ├── physics_notes.md          # Derivations, conventions, subtleties
    ├── toolchain_plan.md
    └── gaps_table.md             # Computed gaps for each fraction
```

## Phase A: Source Acquisition

### A1. Papers (priority order)

| ID | Reference | Content needed |
|----|-----------|---------------|
| P1 | Haldane, "Fractional quantization of the Hall effect" PRL 51, 605 (1983) | Sphere geometry, pseudopotentials |
| P2 | Laughlin, "Anomalous quantum Hall effect" PRL 50, 1395 (1983) | Laughlin wavefunction |
| P3 | Jain, "Composite-fermion approach for FQHE" PRL 63, 199 (1989) | CF theory, filling map |
| P4 | Fano, Ortolani & Colombo, "Configuration-interaction calculations on the FQHE" PRB 34, 2670 (1986) | Early ED gaps, reference values |
| P5 | Morf & d'Ambrumenil, "Disorder in the QHE" PRL / PRB (various) | Disorder + transport |
| P6 | Tsui, Stormer & Gossard, "Two-dimensional magnetotransport in the extreme quantum limit" PRL 48, 1559 (1982) | Original FQHE discovery data |
| P7 | Willett et al., "Observation of an even-denominator quantum number" PRL 59, 1776 (1987) | ν=5/2, extended data |
| P8 | Jain, *Composite Fermions* (Cambridge, 2007) | THE textbook. Chapters 4-7 are essential |
| P9 | Chakraborty & Pietiläinen, *The Quantum Hall Effects* (Springer, 1995) | ED methodology, gap tables |
| P10 | Dykhne & Ruzin, "Theory of the fractional QHE: the two-phase model" PRB 50, 2369 (1994) | Semicircle law for transport |
| P11 | Girvin, "The QHE: novel excitations and broken symmetries" (Les Houches lectures) arXiv:cond-mat/9907002 | Best pedagogical intro |
| P12 | Störmer, "Nobel lecture: The FQHE" Rev. Mod. Phys. 71, 875 (1999) | Experimental context + iconic plots |
| P13 | Peterson, Park & Jolicoeur, "Finite-layer thickness corrections to the CF theory" PRB 78, 155308 (2008) | Finite-width corrections to gaps |
| P14 | Balram, Barkeshli, Rudner & Jain, "Detailed study of the Jain fractions" PRB (various) | Modern gap calculations |
| P15 | Haldane & Rezayi, "Finite-size studies of the incompressible state..." PRL 54, 237 (1985) | Sphere ED, gap extrapolation |

### A2. Key numerical data to extract or download

| ID | Data | Source |
|----|------|--------|
| D1 | Haldane pseudopotentials V_m for Coulomb in LLL on sphere | Derive from formula; validate against P1, P9 |
| D2 | Published ED gaps at ν=1/3, 2/5, 3/7, 2/3, 1/5 | P4, P9, P15 — need N-dependence for extrapolation |
| D3 | GaAs material parameters: m*=0.067 m_e, ε=12.9, g*=-0.44 | Standard, e.g. Ando, Fowler & Stern Rev. Mod. Phys. 54 (1982) |
| D4 | Experimental R_xx, R_xy traces for comparison overlay | P6, P12 — digitise from figures if necessary |
| D5 | CF gap predictions (analytic) | P3, P8 Chapter 7 |

### A3. Retrieval

Same procedure as Th-229 prompt: arXiv first, TIB VPN for paywalled.
Girvin lectures (P11) are freely on arXiv. Jain's textbook (P8) is the one
item you may need to buy or access via TIB. The Les Houches lectures are
actually sufficient for implementation if the book is unavailable.

Citation graph: run Semantic Scholar forward-citation on P1, P2, P3.
Identify any post-2020 papers with improved gap values or transport models.
Cap at ~40 papers total.

## Phase B: Toolchain Planning

### B1. The Physics Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│ INPUT: n_e, T, B_min, B_max, material params                       │
│                                                                     │
│ Step 1: ν(B) = n_e h / (eB)           [landau.jl — trivial]       │
│                                                                     │
│ Step 2: For each target fraction ν = p/q:                           │
│   2a. Map to sphere: N electrons, 2S = N/ν - S_shift   [sphere.jl]│
│   2b. Build LLL Hilbert space at (N, 2S, L_z=0)  [hilbert_space.jl]│
│   2c. Compute Haldane pseudopotentials V_m      [pseudopotentials.jl]│
│   2d. Assemble sparse H in Fock basis            [hamiltonian.jl]   │
│   2e. Lanczos → E_0 (ground state energy)         [exact_diag.jl]  │
│   2f. Compute gap: Δ = E_0(N, 2S+1) - E_0(N, 2S) [exact_diag.jl] │
│   2g. Finite-size extrapolation: Δ(N) → Δ(∞)      [exact_diag.jl] │
│                                                                     │
│ Step 3: Gaps in natural units → physical units:                     │
│   Δ_phys = Δ × (e²/4πεε₀ℓ_B) where ℓ_B = √(ℏ/eB)  [landau.jl]  │
│                                                                     │
│ Step 4: Build R_xy(B), R_xx(B)                    [transport.jl]   │
│   - Plateau at each ν=p/q with width ∝ Δ                           │
│   - R_xx ~ exp(-Δ/2kT) at plateau centers                          │
│   - Transition regions via semicircle law or phenomenological model │
│                                                                     │
│ Step 5: Plot                                       [make_plot.jl]  │
└─────────────────────────────────────────────────────────────────────┘
```

### B2. Julia packages needed

| Need | Package | Notes |
|------|---------|-------|
| Sparse matrices | `SparseArrays` (stdlib) | For Hamiltonian |
| Eigensolvers | `Arpack.jl` or `KrylovKit.jl` | Lanczos for lowest eigenvalues. KrylovKit preferred (pure Julia) |
| Wigner/CG coefficients | `WignerSymbols.jl` or `CGcoefficient.jl` | For monopole harmonics and matrix elements |
| Combinatorics | `Combinatorics.jl` | Enumerate Fock states |
| Plotting | `CairoMakie.jl` | Publication-quality two-panel figure |
| Units | `Unitful.jl` | Optional but nice for dimensional analysis |
| Notebooks | `Pluto.jl` | Interactive version |
| Parallelism | `Threads` or `Distributed` | ED for different N, ν can run in parallel |

### B3. What must be built

**1. Monopole harmonics on Haldane sphere** (`monopole_harmonics.jl`)

The single-particle states in the LLL on a sphere with monopole strength S are:
  Y_{S,S,m}(θ,φ) for m = -S, -S+1, ..., S
giving 2S+1 orbitals. For N electrons at filling ν = N/(2S+1) on the sphere,
the relation between N and S encodes the shift:
  2S = N/ν - S_shift
where S_shift is a topological quantum number (S_shift = q for Laughlin 1/q).

For the LLL, we don't actually need the explicit harmonics for ED. We need
the two-body matrix elements, which reduce to Haldane pseudopotentials.
But we need the angular momentum algebra to construct Fock states with
definite total L_z (and ideally total L for identification of the ground state).

**2. Fock space construction** (`hilbert_space.jl`)

- Enumerate all ways to place N fermions in 2S+1 orbitals (LLL)
- Each orbital labeled by m = -S, ..., S
- Fock state: sorted tuple of occupied m-values
- Total L_z = sum of occupied m-values
- Restrict to L_z = 0 sector (ground state lives here for rotationally invariant H)
- Store as integer bitstrings for efficiency
- Dimension check: C(2S+1, N) total, restricted to L_z=0 sector

Typical dimensions for target calculations:

| ν | N | 2S | Full dim | L_z=0 dim | Feasible? |
|---|---|-----|----------|-----------|-----------|
| 1/3 | 6 | 15 | 8,008 | ~500 | Easy |
| 1/3 | 8 | 21 | 203,490 | ~9,000 | Easy |
| 1/3 | 10 | 27 | 3,124,550 | ~100,000 | Moderate |
| 1/3 | 12 | 33 | ~900M | too big without L | Need total L |
| 2/5 | 6 | 13 | 1,716 | ~130 | Easy |
| 2/5 | 8 | 18 | 43,758 | ~2,400 | Easy |
| 2/5 | 10 | 23 | 1,144,066 | ~50,000 | Easy |
| 3/7 | 6 | 12 | 924 | ~80 | Trivial |

So N=8-10 is the sweet spot: large enough for reasonable extrapolation,
small enough for laptop-scale Lanczos. N=12 at ν=1/3 needs total-L
decomposition or more memory.

**3. Haldane pseudopotentials** (`pseudopotentials.jl`)

For Coulomb interaction in the LLL on a sphere of radius R = √S ℓ_B:

  V_m = ∫ |ψ_m(r₁₂)|² V(r₁₂) d²r₁₂

where ψ_m is the relative wavefunction of two particles with relative
angular momentum m (m odd for fermions).

For Coulomb V(r) = e²/(εr) on the sphere, the pseudopotentials are:

  V_m = (e²/εℓ_B) × Σ_l [ P_l(cos θ) × |⟨S,S,m|l⟩|² ] / (2√S)

More explicitly, using the known formula (Fano et al. 1986):

  V_m = Σ_{l=0}^{2S} (2l+1) [ C^{S,S,l}_{m/2,-m/2,0} ]² / (2√(S + 1/2))

where C are Clebsch-Gordan coefficients. This is eq. (3.43) in Jain's book.

- Compute V_m for m = 1, 3, 5, ..., 2S-1 (odd only, fermions)
- Validate against published tables (P4, P9, or Jain Table 4.1)
- These are the ONLY input to the many-body Hamiltonian

**4. Hamiltonian assembly** (`hamiltonian.jl`)

The two-body Hamiltonian in second quantisation:

  H = (1/2) Σ_{m1,m2,m3,m4} ⟨m1,m2|V|m3,m4⟩ c†_m1 c†_m2 c_m4 c_m3

The matrix elements decompose via CG coefficients into pseudopotentials:

  ⟨m1,m2|V|m3,m4⟩ = δ_{m1+m2, m3+m4} × Σ_m V_m ×
    C^{S,S,2S-m}_{m1,m2,...} × C^{S,S,2S-m}_{m3,m4,...}

(exact formula in Haldane 1983 or Jain Ch. 4)

Assembly strategy:
- Loop over pairs of Fock states |α⟩, |β⟩ in the L_z=0 sector
- For each pair, find the two orbitals that differ
- Compute matrix element using pseudopotentials + CG
- Store in sparse CSC format
- H is real symmetric → use Hermitian wrapper for eigensolver

**5. Exact diagonalisation** (`exact_diag.jl`)

- Use `KrylovKit.eigsolve` for lowest k eigenvalues of sparse H
- k=3-5 is sufficient (need ground state + first excitation)
- The ground state at the Laughlin/Jain filling has total L=0
- The gap is defined as: Δ = E_0(quasihole) - E_0(ground)
  In practice: compute E_0 at 2S (ground) and 2S+1 (one quasihole added),
  or identify the magnetoroton gap in the excitation spectrum at fixed 2S.
- Multiple approaches exist; simplest is the "neutral gap":
  Δ_neutral = E_1(L>0) - E_0(L=0) at fixed N, 2S
- Thermodynamic extrapolation: compute for N = 6, 8, 10, (12 if feasible)
  and fit Δ(N) = Δ(∞) + c/N

**6. Laughlin wavefunction overlap** (`laughlin.jl`, for validation)

For ν=1/3, the Laughlin state on the sphere is known analytically.
Compute its expansion in the Fock basis and check overlap with ED ground state.
Should be > 0.99 for N ≤ 10. This is the smoking-gun validation that the
ED code is correct.

Implementation: the Laughlin state Ψ = Π_{i<j} (u_i v_j - u_j v_i)³
where u = cos(θ/2)e^{iφ/2}, v = sin(θ/2)e^{-iφ/2}.
Expand by Monte Carlo sampling or by direct polynomial expansion for small N.
Monte Carlo is easier: sample random configurations, evaluate Ψ and the Fock
state, accumulate overlap. Or: Jack polynomial representation (Bernevig & Haldane).

For a first pass, skip this and validate via gap values instead.
Implement Laughlin overlap only if gap values disagree with published results.

**7. Composite fermion gap predictions** (`composite_fermion.jl`)

As a cross-check independent of ED:

  Δ_CF(ν = p/(2p+1)) ≈ ℏω_c* = ℏeB*/(m_CF c)

where B* = B - 2n_e φ₀ is the effective field and m_CF is the CF effective mass.
Jain gives m_CF ≈ 0.6-1.0 √B[T] × m_e for GaAs (fitted from numerics).

This gives a quick estimate of all gaps without doing ED. Use it to:
- Sanity-check ED gaps (should be same order of magnitude)
- Fill in gaps for fractions where ED is too expensive
- Estimate gap hierarchy (which fractions are most visible)

**8. Transport model** (`transport.jl`)

This is the semi-phenomenological layer. Be honest about it in the output.

Model: at each B, compute ν(B) = n_e h/(eB). The transport coefficients are:

**At plateau centers (ν near p/q):**
  R_xy = h/(ν_plateau × e²)
  R_xx = R_xx^peak × exp(-Δ(ν)/(2kT))

**Between plateaus (transition regions):**
  Option A (simplest): linear interpolation in R_xy, Gaussian peaks in R_xx
  Option B (better): semicircle law (Dykhne-Ruzin):
    σ_xx² + (σ_xy - (n+1/2)e²/h)² = (e²/(2h))²
    in the transition region between plateaus at ν=n and ν=n+1

  Convert σ → R via R_xy = σ_xy/(σ_xx² + σ_xy²), etc.

**Plateau width:**
  The plateau at filling ν extends over a range ΔB ∝ Δ(ν)/ℏω_c × B_center
  Sharper plateaus (wider ΔB) for larger gaps.
  Phenomenological: plateau exists where the Fermi level sits in the gap,
  i.e. for |ν(B) - p/q| < Δ/(ℏω_c) approximately.
  Disorder broadens Landau levels by Γ ~ ℏ/τ; plateau width also depends on Γ.
  Use mobility μ_e to set Γ = ℏ/(m* μ_e × e B / m*) ... simplify.

For the plot to look right, the key parameters are:
- n_e ≈ 1-3 × 10¹¹ cm⁻² (typical GaAs 2DEG)
- T ≈ 50-300 mK (FQHE only visible at low T)
- μ_e ≈ 10⁵-10⁷ cm²/Vs (high-mobility sample)
- B range: 0 to ~15 T

**Target fractions to include:**

Principal Jain sequence ν = p/(2p+1):
  1/3, 2/5, 3/7, 4/9, 5/11, ...

Conjugate sequence ν = 1 - p/(2p+1):
  2/3, 3/5, 4/7, ...

Second Jain sequence ν = p/(2p-1):
  (maps to ν > 1/2: 2/3, 3/5 again — particle-hole conjugates)

Also: ν = 1/5 (Laughlin), ν = 2/3, and integers ν = 1, 2, 3, ...

Minimum set for a convincing plot: ν = 1, 2, 3 (integer) plus
1/3, 2/5, 3/7, 2/3, 3/5, 4/7, 1/5. That's 10 fractions.
Each fraction except integers needs an ED gap computation.

## Phase C: Implementation

### C1. Order of implementation

```
Step 1: materials.jl
         → GaAs constants: m*, ε, g*. Physical constants. ℓ_B(B), ℏω_c(B).
         → Test: ℓ_B(1 T) ≈ 25.6 nm, ℏω_c(1 T)/k_B ≈ 20 K for GaAs.

Step 2: landau.jl
         → ν(B, n_e), B_ν(ν, n_e), cyclotron energy, Coulomb energy scale e²/(εℓ_B).
         → Test: ν=1/3 at B ≈ 10 T for n_e = 1.5×10¹¹ cm⁻².

Step 3: pseudopotentials.jl
         → V_m for Coulomb on sphere. Needs CG coefficients.
         → Test: reproduce Table 4.1 of Jain or Table in Fano et al. (P4).
         → This is the critical validation: if V_m are wrong, everything downstream is wrong.

Step 4: hilbert_space.jl
         → Enumerate Fock states, L_z sectors, bit-string representation.
         → Test: dimension of L_z=0 sector matches combinatorial formula.
         → Test: for N=3, 2S=5 (ν=1/3 on sphere), dim(L_z=0) should be small
           enough to verify by hand.

Step 5: hamiltonian.jl
         → Sparse H from pseudopotentials + Fock states.
         → Test: H is Hermitian. Trace matches expected value.
         → Test: for N=3 at ν=1/3, the ground state is non-degenerate at L=0
           and the gap is known analytically (3-body problem).

Step 6: exact_diag.jl
         → Lanczos via KrylovKit. Ground state energy + gap.
         → Test: ν=1/3, N=6: gap should be ≈ 0.1 e²/(εℓ_B) (order of magnitude).
         → Test: compare to published values in P4, P9 for multiple N.

Step 7: composite_fermion.jl
         → CF gap estimates for all target fractions.
         → Test: hierarchy agrees with ED (1/3 > 2/5 > 3/7 ...).

Step 8: integer_qhe.jl
         → Analytic: Δ = ℏω_c for integer filling.
         → Trivial but needed for the complete plot.

Step 9: transport.jl
         → Gaps → R_xx(B), R_xy(B).
         → Test: R_xy = h/e² exactly at ν=1 plateau center.
         → Test: R_xx → 0 at plateau centers, peaked between plateaus.
         → Test: R_xy is monotonically increasing between plateaus.

Step 10: 04_make_plot.jl
         → The final figure. Two panels, shared B axis.
         → Annotate each plateau with its filling fraction.
         → Overlay experimental data from P6/P12 if digitised.
```

### C2. The ED inner loop (performance-critical)

The bottleneck is Step 5+6 for the larger Hilbert spaces. Key optimisations:

- **Bit-string Fock states:** represent each state as a UInt64 or UInt128
  (sufficient for 2S+1 ≤ 128 orbitals). Use bit operations for creation/annihilation.
- **Sparse assembly:** only store nonzero matrix elements. For 2-body H,
  the number of nonzeros per row is O(N² × N_orb) which is manageable.
- **Symmetry:** exploit L_z conservation (already done by sector restriction).
  For extra speed, decompose into total-L sectors using a projection operator
  or the method of Haldane & Rezayi. This is optional for N ≤ 10.
- **KrylovKit.eigsolve:** use `:SR` (smallest real) to get ground state.
  Typically converges in < 100 iterations for these matrix sizes.

Expected timings on a modern laptop (single core):

| ν | N | L_z=0 dim | Assembly | Lanczos | Total |
|---|---|-----------|----------|---------|-------|
| 1/3 | 8 | ~9,000 | < 1s | < 1s | seconds |
| 1/3 | 10 | ~100,000 | ~10s | ~30s | minute |
| 2/5 | 8 | ~2,400 | < 1s | < 1s | seconds |
| 2/5 | 10 | ~50,000 | ~5s | ~10s | seconds |

All target fractions at N=8 should complete in under a minute total.
N=10 adds maybe 10 minutes. This is laptop-scale.

### C3. Finite-size extrapolation

For each fraction, compute the gap at N = 6, 8, 10 (and 12 if feasible).
The thermodynamic gap is:

  Δ(N) = Δ(∞) + a/N + b/N²

Fit a 2-parameter model (Δ(∞) and a) using the 3 data points.
Report Δ(∞) in units of e²/(εℓ_B) and also in Kelvin at B = 10 T.

Write results to `docs/gaps_table.md`:

```markdown
| ν | N=6 | N=8 | N=10 | Δ(∞) [e²/(εℓ_B)] | Δ(∞) at 10T [K] | Published |
|---|-----|-----|------|-------------------|-----------------|-----------|
| 1/3 | ... | ... | ... | ~0.1 | ~5 | 0.1 (various) |
| 2/5 | ... | ... | ... | ~0.06 | ~3 | 0.06 |
```

### C4. The transport model (be honest)

The transport model converts gaps into the resistance curves. Document clearly
in the output what is ab initio and what is phenomenological:

**Ab initio:**
- Plateau values R_xy = h/(νe²) — exact, topological
- Gap magnitudes Δ(ν) — from ED, no free parameters
- Gap hierarchy — which fractions are visible, from ED
- Temperature dependence R_xx ∝ exp(-Δ/2kT) — activated transport, standard

**Phenomenological (state honestly):**
- Plateau widths — depend on disorder, parametrised by mobility
- Transition region shape — semicircle law or interpolation
- Overall R_xx scale — depends on sample geometry

The plot should have a caption or annotation that says:
"Plateau positions and gap hierarchy from ab initio exact diagonalisation.
Plateau widths and transition shapes are semi-phenomenological (see docs/physics_notes.md)."

## Phase D: Output

### D1. The figure

Two-panel plot, vertically stacked, shared horizontal B-axis (0 to 15 T).

**Top panel (R_xy):**
- y-axis: R_xy in units of h/e² (so plateaus sit at 1, 1/2, 1/3, 2/5, ...)
- Staircasing plateaus with fractional labels
- Classical Hall line R_xy = B/(n_e e) as dashed reference
- Colour: blue or black

**Bottom panel (R_xx):**
- y-axis: R_xx in Ω (or arbitrary units)
- Shubnikov-de Haas oscillations at high ν
- Deep minima (→ 0) at FQHE fractions
- Peaks between plateaus
- Colour: red or black

**Annotations:**
- Label each visible fraction: 1/3, 2/5, 3/7, 2/3, 3/5, 1, 2, 3
- Temperature and density in corner: T = 100 mK, n_e = 1.5 × 10¹¹ cm⁻²
- "Theoretical prediction" clearly stated

Style: clean, publication-quality. Use CairoMakie with LaTeX labels.
Export as PDF and PNG.

### D2. Sample parameters for the canonical plot

Use parameters close to original Tsui-Stormer experiment:
- n_e = 1.5 × 10¹¹ cm⁻² (→ ν=1/3 at B ≈ 12.4 T)
- T = 100 mK
- μ_e = 10⁶ cm²/Vs (high mobility)
- B range: 0 to 15 T

Also generate a version at T = 50 mK to show more fractions appearing.

## Validation Criteria

1. **Pseudopotentials V_m** match published tables to 6+ significant figures.
2. **ED gap at ν=1/3** matches published values (≈ 0.1 e²/(εℓ_B)) to within
   finite-size error bars across N = 6, 8, 10.
3. **Laughlin overlap** > 0.99 at ν=1/3 for N ≤ 8 (if implemented).
4. **Gap hierarchy** correct: Δ(1/3) > Δ(2/5) > Δ(3/7) > Δ(2/3) > Δ(1/5).
5. **Integer plateaus** at exact values R_xy = h/(ne²).
6. **Plot is visually recognisable** as the FQHE plot to any condensed matter physicist.
7. **Every number traceable** to a source or to our own ED calculation.

## Anti-patterns

- Do NOT implement a full Lanczos from scratch. Use KrylovKit.jl. It is pure Julia,
  well-tested, and faster than anything you'll write in a day.
- Do NOT try to go beyond N=12 for ν=1/3. That's a supercomputer calculation.
  N=10 is sufficient for a convincing extrapolation.
- Do NOT model disorder microscopically. The phenomenological approach is standard
  and honest. A microscopic disorder model would require weeks and add little to the plot.
- Do NOT attempt ν=5/2 (Moore-Read / Pfaffian state). This is a separate, harder problem
  requiring even-denominator physics. Stick to odd-denominator Jain fractions.
- Do NOT optimise the Hamiltonian assembly before validating at N=6. Get it right first.
- Do NOT use Float32 for the Hamiltonian. The energy differences that give the gap
  are small compared to the total energy. Use Float64.

## Time Budget (agentic coding)

| Phase | Estimate |
|-------|----------|
| A: Sources + citation graph | 1-2 hours |
| B: Toolchain plan | 30 min |
| C1-C3: Pseudopotentials + Fock space + validation | 2-3 hours |
| C4-C6: Hamiltonian + ED + gap extraction | 2-3 hours |
| C7-C8: CF cross-check + integer QHE | 30 min |
| C9: Transport model | 1-2 hours |
| C10: Final plot | 1 hour |
| D: Documentation + cleanup | 30 min |
| **Total** | **~1 day** |

## Success Metric

```
julia --project=. scripts/04_make_plot.jl
```

produces `fqhe_plot.pdf`: the iconic fractional quantum Hall effect plot,
computed from nothing but the Coulomb interaction and fundamental constants,
with annotated fractions, and a caption stating clearly what is ab initio
and what is phenomenological.
