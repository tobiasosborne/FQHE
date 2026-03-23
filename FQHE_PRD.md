# PRD: FQHE From First Principles

**Project:** Fractional Quantum Hall Effect — Iconic Plot from First Principles
**Author:** Tobias Osborne / Claude
**Date:** 2026-03-23
**Status:** Draft v1

---

## 1. Vision

Compute the Hall resistance R_xy(B) and longitudinal resistance R_xx(B) for a GaAs/AlGaAs 2DEG from nothing but the Coulomb interaction and fundamental constants. Produce the iconic two-panel FQHE plot that any condensed matter physicist would recognise.

The deeper goal: build a modular, extensible pipeline that functions as a **test harness for quantum many-body methods**. Each computational stage has a typed contract. Researchers can swap in new methods (tensor networks, improved pseudopotentials, microscopic disorder models) and immediately see the effect on a physically measurable observable.

## 2. Success criteria

### 2.1 Minimum viable plot (MVP)

A two-panel figure (R_xy and R_xx vs B, 0–15 T) showing:

- Integer plateaus at ν = 1, 2, 3 at exact quantised values R_xy = h/(νe²)
- Fractional plateaus at ν = 1/3, 2/5, 3/7, 2/3, 3/5, 1/5 with correct hierarchy
- R_xx minima at plateau centres, peaks between plateaus
- Shubnikov–de Haas oscillations at high filling
- Clear annotation of each fraction

The staircase shape is recognisable. Absolute gap values within ~20% of published ED results. Plateau widths are phenomenological (acknowledged in caption).

### 2.2 Quantitative validation targets

| Observable | Target | Source |
|---|---|---|
| Pseudopotential V_J at ν = 1/3, N = 6 | Match Fano et al. Table to 6 sig. fig. | P04 |
| ED gap Δ(1/3), N = 8 | ≈ 0.1 e²/(εℓ) | P04, P15 |
| Extrapolated Δ(1/3, ∞) | 0.1036 ± 0.005 e²/(εℓ) | P04 p.2678 |
| Gap hierarchy | Δ(1/3) > Δ(2/5) > Δ(3/7) | ED output |
| Laughlin overlap (if implemented) | > 0.99 for N ≤ 8 | P04 Table I |
| R_xy at ν = 1 | h/e² = 25812.807 Ω exactly | Topological |
| ℓ_B at 1 T | 25.66 nm | CODATA |
| ℏω_c/k_B at 1 T (GaAs) | ≈ 20 K | Standard |

### 2.3 Extended success (post-MVP)

- Overlay of digitised experimental data (Störmer Figure 18) on computed curves
- Compiled activation gap table compared with Du et al. / Pan et al. measurements
- Community benchmark dataset released (gaps + digitised experimental traces)
- ν = 5/2 via tensor network methods (ITensors.jl / iDMRG on cylinder)

## 3. Architecture

### 3.1 Pipeline stages

The computation is a pipeline of five stages. Each stage has a typed input/output contract. Methods can be swapped at any stage without affecting downstream stages.

```
Stage 1: (material, geometry, B) → single-particle basis
Stage 2: (basis, interaction) → pseudopotentials / matrix elements
Stage 3: (matrix elements, Fock space) → many-body Hamiltonian
Stage 4: (Hamiltonian) → spectrum → gap(ν)
Stage 5: (gaps, disorder params, T, B_range) → R_xy(B), R_xx(B)
```

### 3.2 Critical interface: GapTable

The GapTable is the central data structure. Everything upstream (however computed) reduces to a dictionary {ν ↦ Δ} that the transport model consumes. This is where the "multiscale" character lives.

```julia
struct GapTable
    gaps::Dict{Rational{Int}, Float64}  # ν => Δ in units of e²/(εℓ)
    method::String                       # provenance tag
    N_particles::Dict{Rational{Int}, Vector{Int}}  # system sizes used
    extrapolated::Bool                   # was finite-size extrapolation applied?
end
```

### 3.3 Method swap patterns

Three kinds of method swap must be supported:

**Single-box replacement.** Swap one stage, keep everything else. Example: replace ED (stage 4) with DMRG. Same input (pseudopotentials + geometry), same output (GapTable).

**Multi-box collapse.** A method replaces several stages. Example: composite fermion theory goes from (ν, B, material) directly to an approximate gap, bypassing stages 2–4. The output is still a GapTable.

**New box insertion.** An extra computation step between existing stages. Example: finite-well-width corrections (Peterson et al.) insert between stages 1 and 2, modifying pseudopotentials. Or Chalker-Coddington network model inserts into stage 5, replacing phenomenological plateau widths with a microscopic localisation calculation.

### 3.4 Pipeline registry

Each pipeline configuration is named and produces a plot. The baseline is the unit test.

```julia
pipelines = Dict(
    :baseline    => Pipeline(HaldaneSphere(), BareCoulomb(), ED(Nmax=10),
                             NeutralGap(), Semicircle()),
    :ll_mixing   => Pipeline(HaldaneSphere(), LLMixingCoulomb(n_ll=3),
                             ED(Nmax=10), NeutralGap(), Semicircle()),
    :cf_only     => Pipeline(Planar(), BareCoulomb(), CFTheory(),
                             CFGap(), Semicircle()),
    :dmrg_5half  => Pipeline(Cylinder(L=20), EffectiveSecondLL(),
                             DMRG(χ=8000), TransferMatrixGap(), Semicircle()),
    :chalker     => Pipeline(HaldaneSphere(), BareCoulomb(), ED(Nmax=10),
                             NeutralGap(), ChalkerCoddington(M=64)),
)
```

## 4. Physics pipeline: derivation chain

The derivation proceeds through 5 layers, each a controlled approximation. The full chain is documented in the AF proof tree (42 nodes in `proof/ledger/`). Summary below.

### 4.1 Layer 0→1: Crystal → Clean 2DEG

**Approximation A1 (effective mass):** Born-Oppenheimer + Bloch theorem → E(k) = ℏ²k²/(2m*), m* = 0.067 mₑ. Error: < 10⁻⁶.

**Approximation A2 (2D confinement):** Quantum well width ~10 nm, subband spacing ~10³ K. Project to lowest subband. Error: modifies pseudopotentials by 10–20% (Peterson et al.).

**Approximation A3 (drop disorder, phonons, spin-orbit):** Disorder → phenomenological via mobility. Phonons irrelevant at mK. Zeeman → full spin polarisation at FQHE fields. g* = −0.44.

**Result E1:** H_2DEG = Σᵢ (1/2m*)(pᵢ + eA(rᵢ))² + (1/2) Σ_{i≠j} e²/(ε|rᵢ − rⱼ|).

### 4.2 Layer 1→2: 2DEG → LLL Hamiltonian

**Exact (Landau quantisation):** εₙ = (n + 1/2)ℏωc. Filling ν = nₑh/(eB). LLL wavefunctions φₘ(z) = zᵐ exp(−|z|²/4).

**Key result (two-body problem in LLL):** Analyticity of LLL wavefunctions constrains the two-body state at fixed relative angular momentum m to be unique. Energy eigenvalue vₘ (the m-th Haldane pseudopotential) is independent of centre-of-mass momentum. The interaction is completely specified by {v₁, v₃, v₅, …}.

**Approximation A4 (LLL projection):** Project out n > 0 Landau levels. At B = 10 T: ℏωc ~ 200 K, e²/(εℓ) ~ 160 K. Ratio ~ 1.3. This is the worst approximation in the chain. LL mixing shifts gaps by 10–20% (Faugno et al. 2023).

**Result:** H_LLL = Π₀ V_ee Π₀. Kinetic energy is a constant, dropped. Entire physics determined by Coulomb interaction alone.

### 4.3 Layer 2→3: LLL → Haldane sphere + pseudopotentials

**Approximation A5 (plane → sphere):** Sphere S² of radius R = √S ℓ, threaded by 2S flux quanta. LLL has 2S + 1 orbitals (finite Hilbert space). Full SO(3) symmetry. Shift: 2S = N/ν − S_shift. Exact in N → ∞ limit.

**Exact (pseudopotential decomposition on sphere):**

Π_S V Π_S = Σ_J V_J^(S) P_J

Two-body matrix elements factorise via CG coefficients:

⟨m₁,m₂|V|m₃,m₄⟩ = δ_{m₁+m₂,m₃+m₄} Σ_J ⟨S,m₁;S,m₂|J,M⟩ V_J^(S) ⟨J,M|S,m₃;S,m₄⟩ / R

**Exact (closed-form Coulomb pseudopotentials, Fano et al. eq. 25):**

V_J^(S) = 2 C(4S−2J, 2S−J) C(4S+2J+2, 2S+J+1) / C(4S+2, 2S+1)²

in units of 1/R. Fermions: only 2S − J odd contribute.

NOTE: The current `pseudopotentials.jl` does NOT implement this formula. It has a CG-coefficient sum with a missing l-dependent Coulomb coefficient. This must be fixed before any downstream computation.

### 4.4 Layer 3→4: Pseudopotentials → ED → gaps

**Fock space enumeration:** C(2S+1, N) states as UInt64 bitstrings, restricted to Lz = 0 sector. Dimensions: ν = 1/3, N = 8 → ~9000; N = 10 → ~100,000.

**Hamiltonian assembly:** Sparse matrix from CG coefficients + pseudopotentials. Real symmetric, CSC format. Fermion signs from anticommutation (count occupied orbitals hopped over).

**Diagonalisation:** KrylovKit.eigsolve, lowest 3–5 eigenvalues. Ground state has L = 0. Neutral gap: Δ = E₁(L > 0) − E₀(L = 0).

**Approximation A6 (finite-size extrapolation):** Δ(N) = Δ(∞) + a/N from N = 6, 8, 10. Moderate accuracy (~5%).

**Physical gap conversion:** Δ_phys = Δ_dimless × e²/(4πε₀ε_r ℓ_B). Numerically: E_C/k_B ≈ 50.8 √(B[T]) K for GaAs.

### 4.5 Layer 4→5: Gaps → transport

**Exact (plateau values):** R_xy = h/(νe²). Topological, from gauge invariance + gap existence.

**Semi-phenomenological (activated transport):** R_xx ∝ exp(−Δ/(2k_BT)). Boltzmann weight for thermally excited quasiparticle pairs. Gap Δ is ab initio. Prefactor is phenomenological.

**Phenomenological (plateau width):** Width in B depends on Γ/Δ where Γ = ℏe/(m*μₑ). Requires mobility μₑ as input. Sample-dependent.

**Exact (Dykhne–Ruzin semicircle law):** Transition shape between plateaus is an exact duality result for a two-phase system. (σ*_xy − σ⁰_xy)² + (σ*_xx)² = (σ⁰_xx)².

### 4.6 Approximation severity ranking

| Rank | ID | Approximation | Error |
|---|---|---|---|
| 1 | A4 | LLL projection | Ratio ℏωc/E_C ~ 1.3; LL mixing 10–20% |
| 2 | A6 | Finite-size extrapolation | ~5% from N = 10 |
| 3 | A5 | Sphere geometry | ~1/N artefacts |
| 4 | A2 | 2D confinement | Finite width 10–20% on pseudopotentials |
| 5 | A3 | No disorder | Affects widths only |
| 6 | A1 | Effective mass | ~10⁻⁶ |

### 4.7 Honesty boundary

**Ab initio (zero free parameters):** Plateau positions, gap magnitudes, gap hierarchy, temperature dependence form.

**Phenomenological (requires disorder input):** Plateau widths, transition shapes, R_xx scale.

Everything vertical on the plot (plateau heights, R_xx minima depths) is ab initio. Everything horizontal (widths in B) is phenomenological.

## 5. Known bugs and issues in current codebase

### 5.1 Critical: `pseudopotentials.jl` formula is wrong

The function `coulomb_pseudopotentials` computes a sum over CG coefficients with an l-independent prefactor `1/(2√(S + 1/2))`. The correct formula is a closed-form expression in binomial coefficients (Fano et al. eq. 25, node 1.3.4). No CG sum needed. As written, the pseudopotentials will be wrong, and everything downstream depends on them.

**Fix:** Replace the function body with the binomial formula from node 1.3.4.

### 5.2 Critical: `sphere.jl` integer division fragility

`sphere_flux` uses `N ÷ numerator(ν) * denominator(ν) - S_shift`, which relies on integer division working out by accident. For some (N, ν) combinations it will silently truncate.

**Fix:** Use `twoS = Int(N // ν) - S_shift` with rational arithmetic and assert integrality.

### 5.3 Moderate: `monopole_harmonics.jl` / `hilbert_space.jl` convention clash

`lll_orbitals` returns half-integer m-values (stored as integers 2m), but `hilbert_space.jl` maps orbital indices to m via `m_i = index_i - S` with `S = twoS ÷ 2`. For odd twoS, the `÷ 2` truncation silently breaks Lz counting.

**Fix:** Decide on one convention. Recommend: always work with integer orbital indices 0…2S, compute Lz as Σ(index − S), and handle the half-integer case when twoS is odd by working in units of ℏ/2 throughout.

### 5.4 Minor: `shift` function incomplete

Does not handle particle-hole conjugate fractions (e.g. ν = 2/3) or fractions outside the principal Jain sequences.

**Fix:** Add particle-hole conjugate logic: shift(1 − ν) relates to shift(ν).

## 6. Implementation plan

### Phase C1: Fix foundations (1–2 hours)

1. Fix `pseudopotentials.jl` — implement Fano eq. 25 (binomial formula)
2. Fix `sphere.jl` — rational arithmetic for sphere_flux
3. Reconcile angular momentum conventions
4. Fix `shift` for conjugate fractions
5. Write tests: compare V_J against Fano Table values to 6 sig. fig.

### Phase C2: Hamiltonian + ED (2–3 hours)

1. Implement `hamiltonian.jl` — sparse assembly from CG + pseudopotentials
2. Implement `exact_diag.jl` — KrylovKit wrapper
3. Validate at N = 3 (ν = 1/3): ground state non-degenerate, L = 0
4. Validate at N = 6: gap ≈ 0.1 e²/(εℓ)
5. Run for N = 6, 8, 10 at ν = 1/3, 2/5, 3/7, 2/3, 1/5
6. Finite-size extrapolation → GapTable

### Phase C3: Transport + plot (1–2 hours)

1. Implement `transport.jl` — semicircle law + Arrhenius
2. Implement `composite_fermion.jl` — CF gap estimates for cross-check
3. Generate two-panel plot with CairoMakie
4. Annotate fractions, state what is ab initio vs phenomenological

### Phase C4: Validation data (1 hour)

1. Hard-code reference gap table from Fano et al. + Balram et al.
2. Digitise Störmer Figure 18 (WebPlotDigitizer → CSV)
3. Compile activation gap table from Du et al. / Pan et al.
4. Write comparison tests and overlay plot

### Phase C5: Release as benchmark (optional)

1. Package digitised experimental data + computed gaps as `sources/data/`
2. Document provenance of each number
3. This becomes the `fqhe-benchmarks` community resource

## 7. Future extensions

### 7.1 ν = 5/2 via tensor networks

- Cylinder geometry, iDMRG via ITensors.jl
- Effective second-LL pseudopotentials with LL mixing (Sodemann-MacDonald)
- Gap from transfer matrix spectrum
- Interface: produces a GapTable entry, transport model unchanged

### 7.2 Chalker-Coddington network model

- Microscopic disorder → localisation length → plateau width
- Replaces phenomenological Γ = ℏe/(m*μₑ) with computed localisation transition
- ~200 lines of Julia: random unitary network, transfer matrix, Lyapunov exponents
- Side quest, self-contained, computationally cheap

### 7.3 Finite well width corrections

- Peterson et al. (P13): modify pseudopotentials using quantum well wavefunction
- Inserts between stages 1 and 2
- Reduces systematic error from A2 (~10–20% → ~5%)

### 7.4 LL mixing corrections

- Faugno et al. (P20): perturbative LL mixing shifts gaps 10–20%
- Modifies effective pseudopotentials
- Same insertion point as 7.3

## 8. Experimental data strategy

### 8.1 Current state of available data

No open-access raw FQHE transport data exists in any repository. The iconic experimental traces (Tsui 1982, Störmer 1999) predate open-data mandates. The community is small enough that everyone either runs their own ED or knows the published gap values. There is no canonical `fqhe-benchmarks` repo.

### 8.2 What to compile

**Computational reference (type from papers):**

| ν | N = 6 | N = 8 | N = 10 | Δ(∞) | Source |
|---|---|---|---|---|---|
| 1/3 | | | | 0.1036 | Fano 1986 |
| 2/5 | | | | ~0.06 | Balram 2018 |
| 3/7 | | | | ~0.03 | Balram 2018 |
| 1/5 | | | | 0.0244 | Fano 1986 |

Fill N-dependent values from our own ED runs.

**Experimental activation gaps (type from papers):**

Sources: Du et al. PRL 70 (1993), Pan et al. PRL 83 (1999), Boebinger et al. PRB 36 (1987). These give Δ in kelvin at specific B-fields — directly comparable to our ED output without any transport model.

**Digitised transport traces:**

Störmer Figure 18 (Rev. Mod. Phys. 71, 875, 1999) via WebPlotDigitizer. Two curves (R_xy and R_xx vs B), ~500–1000 points each. Save as CSV in `sources/data/`.

### 8.3 Validation strategy

Two independent tests:

1. **Gap test (ab initio only):** Compare computed Δ(ν) against published ED values and experimental activation gaps. No transport model involved. This validates stages 1–4.

2. **Plot test (full pipeline):** Overlay computed R(B) curves on digitised Störmer data. This validates the whole pipeline including phenomenological transport. Visual comparison — the "does it look right" test.

## 9. Constraints

- 100% Julia. No Python, no Fortran, no wrappers.
- Julia 1.10+. All deps pinned in Project.toml.
- No single source file > 250 lines.
- Every physical constant and formula cites a source in `./sources/`.
- KrylovKit for Lanczos (pure Julia, well-tested). Do not write your own.
- Do not exceed N = 12 for ED. N = 10 is sufficient.
- Do not attempt ν = 5/2 in MVP. Separate extension.
- Do not model disorder microscopically in MVP. Phenomenological via mobility.
- Use Float64 throughout. Energy differences that give the gap are small.

## 10. Dependencies

| Package | Purpose | Notes |
|---|---|---|
| SparseArrays (stdlib) | Hamiltonian storage | |
| LinearAlgebra (stdlib) | Matrix operations | |
| KrylovKit.jl | Lanczos eigensolver | Pure Julia, preferred over Arpack |
| WignerSymbols.jl | CG coefficients for matrix elements | |
| Combinatorics.jl | Fock state enumeration | |
| CairoMakie.jl | Publication-quality figures | |

Future (post-MVP): ITensors.jl (tensor networks), Unitful.jl (dimensional analysis), Pluto.jl (interactive notebook).

## 11. Risks

| Risk | Severity | Mitigation |
|---|---|---|
| Pseudopotential formula wrong | Critical | Validate against Fano tables before any ED |
| Finite-size gaps for higher fractions poorly converged | High | Cross-check with CF theory; accept 20% uncertainty |
| Transport model tuning dominates visual quality | Medium | Accept phenomenological widths; focus validation on gap values |
| LL mixing makes absolute gap comparison with experiment poor | Medium | Note systematic uncertainty; implement corrections in extension |
| Angular momentum convention bugs | Medium | Resolve conventions in Phase C1 before proceeding |
| N = 10 ED for some fractions too slow | Low | Expected < 10 min per fraction; profile first |
