# Tensor Network Methods for the Fractional Quantum Hall Effect: A Survey of Recent Progress (2020--2026)

**Compiled: 2026-03-24**

---

## 1. iDMRG and Infinite DMRG for FQHE

### State of the Art

The infinite density matrix renormalization group (iDMRG) is the leading tensor network
method for studying FQHE on an infinite cylinder. The method works in the Landau gauge,
mapping the 2D problem onto a 1D chain of orbitals, and exploits the matrix product state
(MPS) ansatz to find the ground state variationally.

**Key capabilities of iDMRG for FQHE:**

- Finds the full set of topologically degenerate ground states on an infinite cylinder.
- Extracts topological invariants directly: quantum dimensions, topological spins,
  quasiparticle charges, chiral central charge, and Hall viscosity.
- Handles long-range Coulomb interactions via matrix product operator (MPO) compression.

**Bond dimensions and system sizes (state of the art):**

| Study | Bond dim (chi) | Circumference L | Filling | Notes |
|-------|---------------|-----------------|---------|-------|
| Zaletel, Mong, Pollmann (2015) | up to 15,000 | L = 17--25 l_B | nu=5/2, 7/3 | 3 Landau levels included |
| Hu, Zhu, Zhu (2012) | up to 5,000 | N/A (sphere) | nu=1/3, 5/2 | Up to 20e (1/3), 30e (5/2) |
| Feiguin et al. (2007) | up to 5,000 | N/A (sphere) | nu=1/3, 5/2 | Up to 20e (1/3), 26e (5/2) |
| Recent FCI studies (2024--2025) | up to 16,000 | YC5, YC7 | Various | Hofstadter-Hubbard models |
| Twisted MoTe2 DMRG (2025) | up to 4,096 | XC4--XC6 | nu=1/3, 2/3 | Fractional Chern insulators |

For multicomponent systems (spin + valley), iDMRG memory scales polynomially (as N^3 in the
number of components), versus exponentially for ED. This allowed the first numerical resolution
of the Pfaffian vs anti-Pfaffian debate at nu=5/2 with Landau level mixing: at kappa=1.38,
keeping three Landau levels, the ground state is anti-Pfaffian (Zaletel et al. 2015).

### Geometry

iDMRG on an infinite cylinder is the dominant geometry. The cylinder circumference L (in
units of magnetic length l_B) is the key parameter controlling accuracy. Finite cylinders
(open boundary in x) are also used but are less common for iDMRG.

The mixed real-momentum space approach (Motruk, Zaletel, Mong, Pollmann 2016) uses real
space along the cylinder and momentum space around the circumference, providing significant
speedup and memory reduction.

---

## 2. MPS Representations of FQH States

### Exact MPS from Conformal Field Theory

A landmark result (Estienne, Bernevig, Regnault 2012--2013) showed that model FQH
wavefunctions can be written as **exact** matrix product states when the wavefunction
is expressible as a correlator in a 1+1D conformal field theory (CFT). This applies to:

- **Laughlin states** (all 1/q fillings)
- **Moore-Read (Pfaffian) state** and all paired k=2 Jack polynomial states
- **Gaffnian state**
- **Read-Rezayi k=3 state**
- **Halperin multicomponent states** (exact MPS derived 2018)

The auxiliary Hilbert space of the MPS is precisely the CFT Hilbert space, and the
MPS bond dimension (after truncation) controls the accuracy. The entanglement spectrum
of these MPS exactly reproduces the edge CFT counting.

### Accuracy Achieved

For model Hamiltonians (e.g., V_1 pseudopotential for Laughlin at nu=1/3):
- The overlap between DMRG ground state and Laughlin wavefunction is 0.992--0.994 for
  systems up to 20 electrons.
- Ground state energy accuracy: errors of order 10^{-4} e^2/l_B, which is 2--3 orders
  of magnitude below the excitation gap.

For Coulomb Hamiltonian at nu=1/3:
- Thermodynamic limit energy: E_0 = -0.410164(4) e^2/l_B per electron (iDMRG benchmark).
- Earlier DMRG on sphere: E_0 = -0.41010 e^2/l_B per electron.
- Laughlin variational energy: E_0 = -0.40977 e^2/l_B per electron.

The exact MPS representations are particularly powerful for computing entanglement spectra,
correlation lengths, and topological entanglement entropy with arbitrary precision (limited
only by the CFT truncation level).

### Entanglement Spectrum

The real-space entanglement spectrum obeys a scaling form dictated by the edge CFT,
allowing extraction of the two entanglement velocities of the Moore-Read state. The low-lying
entanglement spectrum serves as a "fingerprint" for topological order identification.

---

## 3. DMRG on Cylinders

### The Cylinder Advantage

The cylinder geometry is the natural geometry for DMRG applied to FQHE because:
1. Translation invariance both along and around the cylinder.
2. Zero curvature effects (unlike the sphere).
3. Absence of gapless edge excitations in the bulk (unlike open geometries).
4. Remarkably rapid convergence compared to sphere geometry.

In the Landau gauge, periodic boundary conditions in y map the 2D electron problem to a
1D chain of Landau-level orbitals. The entanglement of a bipartition across the cylinder
satisfies an area law: S ~ alpha * L, where L is the cylinder circumference. This means
the required bond dimension grows exponentially with L: chi ~ exp(alpha * L).

### Reachable Circumferences

The maximum achievable circumference depends on the filling fraction and the topological
state:

| Filling | Typical L_max (l_B) | Bond dim required | Notes |
|---------|--------------------|--------------------|-------|
| nu=1/3 (Laughlin) | ~20--25 l_B | chi ~ 5,000--15,000 | Area-law prefactor is small |
| nu=5/2 (Moore-Read) | ~20 l_B | chi ~ 6,300--15,000 | Non-Abelian state, larger entanglement |
| nu=7/3 | ~24 l_B | chi ~ 15,000 | 3.5 days, 60 GB on 16-core node |
| FCI (lattice) | YC5--YC7 (lattice units) | chi ~ 6,000--16,000 | Hofstadter/moire models |

For the Laughlin nu=1/3 state, the entanglement is relatively modest: even chi=1000--2000
gives excellent convergence on a cylinder. For more exotic states (Moore-Read, Read-Rezayi),
the entanglement grows faster with L.

### Best Energy-per-Particle Results

**nu = 1/3 Coulomb (thermodynamic limit):**
- E_0 = -0.410164(4) e^2/epsilon*l_B per electron (iDMRG, Zaletel et al. 2015)
- E_0 = -0.41010 e^2/epsilon*l_B per electron (finite DMRG on sphere, Feiguin et al. 2007)
- These are consistent within error bars and represent the best numerical estimates.

**nu = 5/2 Coulomb (thermodynamic limit):**
- E_0 ~ -0.362 to -0.366 e^2/epsilon*l_B per electron (depending on LL mixing treatment)
- These values have been confirmed by multiple groups.

**Excitation gaps (thermodynamic limit extrapolations):**
- nu=1/3 neutral gap: Delta ~ 0.1010 e^2/epsilon*l_B (DMRG), 0.1012 (prior estimates)
- nu=5/2 gap: Delta ~ 0.025--0.030 e^2/epsilon*l_B (upper bound from DMRG)

**Experimental validation:** Yang, Zibrov, Zaletel et al. (PRL 2021) measured energy per
particle in partially filled Landau levels in graphene and found "exceptional agreement with
numerical calculations over the whole Landau level without adjustable parameters" when
screening by filled Landau levels is included.

---

## 4. Tensor Network Approaches Beyond MPS

### PEPS for FQH States

Projected Entangled Pair States (PEPS) are the natural 2D generalization of MPS. Progress
on PEPS for chiral topological states has been historically slow due to fundamental
difficulties: chiral PEPS necessarily have algebraically decaying "gossamer" correlation
tails, and the sign problem in contraction makes optimization challenging.

**Recent breakthroughs (2023--2025):**

- **Variational iPEPS for Laughlin states** (2023--2024): Variational optimization of
  infinite PEPS successfully identified fractional Hall states in the bosonic
  Harper-Hofstadter model, with exponentially decaying bulk correlations and chiral edge
  modes visible in the entanglement spectrum. Published in Phys. Rev. B 109, L241117 (2024).

- **Bulk gap simulation in chiral PEPS** (2025): A new framework demonstrates that PEPS can
  reliably simulate gapped bulk excitations in chiral topological phases, resolving concerns
  about whether gossamer tails make the bulk falsely gapless. For Z_2 chiral states, PEPS
  achieves ground state energies matching exact results to 4 significant figures
  (E_0 = -1.6979J vs exact -1.6980J). For Z_3 states, bulk gap Δ = 0.71J was obtained.

- **Fermionic isometric tensor network states (isoTNS)** (Dai, Wu, Wang, Zaletel 2022,
  published PRL 2025): A new class of 2D tensor networks with isometric constraints that
  allow efficient contraction. Algorithms for both bosonic and fermionic systems have been
  developed and are being evaluated for characterizing fractional quantum Hall states.

- **Alternating isoTNS** (2025): An improved variant where isometric arrows alternate
  direction on columns, providing significantly better representation of 2D ground states
  than the original isoTNS construction.

### MERA

The Multi-scale Entanglement Renormalization Ansatz (MERA) has been less directly applied to
FQHE in continuum. MERA can extract central charges, scaling dimensions, and OPE coefficients
for critical systems, but its application to gapped 2D topological phases in practice has
been limited compared to DMRG. No major FQHE-specific MERA results have appeared in the
2020--2026 period.

---

## 5. Finite-Size Extrapolation from Tensor Network Data

### Methods for Thermodynamic Limit Extrapolation

**From finite DMRG (sphere or cylinder):**
- Ground state energies E(N) are computed for a sequence of particle numbers N.
- Extrapolation to N -> infinity uses fits to E(N)/N = e_inf + a/N + b/N^2 + ...
- Magnetic length rescaling (replacing l_B with l_B' to account for curvature on sphere)
  improves convergence.
- Data from both ED and DMRG at overlapping system sizes validate consistency.

**From iDMRG (infinite cylinder):**
- The system is already infinite in one direction; extrapolation is in cylinder
  circumference L.
- Physical quantities Q(L) are extrapolated: Q(L) = Q_inf + a * exp(-L/xi) + ...
- The topological entanglement entropy gamma is extracted from: S(L) = alpha*L - gamma,
  fitting the entanglement entropy vs circumference.
- Accuracy of ~10^{-3} for gamma requires L ~ 10 * xi (ten times the correlation length).

**Bond dimension extrapolation:**
- For any finite chi, there is a systematic error. Results for chi = 50, 100, ..., 500, ...
  are extrapolated to chi -> infinity.
- The truncation error (sum of discarded singular values) provides an error estimate.

**Challenges:**
- The entanglement length scale governing convergence of topological properties is not
  necessarily the physical correlation length. Pathological examples exist where the
  entanglement length diverges while the physical correlation length stays finite.
- For gap extrapolation, accessible system sizes are often too limited for reliable
  extrapolation, especially at nu=5/2.

---

## 6. Comparison with Exact Diagonalization

### System Size Reach

| Method | nu=1/3 max N_e | nu=5/2 max N_e | Geometry |
|--------|---------------|---------------|----------|
| ED (sphere/torus) | ~16--18 | ~18--20 | Any |
| Finite DMRG (sphere) | ~20--25 | ~26--30 | Sphere |
| Finite DMRG (cylinder) | ~20+ | ~20+ | Cylinder |
| iDMRG (infinite cylinder) | Infinite x L | Infinite x L | Cylinder (1D infinite) |

**DMRG typically handles systems 1.5--2x larger than ED** for the same filling, and
iDMRG handles infinitely long systems at fixed circumference.

### Where Tensor Networks Improve Over ED

1. **System size**: DMRG reaches N_e ~ 25--30, vs ED's N_e ~ 16--20 for typical FQHE
   fillings.
2. **Entanglement diagnostics**: iDMRG directly provides the entanglement spectrum and
   topological entanglement entropy from the MPS representation.
3. **Thermodynamic limit**: iDMRG gives direct access to the infinite-length limit;
   circumference extrapolation is the remaining finite-size effect.
4. **Multicomponent systems**: Polynomial (not exponential) scaling with number of
   internal degrees of freedom (spin, valley, Landau level index).

### Where ED Still Wins

1. **Any geometry**: ED works on sphere, torus, disk, etc. DMRG is most efficient on
   cylinders.
2. **Full spectrum**: ED provides the entire energy spectrum; DMRG targets ground state
   and a few low-lying excitations.
3. **Small system benchmarks**: For N_e < 14, ED is fast and exact -- ideal for benchmarking.
4. **Landau level mixing**: Both ED and DMRG are limited by Landau level truncation.
   Neural network methods (see below) are beginning to surpass both for strong LL mixing.

---

## 7. Key Papers and Groups (2020--2026)

### Leading Groups

| Group / PI | Institution | Focus |
|-----------|-------------|-------|
| **Michael Zaletel** | UC Berkeley | iDMRG for FQHE, isoTNS, moire DMRG, TeNPy |
| **Roger Mong** | U. Pittsburgh | iDMRG for FQHE, topological characterization |
| **Frank Pollmann** | TU Munich | iDMRG, TeNPy development, isoTNS |
| **Nicolas Regnault** | ENS Paris / Flatiron | Exact MPS for FQH, entanglement spectroscopy |
| **B. Andrei Bernevig** | Princeton | MPS/CFT for FQH, FCI theory |
| **Benoit Estienne** | Sorbonne (LPTHE) | Exact MPS from CFT, Halperin states |
| **Norbert Schuch** | U. Vienna | PEPS for topological order |
| **Didier Poilblanc** | U. Toulouse | Chiral PEPS, chiral spin liquids |
| **Johannes Hauschild** | UC Berkeley / TeNPy | TeNPy development and maintenance |

### Landmark Papers (2020--2026)

**iDMRG methodology and FQHE:**
- Zaletel, Mong, Pollmann, "Topological characterization of FQH ground states from
  microscopic Hamiltonians", PRL 110, 236801 (2013). [Foundational for iDMRG+FQHE]
- Zaletel, Mong, Pollmann, Rezayi, "Infinite DMRG for multicomponent quantum Hall
  systems", PRB 91, 205106 (2015). [chi up to 15,000; Pfaffian vs anti-Pfaffian]
- Motruk, Zaletel, Mong, Pollmann, "DMRG on a cylinder in mixed real and momentum
  space", PRB 93, 155139 (2016). [Computational speedup for cylinders]

**Exact MPS for FQH states:**
- Estienne, Regnault, Bernevig, "Fractional quantum Hall matrix product states for
  interacting conformal field theories", PRB 92, 075121 (2013/2015).
- Estienne, Regnault, Bernevig, "MPS description of Halperin states",
  PRB 97, 165136 (2018).

**DMRG for moire materials and FCIs (2020--2025):**
- Soejima, Parker, Bultinck, Hauschild, Zaletel, "Efficient simulation of moire
  materials using DMRG", PRB 102, 205111 (2020). [MPO compression for moire Hamiltonians]
- Dong, Wang, Wang, Soejima, Zaletel et al., "Anomalous Hall Crystals in Rhombohedral
  Multilayer Graphene I", PRL 133, 206503 (2024). [FQHE at zero magnetic field]
- Assouline, Wang, Zhou et al., "Energy gap of the even-denominator FQH state in bilayer
  graphene", PRL 132, 046603 (2024). [DMRG + experiment comparison]
- "Fractional Quantum Hall Wedding Cakes" (2025). [MPS/iDMRG for Hofstadter-Bose-Hubbard]

**Experimental validation:**
- Yang, Zibrov, Bai, Taniguchi, Watanabe, Zaletel, Young, "Experimental determination of
  energy per particle in partially filled Landau levels", PRL 126, 156802 (2021).
  [Exceptional agreement between DMRG and experiment]

**Tensor networks beyond MPS:**
- Dai, Wu, Wang, Zaletel, "Fermionic isometric tensor network states in 2D",
  PRL 134, 026502 (2025).
- Wu, Anand, Lin et al., "Two-dimensional isometric tensor networks on an infinite strip",
  PRB 107, 245118 (2023).
- Variational iPEPS for FQH: PRB 109, L241117 (2024). [Laughlin state in Harper-Hofstadter]
- Simulating bulk gap in chiral PEPS (2025). [Resolves gossamer tail concern]

**Neural network competitors:**
- "Solving the fractional quantum Hall problem with self-attention neural network",
  PRB 111, 205117 (2025). [FNN surpasses ED for strong LL mixing at nu=1/3]
- "Taming Landau level mixing in FQH states with deep learning",
  PRL 134, 176503 (2025). [Real-space NNs beat LL-projected ED]

**Software:**
- Hauschild, Unfried, Anand et al., "Tensor Network Python (TeNPy) version 1",
  SciPost Phys. Codebases 41 (2024).

---

## 8. Open-Source Codes

### Tensor Network Libraries with FQHE Capability

| Code | Language | URL | FQHE support |
|------|----------|-----|-------------|
| **TeNPy** | Python | [github.com/tenpy/tenpy](https://github.com/tenpy/tenpy) | General MPS/iDMRG; used extensively for FQH by Zaletel group. No built-in FQH model but framework supports custom Hamiltonians on cylinders. Version 1.0 released 2024. |
| **ITensor** | Julia (also C++) | [itensor.org](https://itensor.org/) | General tensor network library; ITensorMPS.jl for DMRG. Used for Hofstadter-Hubbard FQH studies. Supports Abelian quantum number conservation. |
| **infinite_cylinder** | Python (TeNPy-based) | [github.com/bartandrews/infinite_cylinder](https://github.com/bartandrews/infinite_cylinder) | Specifically designed for 2D DMRG on infinite cylinders. Built on TeNPy 0.5+. Includes tools for saving/loading ground states and computing observables. |
| **QuantumHallED.jl** | Julia | [github.com/mishmash/QuantumHallED.jl](https://github.com/mishmash/QuantumHallED.jl) | Exact diagonalization (not tensor network) specifically for FQH systems. Useful for benchmarking TN results. Legacy Julia v0.6 code. |

### Notes on Code Availability

- The Zaletel group's production iDMRG codes for FQHE are largely based on TeNPy but are
  not separately released as standalone packages. The key innovations (MPO compression for
  long-range Coulomb, mixed real-momentum space DMRG, multi-Landau-level handling) are
  described in papers but implemented in group-internal code built on TeNPy.

- For PEPS-based FQH calculations, there is no widely available open-source code specifically
  for FQHE. Groups typically use custom implementations of iPEPS optimization and CTMRG
  contraction.

- TeMFpy (mean-field + TeNPy) integrates mean-field wavefunctions with MPS algorithms and
  can be used as a starting point for FQH and spin liquid studies.

---

## 9. Summary Assessment: Can Tensor Networks Beat Small-N ED for Thermodynamic-Limit FQHE?

**The answer is a qualified yes, and tensor networks are already the method of choice for
thermodynamic-limit estimates of ground state energies.** Specifically:

### Ground state energies
- iDMRG on infinite cylinders provides the most accurate thermodynamic-limit estimates of
  ground state energies per particle for FQHE states. The benchmark value
  E_0 = -0.410164(4) e^2/l_B at nu=1/3 from iDMRG is more reliable than any ED
  extrapolation from N_e <= 18 systems.

### Transport gaps
- Gap extrapolation remains more difficult. The transport gap requires computing excited
  states (quasihole-quasielectron pair), and finite-circumference effects can be large.
  Current DMRG estimates for the nu=1/3 gap (Delta ~ 0.101 e^2/l_B) are consistent with
  ED but come with systematic uncertainties from circumference extrapolation.
- For nu=5/2, the gap estimate (Delta ~ 0.025--0.030 e^2/l_B) is still regarded as an
  upper bound, and the true thermodynamic value may be smaller.

### Current limitations
1. **Circumference bottleneck**: The exponential growth of entanglement with cylinder
   circumference limits iDMRG to L ~ 20--25 l_B, beyond which chi > 15,000 becomes
   impractical.
2. **Landau level mixing**: Both DMRG and ED are limited by LL truncation. Neural network
   methods (2024--2025) now surpass both for strong LL mixing at nu=1/3 with up to 12
   electrons, but have not yet been scaled to the system sizes accessible to DMRG.
3. **Gap vs energy**: Ground state energies converge much faster than excitation gaps with
   system size. Thermodynamic-limit gap estimates remain challenging for all methods.

### Outlook
- The combination of iDMRG for large-circumference cylinders + finite-circumference
  extrapolation is the current gold standard for thermodynamic FQHE properties.
- PEPS methods are reaching maturity for chiral states (2024--2025 breakthroughs) and may
  eventually provide truly 2D calculations without the cylinder bottleneck.
- Neural network quantum states are emerging as a competitive alternative, especially for
  including Landau level mixing, but have not yet demonstrated advantages over DMRG at
  large system sizes.
- Fractional Chern insulator calculations in moire materials are driving a new wave of
  DMRG applications with bond dimensions up to chi ~ 200,000 for realistic material models.

---

## Sources

- [Zaletel et al., Infinite DMRG for multicomponent quantum Hall systems (arXiv:1410.3861)](https://ar5iv.labs.arxiv.org/html/1410.3861)
- [Hu, Zhu, Zhu, Comparison of DMRG in different geometries (arXiv:1202.4697)](https://ar5iv.labs.arxiv.org/html/1202.4697)
- [Feiguin et al., DMRG Study of Incompressible FQH States (arXiv:0706.4469)](https://ar5iv.labs.arxiv.org/html/0706.4469)
- [Estienne, Regnault, Bernevig, FQH MPS for interacting CFTs (arXiv:1311.2936)](https://arxiv.org/abs/1311.2936)
- [Estienne, Regnault, Bernevig, MPS for trial QH states (arXiv:1211.3353)](https://arxiv.org/abs/1211.3353)
- [Estienne, Regnault, Bernevig, MPS for Halperin states (arXiv:1802.05278)](https://arxiv.org/abs/1802.05278)
- [Motruk et al., DMRG on cylinder in mixed real and momentum space (arXiv:1512.03318)](https://ar5iv.labs.arxiv.org/html/1512.03318)
- [Zaletel, Mong, Topological characterization of FQH (arXiv:1211.3733)](https://arxiv.org/abs/1211.3733)
- [Misguich, Jolicoeur, DMRG study of FQHE in open cylinder (arXiv:2012.03000)](https://arxiv.org/abs/2012.03000)
- [Soejima et al., Efficient simulation of moire materials using DMRG (arXiv:2009.02354)](https://arxiv.org/abs/2009.02354)
- [Hauschild et al., TeNPy version 1, SciPost Phys. Codebases 41 (2024) (arXiv:2408.02010)](https://arxiv.org/html/2408.02010v1)
- [Neural network for FQHE, PRB 111 205117 (2025) (arXiv:2412.00618)](https://arxiv.org/html/2412.00618v2)
- [Taming LL mixing with deep learning, PRL 134 176503 (2025) (arXiv:2412.14795)](https://arxiv.org/abs/2412.14795)
- [Variational iPEPS for FQH, PRB 109 L241117 (2024) (arXiv:2309.12811)](https://arxiv.org/abs/2309.12811)
- [Simulating bulk gap in chiral PEPS (2025) (arXiv:2502.20142)](https://arxiv.org/html/2502.20142)
- [Dai et al., Fermionic isoTNS in 2D, PRL 134 026502 (2025) (arXiv:2211.00043)](https://ui.adsabs.harvard.edu/abs/2025PRXQ....6d0324W/abstract)
- [Wu et al., 2D isoTNS on infinite strip, PRB 107 245118 (2023) (arXiv:2211.14337)](https://ui.adsabs.harvard.edu/abs/2023APS..MARD62008D/abstract)
- [Dong et al., Anomalous Hall Crystals I, PRL 133 206503 (2024)](https://arxiv.org/abs/2311.05568)
- [Assouline et al., Energy gap of even-denom FQH in bilayer graphene, PRL 132 046603 (2024)](https://www.researchgate.net/publication/373047224_Energy_gap_of_the_even-denominator_fractional_quantum_Hall_state_in_bilayer_graphene)
- [Yang et al., Experimental determination of energy per particle, PRL 126 156802 (2021)](https://pubmed.ncbi.nlm.nih.gov/33929240/)
- [Charge gap > neutral gap (arXiv:2410.11645)](https://arxiv.org/html/2410.11645)
- [FQH Wedding Cakes (arXiv:2510.15472)](https://arxiv.org/html/2510.15472)
- [Andrews et al., Stability of FCI with non-LL continuum limit, PRB 109 245111 (2024)](https://arxiv.org/abs/2310.05758)
- [Quantum Hall to chiral spin liquid (arXiv:2510.19907)](https://arxiv.org/html/2510.19907v1)
- [bartandrews/infinite_cylinder GitHub](https://github.com/bartandrews/infinite_cylinder)
- [mishmash/QuantumHallED.jl GitHub](https://github.com/mishmash/QuantumHallED.jl)
- [Zaletel group publications](https://mikezaletel.github.io/zaletelgroupwebsite/publications/)
- [TeNPy documentation](https://tenpy.readthedocs.io/en/latest/)
- [ITensor](https://itensor.org/)
- [The Tensor Network (resource site)](https://tensornetwork.org/reviews_resources.html)
- [Scarola, FQH at nu=5/2 review chapter](https://scarola.phys.vt.edu/content/dam/scarola_phys_vt_edu/papers/3-s2.0-B9780323908009001359-main.pdf)
- [Regnault, Entanglement Spectroscopy and QH (arXiv:1510.07670)](https://arxiv.org/abs/1510.07670)
- [FCI in twisted MoTe2 (2025) (arXiv:2504.07932)](https://arxiv.org/html/2504.07932v1)
- [Assouline et al., Entropy of strongly correlated electrons in partially filled LL (2025) (arXiv:2503.16738)](https://arxiv.org/abs/2503.16738)
