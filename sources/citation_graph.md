# Citation Graph Analysis — FQHE Project

Date: 2026-03-23
Method: Semantic Scholar API, forward citations of P1 (Haldane), P2 (Laughlin), P3 (Jain), filtered post-2020.

## Summary

| Seed paper | Post-2020 citing papers |
|------------|----------------------|
| P1 Haldane PRL 51 605 (1983) | 298 |
| P2 Laughlin PRL 50 1395 (1983) | 858 |
| P3 Jain PRL 63 199 (1989) | 312 |

The post-2020 landscape is dominated by **fractional Chern insulator** work
(twisted MoTe₂, magic-angle graphene) — a new experimental platform for FQHE
physics at zero magnetic field. These are not directly relevant to our GaAs
Haldane-sphere calculation but confirm that the core ED+pseudopotential
methodology remains the standard computational approach.

## Papers added to sources/papers/ (P16–P24)

### Directly useful for ED gaps and validation

| ID | Reference | Why we need it |
|----|-----------|----------------|
| P16 | Balram et al., "Very-High-Energy Collective States of Partons in FQH Liquids", PRX 12, 021008 (2022) [arXiv:2111.10395] | ED on sphere, excitation spectra beyond magnetoroton, modern benchmark gaps |
| P17 | Balram et al., "Multiple Magnetorotons and Spectral Sum Rules in FQH Systems", PRL 128, 246402 (2022) [arXiv:2111.10593] | Magnetoroton gap values at ν=2/7, 2/9; spectral sum rules for validation |
| P18 | Feldman & Halperin, "Fractional charge and fractional statistics in the QHE", Rep. Prog. Phys. 84, 076501 (2021) [arXiv:2102.08998] | Comprehensive 2021 review; modern gap values, transport, interferometry. 100 citations. |
| P19 | Goldman & Fradkin, "Dirac composite fermion theory of general Jain sequences", PRR 3, 033217 (2021) [arXiv:2105.02092] | CF theory for ν=N/(4N±1), Haldane bound, structure factor constraints |
| P20 | Faugno et al., "CF Pairing Induced by Landau Level Mixing", PRL 130, 186302 (2023) [arXiv:2211.07335] | LL mixing corrections to gaps — important for realistic GaAs comparison |
| P21 | Faugno et al., "CF pairing at half-filled and quarter-filled LLL", PRB 109, 035306 (2024) [arXiv:2311.05083] | Modern ED gap calculations, pairing physics at ν=1/2, 1/4 |

### Useful context and cross-checks

| ID | Reference | Why we need it |
|----|-----------|----------------|
| P22 | Kumar et al., "FQHE at ν=2+4/9", PRR 2, 032035 (2020) [arXiv:2003.07038] | Experimental gap for higher-order fraction, GaAs data |
| P23 | Lian & Zhang, "Dual Haldane sphere and quantized band geometry", PRB 103, L081103 (2021) [arXiv:2004.03609] | Modern treatment of Haldane sphere geometry |
| P24 | Ortiz et al., "Low-complexity eigenstates of a ν=1/3 FQH system", J.Phys.A 53, 505302 (2020) [arXiv:2006.00300] | Analytical structure of ν=1/3 eigenstates, validation approach |

## High-impact post-2020 papers NOT downloaded (outside scope)

These are highly cited but target different systems (not GaAs 2DEG):

- Cai et al., "Signatures of fractional quantum anomalous Hall states in twisted MoTe₂" Nature (2023), 689 cit — **zero-field FCI**
- Zeng et al., "Observation of fractionally quantized anomalous Hall effect" Nature (2023), 648 cit — **zero-field FCI**
- Xie et al., "Fractional Chern insulators in magic-angle TBG" Nature (2021), 462 cit — **FCI in graphene**
- Nakamura et al., "Direct observation of anyonic braiding statistics" Nature Physics (2020), 431 cit — **interferometry**
- Bartolomei et al., "Fractional statistics in anyon collisions" Science (2020), 380 cit — **anyon collision experiment**

These confirm the vitality of FQHE physics but don't provide gap values or
transport models for our specific GaAs Haldane-sphere calculation.

## Key findings for implementation

1. **No fundamentally new gap values**: The Coulomb gaps at standard Jain fractions
   on the Haldane sphere have not been significantly updated since P4 (Fano 1986),
   P14a/P14b (Balram 2018), and P13 (Peterson 2008). Our ED will reproduce these.

2. **LL mixing matters for quantitative GaAs comparison**: P20 (Faugno 2023) shows
   that Landau level mixing shifts gaps by ~10-20% at typical GaAs fields. We should
   note this as a systematic uncertainty in our transport model.

3. **Magnetoroton structure**: P17 shows multiple magnetoroton branches exist for
   higher-order fractions. For ν=1/3 the single-magnetoroton picture is adequate.

4. **Transport model**: No major new transport models found post-2020. The semicircle
   law (P10, Dykhne-Ruzin 1994) and activated transport remain standard. P18
   (Feldman-Halperin 2021 review) confirms this.

5. **CF effective mass**: Modern values for m_CF consistent with Jain's estimates.
   P19 (Goldman-Fradkin) provides theoretical constraints.

## Paper count

| Category | Count |
|----------|-------|
| Original plan (P1–P15) | 15 |
| Downloaded (P1–P15, excl. textbooks P8, P9) | 13 |
| Added from citation graph (P16–P24) | 9 |
| **Total downloaded** | **24** |
| Textbooks not obtained (P8, P9) | 2 |

Well within the ~40 paper cap. The Girvin lectures (P11, 124 pages) serve as
textbook substitute per the original plan.
