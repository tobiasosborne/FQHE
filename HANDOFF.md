# HANDOFF — FQHE Project

**Date:** 2026-03-23
**Status:** Phase A (Source Acquisition) COMPLETE. Phase B (equation graph) COMPLETE. Phase C (implementation) NOT STARTED.

## What was done this session

### 1. Paper acquisition (24 papers, all local)

**All 15 papers from the plan (P1–P15) were acquired except textbooks P8/P9:**

| ID | Source | How obtained |
|----|--------|-------------|
| P01–P04, P06, P07, P10, P12, P15 | APS journals (PRL, PRB, RMP) | **Playwright headed browser via TIB VPN** |
| P05a, P05b, P11, P13, P14a, P14b | arXiv | Direct curl |
| P16–P24 | arXiv (post-2020 from citation graph) | Direct curl |
| P08 (Jain textbook), P09 (Chakraborty) | Not obtained | Textbooks — P11 (Girvin, 124pp) substitutes |

**Papers are in `sources/papers/` but gitignored (copyrighted PDFs).**

### 2. How to re-obtain the papers on a new machine

**arXiv papers (no paywall, just curl):**
```bash
# Run from project root:
curl -sL -o sources/papers/P05a_Morf_cond-mat_0202407.pdf https://arxiv.org/pdf/cond-mat/0202407
curl -sL -o sources/papers/P05b_dAmbrumenil_1008.0969.pdf https://arxiv.org/pdf/1008.0969
curl -sL -o sources/papers/P11_Girvin_cond-mat_9907002.pdf https://arxiv.org/pdf/cond-mat/9907002
curl -sL -o sources/papers/P13_Peterson_0801.4819.pdf https://arxiv.org/pdf/0801.4819
curl -sL -o sources/papers/P14a_Balram_1803.10427.pdf https://arxiv.org/pdf/1803.10427
curl -sL -o sources/papers/P14b_Balram_1807.02997.pdf https://arxiv.org/pdf/1807.02997
curl -sL -o sources/papers/P16_Balram_VeryHighEnergy_2111.10395.pdf https://arxiv.org/pdf/2111.10395
curl -sL -o sources/papers/P17_Balram_Magnetorotons_2111.10593.pdf https://arxiv.org/pdf/2111.10593
curl -sL -o sources/papers/P18_Feldman_FracChargeStats_2102.08998.pdf https://arxiv.org/pdf/2102.08998
curl -sL -o sources/papers/P19_Goldman_DiracCF_Jain_2105.02092.pdf https://arxiv.org/pdf/2105.02092
curl -sL -o sources/papers/P20_Faugno_CFPairing_LLMixing_2211.07335.pdf https://arxiv.org/pdf/2211.07335
curl -sL -o sources/papers/P21_Faugno_CFPairing_quarter_2311.05083.pdf https://arxiv.org/pdf/2311.05083
curl -sL -o sources/papers/P22_Kumar_FQHE_2plus4over9_2003.07038.pdf https://arxiv.org/pdf/2003.07038
curl -sL -o sources/papers/P23_Lian_DualHaldane_2004.03609.pdf https://arxiv.org/pdf/2004.03609
curl -sL -o sources/papers/P24_Ortiz_LowComplexity_nu13_2006.00300.pdf https://arxiv.org/pdf/2006.00300
```

**APS papers (paywalled, need TIB VPN + headed browser):**

APS uses Cloudflare bot protection. `curl` gets 403 even with cookies — Cloudflare
fingerprints the TLS stack. The method that works:

1. **Connect to TIB VPN** (must be active — IP authentication)
2. **Open headed Playwright browser** (headless gets blocked):
   ```bash
   playwright-cli open --browser=chromium --headed --persistent
   ```
3. **Navigate to any APS abstract** — a Cloudflare challenge page appears.
   **Manually click the challenge** in the browser window (WSLg display required).
4. **Once past Cloudflare**, use `page.request.get()` to fetch PDFs through the
   authenticated browser session. The script `scripts/fetch_via_browser.sh` automates
   this — it base64-encodes the PDF from `page.request.get()` and pipes to disk:
   ```bash
   # After browser is open and Cloudflare is passed:
   bash scripts/fetch_via_browser.sh
   ```
5. The standalone Node.js script `scripts/fetch_aps_papers.mjs` does NOT work — a
   fresh headless browser gets its own Cloudflare challenge. Must use the existing
   browser session.

**APS paper URLs (DOI → PDF):**
| Paper | PDF URL |
|-------|---------|
| P01 Haldane | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.51.605` |
| P02 Laughlin | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.50.1395` |
| P03 Jain | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.63.199` |
| P04 Fano | `https://journals.aps.org/prb/pdf/10.1103/PhysRevB.34.2670` |
| P06 Tsui | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.48.1559` |
| P07 Willett | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.59.1776` |
| P10 Dykhne | `https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.2369` |
| P12 Störmer | `https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.71.875` |
| P15 Haldane-Rezayi | `https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.54.237` |

### 3. Citation graph analysis

Forward citations of P1, P2, P3 via Semantic Scholar API. Results in
`sources/citation_graph.md`. 9 post-2020 papers added (P16–P24).
Key finding: no new gap benchmarks since Balram (2018); LL mixing corrections
~10-20% (Faugno 2023); transport model unchanged.

### 4. Equation derivation graph (AF proof tree)

Complete derivation chain from full many-body Hamiltonian to Julia equations,
built as an AF (Adversarial Proof Framework) proof tree in `proof/`.

```bash
cd proof && af status   # View the full 31-node tree
af get 1.3.4            # Example: see the pseudopotential formula node
```

**5 layers, 31 nodes**, every equation tagged with paper reference:
- Layer 0→1: Full H → clean 2DEG (Born-Oppenheimer, eff mass, 2D confinement)
- Layer 1→2: 2DEG → LLL Hamiltonian (Landau levels, LLL projection)
- Layer 2→3: LLL → Haldane sphere + pseudopotentials (P01/P04 equations)
- Layer 3→4: Pseudopotentials → sparse ED → gaps (Fock space, Lanczos, extrapolation)
- Layer 4→5: Gaps → R_xy, R_xx (topological plateaus, activated transport, semicircle law)

6 approximations ranked by severity:
1. A4 (LLL projection) — worst, ratio ~1.3
2. A6 (finite-size extrapolation) — moderate
3. A5 (sphere geometry) — exact in limit
4. A2 (2D confinement) — ~10-20% on pseudopotentials
5. A3 (no disorder) — widths only
6. A1 (effective mass) — negligible

## What to do next (Phase C: Implementation)

Follow the order in CLAUDE_fqhe.md Sec C1. All equations are now grounded in
local papers and organized in the AF proof tree. Start with:

1. `materials.jl` — GaAs constants (trivial, ~30 lines)
2. `landau.jl` — ν(B), ℓ_B(B), energy scales
3. `pseudopotentials.jl` — **V_J^(S) from P04-Fano eq (25)** ← proof node 1.3.4
4. `hilbert_space.jl` — Fock space enumeration ← proof node 1.4.1
5. `hamiltonian.jl` — sparse H via CG coefficients ← proof node 1.4.2
6. `exact_diag.jl` — KrylovKit wrapper ← proof node 1.4.3
7. Continue per CLAUDE_fqhe.md...

## Prerequisites on new machine

- Julia 1.10+
- `af` CLI (go binary, `go install` or copy from `~/go/bin/af`)
- `playwright-cli` (npm: `@playwright/cli`) — only needed for APS paper re-fetch
- TIB VPN — only needed for APS paper re-fetch
- WSLg or X11 display — only needed for headed browser
