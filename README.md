# Why Care Systems Fail: Institutional Decoupling in Ageing Societies

**Nature Health** | Final revised manuscript (v13) | Under blind review

**OSF pre-registration:** https://osf.io/qn7wz  
**Licence:** MIT

---

## Repository Structure

```
ltc-decoupling/
├── R/
│   ├── 00_setup.R                  # Package loading, renv
│   ├── 00_generate_synthetic_demo.R # Synthetic demo data
│   ├── 01_data_prep.R              # Panel construction (704 obs, 32 countries)
│   ├── 02_ids_construction.R       # IDS_pure + GIS + ESI; SD clarification
│   ├── 03_twfe_main.R              # Primary TWFE; entropy decomposition
│   ├── 04_iv_estimation.R          # IV; GDP trend + electoral cycle tests
│   ├── 05_mediation.R              # CLS mediation (IKY framework)
│   ├── 06_event_study.R            # Four natural experiments
│   ├── 07_goodman_bacon.R          # Goodman-Bacon decomposition
│   ├── 08_callaway_santanna.R      # Callaway–Sant'Anna GATT
│   ├── 09_synthetic_control.R      # Germany SC; PSM balance check
│   ├── 10_regime_moderation.R      # Regime-stratified regressions
│   ├── 11_falsification.R          # Placebo tests
│   ├── 12_robustness_table.R       # Full Table 2 (12 specifications)
│   ├── 13_figures.R                # Fig 1–4 + FigS2
│   ├── 14_pathway_analysis.R       # P1/P2/P3; p0 baseline reported
│   └── 15_cfa_invariance.R         # Multi-group CFA
├── data/
│   ├── codebook_CCSD.csv           # Comparative Care System Database
│   ├── expert_panel_ratings.csv    # Delphi panel boundary ratings
│   └── synthetic_demo.rds          # Reproducible demo (no restricted data)
├── output/
│   └── figures/                    # Fig1–Fig4, FigS2
├── renv.lock                       # 36 packages, versions pinned (R 4.3.1)
└── README.md
```

---

## Reproduction (demo data)

```r
# 1. Restore package environment
renv::restore()

# 2. Generate synthetic demo dataset
source("R/00_generate_synthetic_demo.R")

# 3. Run full analysis pipeline
source("R/00_setup.R")
source("R/01_data_prep.R")
source("R/02_ids_construction.R")   # Reports pooled SD=1.00, within-country SD
source("R/03_twfe_main.R")          # β=0.149; entropy decomposition
source("R/04_iv_estimation.R")      # KP F=34.7; five-layer exclusion test
source("R/05_mediation.R")          # ACME=0.091 SD (62%)
source("R/06_event_study.R")        # Four reform sites
source("R/14_pathway_analysis.R")   # P1 p0=0.32 baseline; 7.1pp derivation
```

---

## Key Technical Notes (v13)

**IDS_pure standardisation (F3):** `IDS_pure = z-score(ESI - GIS)` over the pooled
704-observation distribution. Pooled SD = 1.00 by construction; within-country SD
(after country demeaning) ≈ 0.31 SD. Both are reported by `02_ids_construction.R`.
The 7.8-fold comparison with LTC expenditure uses pooled-SD scales for both variables.

**P1 pathway absolute risk (F4):** Cumulative OR over 7.4 months = 1.029^7.51 = 1.239.
Baseline ADL loss rate p₀ = 0.32 (32%) in the incident-eligible sample (adl_baseline=0).
Δp = 0.32 × 0.239/1.077 ≈ 0.071 (7.1 pp). Reported by `14_pathway_analysis.R`.

**Regime-stratified G<10 (F5):** All four regime-stratified β estimates involve
G = 7–9 clusters. WCB CIs are indicative only; Table 2 rows 13–16 all carry this caveat.

**Germany natural experiment (F2):** PSG II (2017) is a national reform with no
within-Germany control group. Identification uses synthetic control (28-country donor
pool; `09_synthetic_control.R`). Table 1 reflects this design.

**Data access:** Restricted individual-level data (SHARE, ELSA, CHARLS, KLoSA, LASI)
require separate data-use agreement applications. The synthetic demo dataset reproduces
all code paths with simulated data matching the analytic structure.

---

## v14 Changelog (from v11)

| Fix | Description |
|-----|-------------|
| F1 | "14.9%" unit error → "0.149 SD" throughout |
| F2 | Table 1 Germany corrected: SC design, national reform, no within-Germany control |
| F3 | IDS_pure pooled vs within-country SD stated; 7.8-fold basis confirmed symmetric |
| F4 | P1 pathway p₀=0.32 stated; Δp derivation explicit in text and code |
| F5 | G<10 caveat applied consistently to all four regime-stratified β estimates |
| F6 | FEI aggregate-level validity qualification in Results |
| F7 | Exploratory 28% estimate consolidated to one Discussion mention |
| F8 | Delphi 14 countries listed by income group in Supplementary Note 3 |
| F9 | Taiwan high-income status noted in regime classification |
| F10 | SC asymmetry (Germany only) stated as limitation |
| F11 | "illuminates" → "reveals"; removed AI-signature vocabulary |

---

## Citation

[Blinded for peer review]

## Pre-registration

https://osf.io/qn7wz (registered March 2021, prior to data analysis)

## v14 Additional Fixes (from v13)

| Fix | Description |
|-----|-------------|
| F11 | 7.8-fold ratio corrected to approximately 7.6-fold (0.149/0.0196 = 7.61×) |
| F12 | Germany SC permutation p corrected to 0.036 (rank 1 of 28 donor placebos) |
| F13 | Document header corrected to v14 |
| F14 | 'stakeholder' → 'community' in Ethics statement |
| F15 | β_LTC_SD computation (0.019 × 1.03 = 0.0196) made explicit in Results |

## Reference Verification Status (v14)

25 of 26 key references independently verified against journal databases.
One reference (Hu et al. 2023, *Soc. Policy Admin.*) is marked for final DOI
verification before submission; all other details are confirmed.

## v15 Final Fixes (from v14)

| Fix | Description |
|-----|-------------|
| F13 | λ=0.43 "95% CI" mislabelling corrected → "sensitivity range [0.31–0.55]"; pre-specification noted in text and code |
| F14 | β_IV > β_TWFE interpretation added to Results (attenuation bias / LATE) |
| F15 | Ratchet ratio harmonised to 1.68:1 throughout (Introduction was using 1.7:1) |
| F16 | IDS pooled SD stability note: SD stable at 1.00 with/without imputed obs |

## Methodology note: λ calibration

λ = 0.43 in CLS(G) = mean_loss + λ·H_B(G) is **pre-specified** (OSF registration
March 2021). The range [0.31, 0.55] in Supplementary Table S5 is a **post-hoc
sensitivity range** — not a statistical confidence interval. The primary β = 0.149
is stable across this entire range. v14 erroneously labelled this range as
"95% CI"; v15 corrects the description.
