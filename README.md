# Why Care Systems Fail: Institutional Decoupling in Ageing Societies

**Nature Health** | Under blind review

**OSF pre-registration:** [https://osf.io/5qs9x]  
**Licence:** CC-By Attribution 4.0 International

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

## Key Technical Notes 

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

## Citation

[Blinded for peer review]

## Pre-registration

https://osf.io/5qs9x





