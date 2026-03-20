# Why Care Systems Fail: Institutional Decoupling in Ageing Societies

**Nature Health** | Submitted manuscript | Under review

**Author:** Peng Li, Huachiew Chalermprakiet University  
**Licence:** MIT  
**Analysis documentation:** https://osf.io/5qs9x

---

## Repository Structure

```
ltc-decoupling/
├── R/
│   ├── 00_setup.R                  # Package loading; renv restoration
│   ├── 00_generate_synthetic_demo.R # Synthetic demonstration dataset
│   ├── 01_data_prep.R              # Panel construction (704 obs, 32 countries)
│   ├── 02_ids_construction.R       # IDS_pure, GIS, ESI; λ sensitivity (Table S5)
│   ├── 03_twfe_main.R              # Primary TWFE; entropy decomposition
│   ├── 04_iv_estimation.R          # IV estimation; five-layer exclusion tests
│   ├── 05_mediation.R              # CLS mediation (Imai–Keele–Yamamoto)
│   ├── 06_event_study.R            # Four natural experiments
│   ├── 07_goodman_bacon.R          # Goodman-Bacon decomposition
│   ├── 08_callaway_santanna.R      # Callaway–Sant'Anna GATT
│   ├── 09_synthetic_control.R      # Germany SC; PSM balance check
│   ├── 10_regime_moderation.R      # Regime-stratified regressions
│   ├── 11_falsification.R          # Placebo tests
│   ├── 12_robustness_table.R       # Full Table 2 (12 specifications)
│   ├── 13_figures.R                # Figures 1–4
│   ├── 14_pathway_analysis.R       # P1/P2/P3; baseline ADL rate
│   └── 15_cfa_invariance.R         # Multi-group CFA invariance
├── data/
│   ├── codebook_CCSD.csv           # Comparative Care System Database codebook
│   ├── expert_panel_ratings.csv    # Delphi panel boundary ratings
│   └── synthetic_demo.rds          # Reproducible demo (no restricted data)
├── renv.lock                       # 36 packages, R 4.3.1
└── README.md
```

---

## Reproduction

```r
# 1. Restore package environment
renv::restore()

# 2. Generate synthetic demonstration dataset
source("R/00_generate_synthetic_demo.R")

# 3. Full analysis pipeline
source("R/00_setup.R")
source("R/01_data_prep.R")
source("R/02_ids_construction.R")   # Reports pooled SD=1.00; within-country SD=0.31
source("R/03_twfe_main.R")          # β=0.149; entropy decomposition
source("R/04_iv_estimation.R")      # KP F=34.7; five-layer exclusion triangulation
source("R/05_mediation.R")          # ACME=0.091 SD (62%)
source("R/06_event_study.R")        # Germany, Japan, Korea, China
source("R/14_pathway_analysis.R")   # P1 p₀=0.32 baseline; Δp=7.1pp derivation
```

---

## Key Statistical Notes

**IDS_pure standardisation:** `IDS_pure = z-score(ESI − GIS)` over the pooled
704-observation distribution. Pooled SD = 1.00 by construction; within-country SD
(after country demeaning) = 0.31 SD. Both reported by `02_ids_construction.R`.
The approximately 7.6-fold comparison with LTC expenditure (0.149/0.0196) uses
pooled-SD scales for both variables.

**λ = 0.43 in CLS(G):** Pre-specified value; the range [0.31–0.55] is a
post-hoc sensitivity range tested in `02_ids_construction.R` (Supplementary
Table S5). Primary β is stable across this range.

**P1 pathway absolute risk:** cumOR = 1.029^7.51 = 1.239 (over 7.4 months).
Baseline ADL loss rate p₀ = 0.32 in the incident-eligible sample.
Δp = 0.32 × 0.239/1.077 ≈ 0.071 (7.1 pp). Reported by `14_pathway_analysis.R`.

**Germany natural experiment:** PSG II (2017) is a national reform; identification
uses synthetic control with 28-country donor pool (`09_synthetic_control.R`).
Permutation p = 0.036 (rank 1 of 28 donor placebos).

---

## Data Access

Restricted individual-level data (SHARE, ELSA, CHARLS, KLoSA, LASI) require
separate data-use agreement applications from the originating cohorts. The
synthetic demonstration dataset reproduces all code paths with simulated data.

---

## Citation

Li, P. (2026). Why care systems fail: Institutional decoupling in ageing societies.
*Nature Health* (under review).
