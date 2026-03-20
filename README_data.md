# Data Sources and Access Instructions

## Primary Data Sources

### Individual-level cohort data (obtain directly from providers)

| Study | Coverage | Key variables used | Access URL |
|-------|----------|--------------------|------------|
| SHARE (waves 1–9) | 28 European countries + Israel | ADL/IADL, hospitalisation, care need level | https://share-project.org/data-access |
| ELSA (waves 1–10) | England | ADL/IADL, hospitalisation | https://ukdataservice.ac.uk |
| CHARLS (waves 1–5) | China | ADL/IADL, care utilisation | https://charls.pku.edu.cn |
| KLoSA (waves 1–9) | South Korea | ADL/IADL, care need classification | https://keis.or.kr |
| LASI (wave 1) | India | ADL/IADL, care utilisation | https://iipsindia.ac.in |

### Country-level governance data (publicly available)

| Source | Variables | DOI/URL |
|--------|-----------|---------|
| OECD LTC Statistics 2023 | ltc_exp_gdp_pct, ltc_recipients_rate | https://doi.org/10.1787/7a7afb35-en |
| World Bank WGI | govt_effectiveness | https://info.worldbank.org/governance/wgi/ |
| ICRG | icrg_quality | https://www.prsgroup.com/icrg |
| EUROBAROMETER | Care attitudes (Rounds 71,76,81,87) | https://europa.eu/eurobarometer |

### CCSD (Comparative Care System Database)

The CCSD is an original harmonised governance panel constructed for this study.
Variable-level documentation is in `codebook_CCSD.csv`.
The full CCSD dataset is available from the corresponding author upon reasonable request.

## Synthetic Demonstration Dataset

`synthetic_demo.rds` is a synthetic dataset generated from multivariate normal
draws calibrated to the descriptive statistics in Table S1. It has the same
column structure as the analytic panel and allows all code to run without
access to the primary data. **It contains no real individual-level data and
results will NOT match the paper.**

To generate a new synthetic demo dataset:
```r
source("R/00_setup.R")
set.seed(42)
# See generate_synthetic_demo.R for full generation code
```
