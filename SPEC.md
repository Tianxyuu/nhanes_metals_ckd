# NHANES metals & CKD (2017-2018)

Goal: Cross-sectional analysis of blood lead/cadmium and kidney function.

- Cycle: 2015-2016 (NHANES I) + 2017-2018 (NHANES J)
- Exposure: blood lead (LBXBPB), blood cadmium (LBXBCD), log-transformed
- Outcome: eGFR from serum creatinine (LBXSCR), age (RIDAGEYR), sex (RIAGENDR)
- Optional: CKD indicator (eGFR < 60)
- Covariates (minimal): age, sex
- Optional covariates: race/ethnicity, education, PIR, BMI, smoking
- Models: linear (eGFR), logistic (CKD)
- Survey weights: SDMVPSU, SDMVSTRA, WTMEC2YR

Artifacts:
- R scripts under `R/`
- R Markdown report under `report/`
- Cached data under `data/` (optional)
