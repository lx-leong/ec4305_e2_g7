# EC4305 Group E2 G7: Medicare at 65

## Research question

Does becoming eligible for Medicare at age 65 increase healthcare use? We measure use with two outcomes:

- **Doctor visit** in the past year (`sawdr`), the main outcome.
- **Overnight hospital stay** in the past year (`inhosp`), the secondary outcome.

We also ask whether the effect is larger for groups that gain more insurance coverage at 65, comparing four race × education groups.

We use a regression discontinuity (RD) design at age 65 with NHIS data for 1992–1996, following Card, Dobkin and Maestas. The **sharp RD** estimates the effect of eligibility. The **fuzzy RD** uses Medicare coverage as the first stage to estimate the effect of having Medicare.

**Main specification:** local linear, ±5-year bandwidth, triangular kernel × survey weights, survey-year fixed effects, HC1 standard errors.

## Folder structure

```
ec4305_e2_g7/
├── data/
│   ├── c9203.sas7bdat                  # Raw replication file, NHIS 1992–2003
│   ├── NHIS_raw.zip                    # Raw NHIS 2000–2003 files (sawdr audit only)
│   ├── nhis_clean_full.csv             # Cleaned 1992–1996 data used by all analyses
│   ├── ec4305_data_cleaning.Rmd
│   ├── ec4305 project EDA.Rmd
│   ├── ec4305 sawdr analysis.Rmd
│   └── output/
│       ├── data_html/                  # HTML reports of the three data Rmds
│       │   ├── ec4305_data_cleaning.html
│       │   ├── ec4305-project-EDA.html
│       │   └── ec4305-sawdr-analysis.html
│       └── sawdr_audit/                # Tables and figure from the sawdr audit
├── src/
│   ├── rd_functions.R                  # Shared RD functions
│   ├── sharp_rd_identification.Rmd
│   ├── sharp_rd_estimation.Rmd
│   └── Fuzzy_RD_1992_1996.Rmd
└── output/
    ├── sharp_rd_identification/        # HTML report + result tables
    ├── sharp_rd_estimation/
    └── fuzzy_rd/
```

## What each file does

### Data (`data/`)

| File | What it does |
|---|---|
| `ec4305_data_cleaning.Rmd` | Reads `c9203.sas7bdat` and keeps 1992–1996. Checks the coding of age, eligibility, weights and binary variables. Adds corrected variables: `poor_clean` (unknown health set to missing) and `lowed_clean` (missing education set to missing, not "below high school"). Keeps rows with missing values, so each model drops only what it needs. Saves `nhis_clean_full.csv`. |
| `ec4305 project EDA.Rmd` | Exploratory analysis of the cleaned data: sample checks, missingness by variable and year, survey-weighted means for ages 60–<65 vs 65–70, age profiles of outcomes and covariates around 65, response rates and sample size by age. |
| `ec4305 sawdr analysis.Rmd` | Explains why the analysis stops at 1996. Doctor-visit answers are 47–59% missing from 1997 on. Using the raw 2000–2003 NHIS files, it shows the supplied coding left out respondents with 10 or more visits, and measures how many answers the correction recovers. Separate from the main analysis. |

### Analysis (`src/`)

| File | What it does |
|---|---|
| `rd_functions.R` | Shared code for the three analysis files: finds and loads `nhis_clean_full.csv`, builds race × education groups, and runs the local linear RD (`fit_rd`) and fuzzy RD (`fit_fuzzy`). It also contains the plotting and table helpers. |
| `sharp_rd_identification.Rmd` | Tests the RD assumptions **without** estimating the outcome effects. It covers eligibility coding, sample flow, whether answers go missing at 65, age heaping, covariate balance, the employment jump at 65, and the insurance gain at 65 by group. |
| `sharp_rd_estimation.Rmd` | Sharp RD estimates of the jump in doctor visits and hospital stays at 65. Robustness checks: bandwidth, functional form, covariates, donut, survey weights and year effects, and placebo cutoffs. It also compares the four subgroups' outcome jumps with their insurance gains. |
| `Fuzzy_RD_1992_1996.Rmd` | Fuzzy RD using Medicare coverage as the first stage. It reports the first stage, reduced form and fuzzy effect, along with inference checks (clustered SEs, rdrobust), covariate balance, bandwidth and specification checks, weight and 1993 sensitivity, placebo cutoffs and subgroup effects. |

## Outputs

Each analysis Rmd saves its HTML report and result tables to its own folder in `output/`. The three data Rmds knit their HTML reports into `data/output/data_html/`.

| Folder | Files |
|---|---|
| `data/output/data_html/` | `ec4305_data_cleaning.html`, `ec4305-project-EDA.html`, `ec4305-sawdr-analysis.html` |
| `output/sharp_rd_identification/` | `sharp_rd_identification.html`, `identification_results.csv` (all identification tests), `sample_flow.csv`, `handoff_subgroups.csv` (insurance gains by group), `model_sample_sizes.csv` |
| `output/sharp_rd_estimation/` | `sharp_rd_estimation.html`, `estimation_results.csv` (main and robustness estimates), `subgroup_results.csv`, `model_sample_sizes.csv` |
| `output/fuzzy_rd/` | `Fuzzy_RD_1992_1996.html`, `fuzzy_results.csv` (all fuzzy estimates), `fuzzy_diagnostics.csv` (balance, placebo and other jump tests), `subgroup_results.csv`, `model_sample_sizes.csv` |
| `data/output/sawdr_audit/` | `original_missingness.csv`, `correction_audit.csv`, `distribution_validation.csv`, `original_sawdr_missingness.png` |

`model_sample_sizes.csv` logs the sample behind every model, for auditing.

## How to run

**R packages:** `haven`, `readr`, `dplyr`, `tidyr`, `ggplot2`, `rdrobust`, `knitr` and `rmarkdown`.

1. Knit `data/ec4305_data_cleaning.Rmd` to rebuild `nhis_clean_full.csv`. This is only needed if the cleaning changes, because the CSV is already in the repo.
2. Knit the analysis files in `src/`. Each one reads `data/nhis_clean_full.csv` directly and doesn't depend on the others, so they can run in any order.
3. The EDA and sawdr audit can be knitted at any time.

The data Rmds use paths relative to `data/`, so knit them from that folder (RStudio's Knit button does this). Their `knit:` header sends the HTML to `data/output/data_html/`; calling `rmarkdown::render()` directly ignores that header and writes the HTML next to the Rmd. The `src/` files find the project root on their own.
