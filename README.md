## Identification of Risk Groups in Pharmacovigilance Using Penalized Regression and Machine Learning (RGP)

`RGP` is an `R` package for analyzing healthcare claims data and 
simulated data using penalized regression and machine learning methods.
This package contains function wrappers to create a simulated cohort,
group predictors based on functional targets (from KEGG and TTD) and
conventional groups (ATC/ICD systems) and analyze the data using various types
of penalized regression (LASSO) and machine learning methods (random forests and
block forests).

### Functionalities

1. cohort simulation (`R/sim_create_cohort.R`)

2. functional target-based grouping - KEGG (`R/ftarget_db_manager.R`)

3. functional target-based grouping - TTD (`R/ftarget_db_manager.R`)

4. penalized regression, group-based analysis and results assessment
(`R/rgp_grpl.R`)

5. classification measures (`R/rgp_classification_measures.R`)

### Structure
* `data`:
  * `ttd_clean_data`: data curated from TTD (drug-target, disease-target.. etc.) as binary matrices
  * `Paper13_events.csv`: file containing ICD list of the ADEs for the project
* `R`:  functions and packages used in the analysis
  * `algorithms-task13.R`: analysis methods wrappers functions
  * `helpers.R`: helper functions
  * `ftsim`: Functional Target Simulation; a full independent package to download, clean and manage TTD and KEGG data. It also use them to create a simple simulated cohort with grouped covariates.
* `README.md`: We are here. We explain the project contents, used packages, and of course the directory structure.
* `analysis`:
  * analysis scripts from 2-8. Run one by one on slave nodes. Steps to obtain healthcare claims data from health insurance databases are omitted for data protection reasons.
  * code to generate descriptive statistics and read performance metrics
  * script to define study variables (inclusion criteria that applies to GePaRD)

See the documentation `?` and `?` for more info.

### Installation

```R
devtools::install_github("bips-hb/rgp")
```

### Usage
```R
```

### Acknowledgements

We gratefully acknowledge the financial support from the innovation fund (“Innovationsfonds”) of the Federal Joint Committee in Germany (grant number: 01VSF16020).

### Contact

Mariam R. Rizkallah\
Leibniz Institute for Prevention Research & Epidemiology - BIPS GmbH
E-mail: rizkallah-issak [at] leibniz-bips [dot] de
