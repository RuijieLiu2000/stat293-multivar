# Multi-subject VAR simulations

This project compares multi-VAR, adaptive multi-VAR, individual LASSO, and pooled LASSO across **16 conditions, with five replications per condition**. Start with [FinalReport.pdf](FinalReport.pdf) for the simulation report.

## File guide

### Setup, analysis, and reports

| File | Purpose |
|---|---|
| [001_requirements.R](001_requirements.R) | Install missing packages and load them. |
| [002_requirements_lite.R](002_requirements_lite.R) | Load already installed packages. |
| [00_load_functions.R](00_load_functions.R) | Load packages and all function files. |
| [10_run_simulations.R](10_run_simulations.R) | Run the full study and save figures and tables. |
| [11_RMD_PreliminarySims.Rmd](11_RMD_PreliminarySims.Rmd) | Fit a small multi-VAR example, then use its estimates to fit adaptive multi-VAR. |
| [21_RMD_16ConditionsComparisons.Rmd](21_RMD_16ConditionsComparisons.Rmd) | Compare all 16 conditions using the saved results. |
| [31_RMD_FLIPDataInventory.Rmd](31_RMD_FLIPDataInventory.Rmd) | Scan the FLIP real data. |
| [flip/](flip/) | Real data from the FLIP daily diary study. |
| [FinalReport.Rmd](FinalReport.Rmd) | Source for the simulation report. |
| [FinalReport.pdf](FinalReport.pdf) | Rendered simulation report. |
| [90_render_report.R](90_render_report.R) | Rebuild the report PDF from saved results. |

### Function files

| File | Purpose |
|---|---|
| [01_funct_GeneratingSimData.R](R/01_funct_GeneratingSimData.R) | Generate coefficients and simulated series. |
| [02_funct_RegressionDesign.R](R/02_funct_RegressionDesign.R) | Construct regression data and penalty weights. |
| [03_funct_MethodImplementation.R](R/03_funct_MethodImplementation.R) | Fit models and construct adaptive weights. |
| [04_funct_RollingValidation.R](R/04_funct_RollingValidation.R) | Select penalties and refit models. |
| [05_funct_ModelEvaluation.R](R/05_funct_ModelEvaluation.R) | Evaluate recovery and forecasts; summarize replications. |
| [06_funct_SimulationStudy.R](R/06_funct_SimulationStudy.R) | Define the 16 conditions and run simulations. |
| [07_funct_ResultsPresentation.R](R/07_funct_ResultsPresentation.R) | Export figures and formatted tables. |

### Results

- [results/figures/](results/figures/): four figures covering relative Frobenius error, one-step RMSFE, sensitivity, and specificity.
- [results/tables/](results/tables/): seven compact tables covering the design, four performance measures, and two paired comparisons against individual LASSO.

## Run the analysis

Set the working directory to `simulations_clean` in R or RStudio.

**Install packages once:**

```r
source("001_requirements.R")
```


**Run all 16 conditions and rebuild the report:**

```r
source("00_load_functions.R")
study <- run_study(profile = "preliminary", reps = 5)
source("90_render_report.R")
```

This replaces the figures and tables in `results/` and refreshes `FinalReport.pdf`. To rebuild only the report from existing results, run `source("90_render_report.R")` alone.

**Knit the analysis Rmd files:**

Rmd 11 generates and fits its small example; Rmd 21 displays saved results without rerunning the full study; Rmd 31 is scan of the FLIP real data.

## Read the results

- **Relative Frobenius error:** coefficient-recovery error relative to the true matrix; lower is better.
- **RMSFE:** forecast error; lower is better.
- **Sensitivity:** proportion of true nonzero coefficients recovered; higher is better.
- **Specificity:** proportion of true zeros correctly excluded; higher is better.

Tables report means with Monte Carlo standard errors. Figure bars show 95% Monte Carlo intervals. In paired tables, negative differences favor the method over individual LASSO.

## Reference

Fisher et al. (2022). *Penalized Estimation and Forecasting of Multiple Subject Intensive Longitudinal Data*. Psychometrika, 87, 403–431. (https://doi.org/10.1007/s11336-021-09825-7).
