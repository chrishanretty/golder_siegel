# Replication of Golder, Golder & Siegel (2012) in R

This is an R replication of the zero-intelligence government formation model in

> Golder, Golder & Siegel (2012), "Modeling the Institutional Foundation of
> Parliamentary Government Formation", *JOP* 74(2).

`ABM_summary.md` summarises the model and how ambiguities were resolved.

`EVALUATION.md` sets out how well the replication matches the published statistics and diagnoses the residual gaps.

## Files

| Path | Purpose |
|---|---|
| `R/model.cpp` | The government-formation loop: Rcpp, using R's RNG. |
| `R/ggs_model.R` | Election stage, diversity and ENPP measures, a pure-R reference version of the loop, and `ggs_one_run()` for a single run. |
| `R/harness.R` | Parallel batch runner. It uses L'Ecuyer RNG streams per chunk, checkpointed chunks, `mclapply` on Unix and PSOCK elsewhere. |
| `R/analysis.R` | ENPP weights and Tables 1–4. Run on the authors' deposited data, it reproduces the published numbers. |
| `R/compare.R` | Compares a batch with the GGS simulated data, overall and conditional on *m* and β. |
| `R/report.R` | Writes `outputs/results.md` and `outputs/fig1-3.png`. |
| `R/experiments.R` | Quick experiments on implementation variants, used to diagnose the residual gaps. |
| `tests/test_model.R` | Checks that the C++ core and R reference are identical, plus classification and invariant tests. |
| `data/original/` | The authors' Dataverse deposit (doi:10.7910/DVN/KIFRVV): simulated output, real-world data and figure script. |

## Reproduce

Run these from the project root:

```sh
Rscript tests/test_model.R
Rscript R/harness.R 50000 outputs/sim_main 8     # ~23 min on 8 cores
Rscript R/report.R outputs/sim_main.rds outputs/
Rscript R/harness.R 20000 outputs/sim_linf 8 Linf   # sensitivity: Chebyshev voting
```

To run a single system interactively:

```r
source("R/ggs_model.R"); load_ggs_cpp()
set.seed(1); res <- ggs_one_run()
```
