# Figure-generation code

This folder contains the R scripts that generate the figures in the
accompanying HAL tutorial manuscript, provided here for transparency and
reproducibility. The code is original to this project but builds directly on
prior work, credited below.

## Files

| File | Purpose |
|------|---------|
| `fig1_hal.R` | **Figure 1**: HAL fitting a complex function at several sectional variation norms (over-regularized / CV-chosen / under-regularized). |
| `fig2_density_width.R` | **Figure 2**: Wald CI width vs. the HAL L1 penalty (lambda) for the average-squared-density parameter; shows the width "plateau". |
| `hal_helpers.R` | Lightweight glmnet-based HAL helpers used by `fig1_hal.R`. |
| `hal_style.R` | Shared plotting style (palette, line types, axes) so all figures match. |

## How to run

From **this** folder:

```sh
Rscript fig1_hal.R          # writes fig1_hal.png / .eps
Rscript fig2_density_width.R # writes fig2_density_width.png / .eps / .rds
```

`fig1_hal.R` needs only `glmnet` (and base R). `fig2_density_width.R` needs the
`TMLEbootstrap` package. See the dependency note below.

## Dependency note for Figure 2

`fig2_density_width.R` uses the `TMLEbootstrap` package. With current
`hal9001` (>= 0.4.0), the upstream `TMLEbootstrap` needs a small fix in
`R/densityHAL.R` (observation `weights` must be passed **inside** `fit_control`,
not as a top-level argument to `hal9001::fit_hal()`). A fork with this fix is
available, so you can install it directly:

```r
remotes::install_github("zbutzin/TMLEbootstrap", ref = "fix-hal9001-weights-and-cleanup")
```

Then run `fig2_density_width.R` as-is (it calls `library(TMLEbootstrap)`).

## Credits / prior work

- **Figure 1** follows David Benkeser's HAL simulation code
  (<https://github.com/benkeser/hal>, `Simulation/`) and the `hal9001`
  package (Coyle, Hejazi, Phillips, et al.,
  <https://github.com/tlverse/hal9001>). The clean recreation used here was
  developed by Yunwen (Wendy) Ji.
- **Figure 2** is based on the `TMLEbootstrap` package by Weixin Cai
  (<https://github.com/wilsoncai1992/TMLEbootstrap>) and
  Cai & van der Laan (2020), *Nonparametric bootstrap inference for the
  targeted highly adaptive LASSO estimator*. Package edits by Yunwen (Wendy) Ji.

The figures were recreated and restyled for this tutorial manuscript.
