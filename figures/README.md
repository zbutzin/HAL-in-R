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
| `hal_style.R` | Shared plotting style (palette, line types, axes, figure size/font, EPS device) so all figures match. |

## Reading the figures

**Figure 1.** The 0th-order HAL basis places a knot `I(x >= x_i)` at every observed
point, so the outermost knots are supported by only a handful of observations and
the fit is unstable in the last percentile at each end. That is an artifact of the
basis, not a property of HAL worth illustrating, so the *fitted curves* are drawn
only over the interior middle 99% of `x` (`TRIM` in `fig1_hal.R`); the data cloud
is still shown over its full range.

**Figure 2.** Three regions, left to right:

- the **plateau** (small lambda) — the width is flat, i.e. insensitive to lambda;
  this is the region the bootstrap lambda is selected from, and it is shaded;
- an **overshoot** just above the plateau, where the width bumps up before falling.
  This is a real, smooth feature of the influence-curve width, not noise (see below);
- a **degenerate shelf** (large lambda) where the over-regularized HAL density fit
  collapses and the Wald interval has *exactly* zero width. These points carry no
  information — they are not unusually precise intervals — so they are drawn as open
  grey symbols to set them apart from the informative (filled) points.

The selected lambda is marked with a dashed line. It comes from the package's own
`select_lambda_pleateau_wald()`, applied to the informative (nonzero-width) points
only — the curve has two flat regions, and the plateau finder will latch onto the
degenerate zero shelf if it is handed the full grid.

Note that the Wald width here is influence-curve based (`TMLEbootstrap` computes it
as `diff(pointTMLE$CI)`), so it is **deterministic given lambda** — there is no Monte
Carlo noise in this curve, and `n_bootstrap` does not enter it. The lambda grid
therefore controls how well the curve is *resolved*, not how noisy it is: on a coarse
grid the overshoot looks like a single stray point, and on the denser grid used here
it resolves into the smooth bump it actually is.

## How to run

From **this** folder:

```sh
Rscript fig1_hal.R          # writes fig1_hal.png / .eps
Rscript fig2_density_width.R # writes fig2_density_width.png / .eps / .rds
```

`fig1_hal.R` needs only `glmnet` (and base R). `fig2_density_width.R` needs the
`TMLEbootstrap` package (and `ggplot2`). See the dependency note below.

Both scripts write EPS via cairo when it is available and fall back to the classic
`postscript` device otherwise. The fallback is probed by actually opening a device,
because `capabilities("cairo")` can report `TRUE` on a machine where `cairo_ps()`
still fails to load (e.g. macOS without XQuartz) — in which case the earlier scripts
silently produced no `.eps` at all.

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
