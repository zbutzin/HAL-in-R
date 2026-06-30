## ============================================================
## Figure 2: Wald confidence-interval WIDTH as a function of the HAL L1
## penalty (lambda), for the average-squared-density parameter.
## Shows the width "plateau" used to select lambda for the nonparametric
## bootstrap. (The derivative panel of the original two-panel figure is
## intentionally omitted; we keep only the width curve.)
##
## Run from this folder:  Rscript fig2_density_width.R
## Outputs: fig2_density_width.png/.eps and fig2_density_width.rds
##
## CREDITS / PRIOR WORK:
##   - Method and simulation machinery from the TMLEbootstrap package by
##     Weixin Cai (https://github.com/wilsoncai1992/TMLEbootstrap);
##     Cai & van der Laan (2020), "Nonparametric bootstrap inference for
##     the targeted HAL-LASSO."
##   - Package edits / fork by Yunwen (Wendy) Ji. Restyled here.
##
## DEPENDENCY NOTE (important):
##   This script requires the TMLEbootstrap package. As written for current
##   hal9001 (>= 0.4.0), TMLEbootstrap needs a small fix in densityHAL.R:
##   observation `weights` must be passed inside `fit_control`, not as a
##   top-level argument to hal9001::fit_hal(). Use a fork that includes this
##   fix. See figures/README.md for details.
## ============================================================

## ---- locate / load TMLEbootstrap ----
## Option A: a locally checked-out (fixed) copy via devtools:
##   PKGDIR <- "/path/to/TMLEbootstrap"; devtools::load_all(PKGDIR)
## Option B: an installed (fixed) copy:
##   library(TMLEbootstrap)
if (!requireNamespace("TMLEbootstrap", quietly = TRUE)) {
  stop("TMLEbootstrap not installed. See the DEPENDENCY NOTE / figures/README.md.")
}
suppressMessages(library(TMLEbootstrap))
suppressMessages(library(ggplot2))
source("hal_style.R")   # shared palette (HAL_COL)

## ---- settings ----
N_SIM <- 2000; N_MODE <- 3; N_BOOTSTRAP <- 200   # 3-modal truth at n=2000 (distinct plateau)
N_LAMBDA <- 20; LAMBDA_MAX <- -1; LAMBDA_MIN <- -7   # wide enough to show the plateau
BIN_WIDTH <- 0.5; EPS_STEP <- 0.01; SEED <- 123

set.seed(SEED)
data_out <- simulate_density_data(n_sim = N_SIM, n_mode = N_MODE)
tuner <- avgDensityTuneHyperparam$new(
  data = data_out, bin_width = BIN_WIDTH,
  epsilon_step = EPS_STEP, n_bootstrap = N_BOOTSTRAP
)
lambda_grid <- 10^seq(LAMBDA_MAX, LAMBDA_MIN, length.out = N_LAMBDA)
tuner$add_lambda(lambda_grid = lambda_grid, to_parallel = FALSE)

df_w <- tuner$get_lambda_df()   # columns: lambda, width, kindCI
print(df_w)

p <- ggplot(df_w, aes(x = lambda, y = width)) +
  geom_line(color = HAL_COL[["cv"]], linewidth = 1) +
  geom_point(color = HAL_COL[["cv"]], size = 2) +
  scale_x_log10() +
  labs(x = expression(lambda~"(log scale)"), y = "Wald CI width") +
  theme_classic(base_size = 13)

ggsave("fig2_density_width.png", p, width = 7, height = 5, dpi = 300)
ggsave("fig2_density_width.eps", p, width = 7, height = 5, device = cairo_ps)
saveRDS(df_w, "fig2_density_width.rds")
cat("Figure 2 saved.\n")
