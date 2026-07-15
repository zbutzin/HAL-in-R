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
##   This script requires the TMLEbootstrap package. Install the fork:
##     remotes::install_github("zbutzin/TMLEbootstrap",
##                             ref = "fix-hal9001-weights-and-cleanup")
##   The fork is needed because upstream errors against hal9001 >= 0.4.0, but it
##   does NOT fit the density correctly on its own: it passes the observation
##   weights inside fit_control, where hal9001::fit_hal() overwrites them with
##   NULL. So we also source tmleboot_patch.R below, which puts the weights back
##   where fit_hal expects them.
##
##   DO NOT REMOVE THAT source() LINE. Without it the fitted density comes out
##   nearly flat, and because the Wald interval is influence-curve based the CI
##   widths collapse to ~1/40th of their true size (plateau 3.6e-4 instead of
##   1.5e-2). The y axis below is auto-scaled, so the figure will still render
##   perfectly and give you no visual warning. The Psi check further down is what
##   catches it. See tmleboot_patch.R and figures/README.md.
##
## NOTE ON THE LAMBDA GRID:
##   The Wald width here is influence-curve based (TMLEbootstrap computes it as
##   diff(pointTMLE$CI)); it is deterministic given lambda, with no Monte Carlo
##   noise. So the grid controls how well the curve is RESOLVED, not how noisy it
##   is. The earlier 20-point grid made the overshoot just left of the plateau
##   look like a single stray point; on the denser grid below it resolves into
##   the smooth bump it actually is.
## ============================================================

## ---- load TMLEbootstrap (see DEPENDENCY NOTE above to install) ----
if (!requireNamespace("TMLEbootstrap", quietly = TRUE)) {
  stop("TMLEbootstrap not installed. Install the fixed fork:\n",
       "  remotes::install_github(\"zbutzin/TMLEbootstrap\", ref = \"fix-hal9001-weights-and-cleanup\")")
}
suppressMessages(library(TMLEbootstrap))
suppressMessages(library(ggplot2))
source("tmleboot_patch.R")   # restores the observation weights in densityHAL$fit
source("hal_style.R")        # shared palette (HAL_COL, HAL_ACCENT) + geometry + theme

## ---- settings ----
N_SIM <- 2000; N_MODE <- 3                      # 3-modal truth at n = 2000
N_LAMBDA <- 50; LAMBDA_MAX <- -1; LAMBDA_MIN <- -8   # denser + wider than before
BIN_WIDTH <- 0.5; EPS_STEP <- 0.01; SEED <- 123

set.seed(SEED)
data_out <- simulate_density_data(n_sim = N_SIM, n_mode = N_MODE)
tuner <- avgDensityTuneHyperparam$new(
  data = data_out, bin_width = BIN_WIDTH, epsilon_step = EPS_STEP
)
lambda_grid <- 10^seq(LAMBDA_MAX, LAMBDA_MIN, length.out = N_LAMBDA)
invisible(capture.output(
  tuner$add_lambda(lambda_grid = lambda_grid, to_parallel = FALSE)
))

df_w <- tuner$get_lambda_df()   # columns: lambda, width, kindCI

## ---- degenerate region -----------------------------------------------------
## At large lambda the HAL density fit collapses and the Wald interval degenerates
## to *exactly* zero width. Those points are not "very precise" estimates -- they
## carry no information -- so they are kept in the data but drawn as open symbols,
## with the connecting line broken across them.
df_w$degenerate <- df_w$width <= 0
df_ok  <- df_w[!df_w$degenerate, ]
df_bad <- df_w[ df_w$degenerate, ]

## ---- plateau + selected lambda ---------------------------------------------
## Select on the informative points only. The curve has TWO flat regions -- the
## real plateau at small lambda and the degenerate zero shelf at large lambda --
## and the plateau finder would happily latch onto the zero shelf if it were fed
## the full grid.
lambda_sel <- tuner$select_lambda_pleateau_wald(df_ok)
plateau_lo <- min(df_ok$lambda)

cat(sprintf("Plateau starts at lambda = %.3g (width = %.3g)\n",
            lambda_sel, df_ok$width[which.min(abs(df_ok$lambda - lambda_sel))]))
cat(sprintf("Degenerate (zero-width) for lambda >= %.3g : %d of %d grid points\n",
            min(df_bad$lambda), nrow(df_bad), nrow(df_w)))

## ---- CV-selected lambda ----------------------------------------------------
## hal9001's own cross-validation (cv_select = TRUE), fit on the same frequency-
## compressed bin-rows the density HAL uses. On this simulation it lands at the
## MAXIMUM of glmnet's lambda path -- a maximally-regularized, essentially null
## density -- so it sits to the RIGHT of the plateau, in the zero-width region.
## This is the "cross-validation choice" of Demo.Rmd: undersmoothing means moving
## lambda LEFT from here until the Wald width stabilizes on the plateau.
ld_cv  <- get("longiData", envir = asNamespace("TMLEbootstrap"))$new(
  x = data_out$x, bin_width = BIN_WIDTH)
dfc_cv <- ld_cv$generate_df_compress(x = data_out$x)
set.seed(SEED)
hal_cv <- hal9001::fit_hal(
  X = dfc_cv[, "box"], Y = dfc_cv$Y, family = "binomial",
  weights = dfc_cv$Freq, fit_control = list(cv_select = TRUE),
  return_lasso = TRUE, return_x_basis = FALSE, yolo = FALSE)
lambda_cv <- hal_cv$lasso_fit$lambda.min
cat(sprintf("CV-selected lambda = %.3g\n", lambda_cv))

## ---- sanity check against the closed-form truth ----------------------------
## The y axis is auto-scaled, so a badly wrong CI scale still renders as a clean
## figure. This is the check that catches it. simulate_density_data() draws from
## an equal-weight Gaussian mixture with modes evenly spaced on [-4, 4] and
## sd = 10 / n_mode / 6, so Psi = Integral p(x)^2 dx has a closed form: with
## p = (1/K) sum_k N(m_k, s), Integral p^2 = (1/K^2) sum_{j,k} N(m_j - m_k; 0, s*sqrt(2)).
## If the observation weights are ever dropped again, Psi halves and this shouts.
modes    <- seq(-4, 4, length.out = N_MODE)
sigma    <- 10 / N_MODE / 6
psi_true <- mean(outer(modes, modes, function(a, b) dnorm(a - b, 0, sigma * sqrt(2))))
psi_hat  <- tuner$dict_boot[[which.min(df_w$lambda)]]$bootOut$pointTMLE$Psi
cat(sprintf("Psi: fitted = %.4f vs truth = %.4f (ratio %.2f)\n",
            psi_hat, psi_true, psi_hat / psi_true))
if (abs(psi_hat / psi_true - 1) > 0.15) {
  warning("Fitted Psi is far from the closed-form truth, so the CI widths are not ",
          "trustworthy. Is source(\"tmleboot_patch.R\") still at the top of this ",
          "script? Without it the observation weights are dropped. See tmleboot_patch.R.")
}

## ---- plot ------------------------------------------------------------------
## U is the y-axis tick spacing, and doubles as the unit the text annotations are
## positioned in. It is set from the data (rather than a fixed 1e-4) so the figure
## renders correctly regardless of the CI-width scale the local package versions
## produce -- some environments give density widths ~1e-2, others ~1e-4. The
## plateau band, selected-lambda marker, and labels all stay proportional. The
## y axis keeps R's default scientific labels (e.g. 0e+00, 3e-03, ...).
U <- max(df_w$width) / 4.45

p <- ggplot(mapping = aes(x = lambda, y = width)) +
  ## the plateau: the whole point of the figure
  annotate("rect", xmin = plateau_lo, xmax = lambda_sel,
           ymin = -Inf, ymax = Inf, fill = HAL_ACCENT[["band"]]) +
  ## two vertical lambda markers: the plateau lambda (used for the bootstrap) and
  ## the CV-selected lambda (to its right). Distinct macaron colors so they read
  ## as two different selectors, not one.
  annotate("segment", x = lambda_sel, xend = lambda_sel, y = -Inf, yend = Inf,
           color = HAL_ACCENT[["lam_plateau"]], linetype = 2, linewidth = 0.7) +
  annotate("segment", x = lambda_cv, xend = lambda_cv, y = -Inf, yend = Inf,
           color = HAL_ACCENT[["lam_cv"]], linetype = 2, linewidth = 0.7) +
  ## the curve runs across the whole grid, including down to the zeros
  geom_line(data = df_w, color = HAL_COL[["cv"]], linewidth = 0.9) +
  ## informative points: filled
  geom_point(data = df_ok, color = HAL_COL[["cv"]], size = 1.9) +
  ## degenerate (zero-width) points: open circles, so they read as "no information
  ## here" rather than as an extremely precise interval
  geom_point(data = df_bad, color = HAL_ACCENT[["dead"]], size = 1.9,
             shape = 1, stroke = 0.7) +
  ## labels (annotation heights are given in units of U, so they stay put)
  annotate("text", x = sqrt(plateau_lo * lambda_sel), y = 0.55 * U,
           label = "plateau", color = HAL_ACCENT[["band_txt"]],
           size = 3.4, fontface = 2) +
  ## marker labels sit ABOVE the curve (two lines each, centred on their line) so
  ## neither crosses the width curve or the other. Values shown to 2 sig figs.
  annotate("text", x = lambda_sel, y = 5.3 * U, hjust = -0.1, vjust = 1,
           label = sprintf("atop(plateau~lambda, %s)",
                           format(lambda_sel, digits = 2, scientific = FALSE)),
           parse = TRUE, color = HAL_ACCENT[["lam_plateau"]], size = 3.0) +
  annotate("text", x = lambda_cv, y = 5.3 * U, hjust = -0.1, vjust = 1,
           label = sprintf("atop(CV*'-'*selected~lambda, %s)",
                           format(lambda_cv, digits = 2, scientific = FALSE)),
           parse = TRUE, color = HAL_ACCENT[["lam_cv"]], size = 3.0) +
  ## degenerate-shelf note, centred over the open circles and right of the CV line
  annotate("text", x = sqrt(min(df_bad$lambda) * max(df_bad$lambda)), y = 0.35 * U,
           hjust = 0.2, label = "CI collapses\nto zero width",
           color = HAL_ACCENT[["dead"]], size = 3.0, lineheight = 0.95) +
  scale_x_log10(
    breaks = 10^seq(LAMBDA_MIN, LAMBDA_MAX),
    labels = parse(text = paste0("10^", seq(LAMBDA_MIN, LAMBDA_MAX)))
  ) +
  ## floor sits just below 0 so the zero-width circles are not bisected by the axis
  scale_y_continuous(limits = c(-0.18 * U, 5.4 * U), labels = scales::label_scientific(),
                     expand = c(0, 0)) +
  labs(x = expression(lambda ~ "(log scale; larger" ~ lambda ~ "= more regularization)"),
       y = "Wald CI width") +
  hal_ggtheme()

ggsave("fig2_density_width.png", p, width = FIG_W, height = FIG_H, dpi = FIG_DPI)
ggsave("fig2_density_width.eps", p, width = FIG_W, height = FIG_H,
       device = hal_eps_device())
saveRDS(df_w, "fig2_density_width.rds")
cat("Figure 2 saved.\n")
