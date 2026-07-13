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
##   This script requires the TMLEbootstrap package. With current hal9001
##   (>= 0.4.0), upstream TMLEbootstrap needs a small fix in densityHAL.R
##   (observation weights passed inside fit_control). Install the fixed fork:
##     remotes::install_github("zbutzin/TMLEbootstrap",
##                             ref = "fix-hal9001-weights-and-cleanup")
##   See figures/README.md for details.
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
source("hal_style.R")   # shared palette (HAL_COL, HAL_ACCENT) + geometry + theme

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
  annotate("segment", x = lambda_sel, xend = lambda_sel, y = -Inf, yend = Inf,
           color = HAL_ACCENT[["pick"]], linetype = 2, linewidth = 0.5) +
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
  annotate("text", x = lambda_sel, y = 4.45 * U, hjust = -0.12,
           label = sprintf("selected~lambda == %.1e", lambda_sel),
           parse = TRUE, color = HAL_ACCENT[["pick"]], size = 3.2) +
  annotate("text", x = max(df_bad$lambda), y = 0.05 * U, hjust = 2,
           label = "CI collapses\nto zero width",
           color = HAL_ACCENT[["dead"]], size = 3.0, lineheight = 0.95) +
  scale_x_log10(
    breaks = 10^seq(LAMBDA_MIN, LAMBDA_MAX),
    labels = parse(text = paste0("10^", seq(LAMBDA_MIN, LAMBDA_MAX)))
  ) +
  ## floor sits just below 0 so the zero-width circles are not bisected by the axis
  scale_y_continuous(limits = c(-0.18 * U, 5 * U), labels = scales::label_scientific(),
                     expand = c(0, 0)) +
  labs(x = expression(lambda ~ "(log scale; larger" ~ lambda ~ "= more regularization)"),
       y = "Wald CI width") +
  hal_ggtheme()

ggsave("fig2_density_width.png", p, width = FIG_W, height = FIG_H, dpi = FIG_DPI)
ggsave("fig2_density_width.eps", p, width = FIG_W, height = FIG_H,
       device = hal_eps_device())
saveRDS(df_w, "fig2_density_width.rds")
cat("Figure 2 saved.\n")

