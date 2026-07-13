## ============================================================
## Figure 1: HAL fits a complex function at several sectional
## variation norms (over-regularized / CV-chosen / under-regularized).
##
## Run from this folder:  Rscript fig1_hal.R
## Outputs: fig1_hal.png and fig1_hal.eps (written to this folder).
##
## CREDITS / PRIOR WORK:
##   - The HAL univariate-fit simulation idea follows David Benkeser's
##     HAL simulation code: https://github.com/benkeser/hal (Simulation/).
##   - HAL is provided by the hal9001 package (Coyle, Hejazi, Phillips, et al.);
##     here we use a lightweight glmnet-based implementation of the 0th-order
##     (indicator) HAL basis for a transparent, dependency-light example.
##   - Clean recreation of the figure by Yunwen (Wendy) Ji; restyled here.
## ============================================================

source("hal_helpers.R")   # build_hal_design_matrix, fit_hal_cv, calculate_variation_norms, ...
source("hal_style.R")     # shared style: HAL_COL/LTY/LWD, hal_* helpers

## ---- truth function (smooth, multi-frequency) ----
truth <- function(x) sin(1.2 * x) + 0.7 * sin(2.7 * x)

## The 0th-order HAL basis is I(x >= knot), with a knot at every observed x.
## The outermost knots are supported by only a handful of points, so the fitted
## step function is wildly unstable in the last percentile at each end. That is a
## boundary artifact of the basis, not a property of HAL worth illustrating here,
## so the FITTED CURVES are drawn only over the interior TRIM quantile range.
## The data cloud is still plotted over its full range.
TRIM <- 0.005   # middle 99% of x (0.5% trimmed from each end)

compute_fig1 <- function(truth, n = 500, xmin = -4, xmax = 4, noise = 0.8,
                         low_mult = 0.15, high_mult = 5, seed = 123) {
  set.seed(seed)
  x  <- runif(n, xmin, xmax)
  y  <- truth(x) + rnorm(n, 0, noise)

  ## fitted curves live on the interior grid only (see TRIM note above)
  lim <- unname(quantile(x, c(TRIM, 1 - TRIM)))
  xg  <- seq(lim[1], lim[2], length.out = 400)

  cv     <- fit_hal_cv(build_hal_design_matrix(x), y)
  norms  <- calculate_variation_norms(cv)
  cvnorm <- norms[which(cv$lambda == cv$lambda.min)]
  tgt    <- find_lambdas_by_norm(cv, c(low_mult, high_mult) * cvnorm)

  dn   <- build_hal_design_matrix(x, xg)
  lams <- c(low = tgt[[1]]$lambda, cv = cv$lambda.min, high = tgt[[2]]$lambda)
  P    <- predict_hal_multiple(cv$glmnet.fit, dn, lams)

  list(x = x, y = y, xg = xg, P = P,
       achieved = c(tgt[[1]]$achieved_norm, cvnorm, tgt[[2]]$achieved_norm))
}

draw_fig1 <- function(d, pts_col = HAL_PTS) {
  hal_par()
  plot(d$x, d$y, pch = 16, col = pts_col, cex = 0.6,
       xlab = "x", ylab = "", bty = "l", xaxt = "n")
  hal_axes(ylab = "Y")
  hal_line(d$xg, truth(d$xg),   "truth")
  hal_line(d$xg, d$P[, "low"],  "lo")
  hal_line(d$xg, d$P[, "high"], "hi")
  hal_line(d$xg, d$P[, "cv"],   "cv")        # CV-chosen drawn last (on top)
  hal_legend("topright", d$achieved)
}

d <- compute_fig1(truth)

hal_png("fig1_hal.png")
draw_fig1(d); dev.off()

## EPS: cairo supports alpha, classic postscript does not, so the point cloud
## falls back to a solid grey when cairo is unavailable.
has_alpha <- hal_eps("fig1_hal.eps")
draw_fig1(d, pts_col = if (has_alpha) HAL_PTS else HAL_PTS_EPS)
dev.off()

cat(sprintf("Figure 1 done. Norms: low=%.1f  CV=%.1f  high=%.1f\n",
            d$achieved[1], d$achieved[2], d$achieved[3]))
