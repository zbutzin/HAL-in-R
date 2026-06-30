## ============================================================
## Locked plotting style for the HAL tutorial figures (Figs 1-3)
## Source this in every figure script for a consistent look.
##
##   keys:  truth | lo (over-regularized) | cv (CV-chosen) | hi (under-regularized)
##   design: grayscale-safe, so each line differs in BOTH color and line type.
## ============================================================

HAL_COL <- c(truth = "grey35", lo = "#d95f02", cv = "#2c7fb8", hi = "#1b9e77")
HAL_LTY <- c(truth = 2,        lo = 3,         cv = 1,         hi = 4)
HAL_LWD <- c(truth = 2,        lo = 2,         cv = 3.5,       hi = 1.8)

## data-point color. PNG uses semi-transparency; EPS (no alpha) passes a solid grey.
HAL_PTS  <- adjustcolor("grey55", 0.45)
HAL_PTS_EPS <- "grey70"

## common margins + axis-title positioning
hal_par <- function() par(mar = c(4.2, 5.2, 1.2, 1), mgp = c(3, 0.8, 0))

## L-shaped axes: x labels with NO tick marks; y title kept clear of the edge
hal_axes <- function(ylab = "Y", yline = 2.6) {
  axis(1, tick = FALSE)
  title(ylab = ylab, line = yline)
}

## draw one named fit ("truth"/"lo"/"cv"/"hi") over grid xg
hal_line <- function(xg, yv, key) {
  lines(xg, yv, col = HAL_COL[[key]], lty = HAL_LTY[[key]], lwd = HAL_LWD[[key]])
}

## standard legend for the 4 elements (pass the achieved norms for lo/cv/hi)
hal_legend <- function(pos = "topright", norms) {
  legend(pos, bty = "n", cex = 0.85, seg.len = 2.6,
         legend = c("Truth",
                    sprintf("HAL, low norm  (%.1f)",  norms[1]),
                    sprintf("HAL, CV-chosen (%.1f)",  norms[2]),
                    sprintf("HAL, high norm (%.1f)",  norms[3])),
         col = HAL_COL[c("truth","lo","cv","hi")],
         lty = HAL_LTY[c("truth","lo","cv","hi")],
         lwd = HAL_LWD[c("truth","lo","cv","hi")])
}
