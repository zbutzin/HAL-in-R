## ============================================================
## Locked plotting style for the HAL tutorial figures (Figs 1-2)
## Source this in every figure script for a consistent look.
##
##   keys:  truth | lo (over-regularized) | cv (CV-chosen) | hi (under-regularized)
##   design: grayscale-safe, so each line differs in BOTH color and line type.
## ============================================================

HAL_COL <- c(truth = "grey35", lo = "#d95f02", cv = "#2c7fb8", hi = "#1b9e77")
HAL_LTY <- c(truth = 2,        lo = 3,         cv = 1,         hi = 4)
HAL_LWD <- c(truth = 2,        lo = 2.6,       cv = 3.5,       hi = 1.8)

## accents used by Fig 2 (plateau band, the two lambda markers, degenerate region).
## `band` is a SOLID light blue rather than a transparent blue: EPS has no alpha
## channel, so an alpha fill would be dropped (or rasterized) on export.
##
## lam_plateau / lam_cv mark the two vertical lambda lines. They are soft
## "macaron" tones (a rose and a lavender), chosen to sit apart from each other
## and from the blue width curve so the two markers never read as one.
HAL_ACCENT <- c(band = "#eaf2f9", band_txt = "#2c7fb8",
                pick = "grey20", dead = "grey55",
                lam_plateau = "#d17a95",   # macaron rose      -> plateau lambda
                lam_cv      = "#8e7cc3")   # macaron lavender   -> CV-selected lambda

## data-point color. PNG uses semi-transparency; EPS (no alpha) passes a solid grey.
HAL_PTS  <- adjustcolor("grey55", 0.45)
HAL_PTS_EPS <- "grey70"

## ---- shared figure geometry ------------------------------------------------
## Both figures export at the same physical size, resolution, and base font so
## they sit side by side in the manuscript without one looking heavier.
FIG_W   <- 7.2   # inches
FIG_H   <- 5.1   # inches
FIG_DPI <- 300
FIG_PT  <- 12    # base font size (pt)

## open a PNG device at the shared geometry
hal_png <- function(file) {
  png(file, width = FIG_W, height = FIG_H, units = "in",
      res = FIG_DPI, pointsize = FIG_PT)
}

## Is the cairo EPS device actually usable?
## capabilities("cairo") is NOT trustworthy here: it can report TRUE while
## cairo_ps() still fails to load its shared library (e.g. macOS without XQuartz),
## silently leaving you with no .eps at all. So probe by really opening a device.
CAIRO_OK <- local({
  if (!isTRUE(capabilities("cairo"))) return(FALSE)
  f <- tempfile(fileext = ".eps")
  ok <- tryCatch({
    suppressWarnings(cairo_ps(f))
    opened <- !identical(names(grDevices::dev.cur()), "null device")
    if (opened) dev.off()
    opened
  }, error = function(e) FALSE)
  unlink(f)
  isTRUE(ok)
})

## open an EPS device at the shared geometry, falling back to classic postscript.
## Returns TRUE if the cairo device (which supports alpha) is in use.
hal_eps <- function(file) {
  if (CAIRO_OK) {
    cairo_ps(file, width = FIG_W, height = FIG_H,
             pointsize = FIG_PT, fallback_resolution = 600)
  } else {
    setEPS()
    postscript(file, width = FIG_W, height = FIG_H, pointsize = FIG_PT)
  }
  CAIRO_OK
}

## ggsave device name matching hal_eps()
hal_eps_device <- function() if (CAIRO_OK) cairo_ps else "eps"

## common margins + axis-title positioning; las = 1 keeps tick labels horizontal
hal_par <- function() {
  par(mar = c(4.2, 4.6, 1.2, 1), mgp = c(3, 0.7, 0), las = 1)
}

## L-shaped axes with tick marks on both axes (x and y treated the same)
hal_axes <- function(ylab = "Y", yline = 3.0) {
  axis(1)
  title(ylab = ylab, line = yline, las = 0)
}

## draw one named fit ("truth"/"lo"/"cv"/"hi") over grid xg
hal_line <- function(xg, yv, key) {
  lines(xg, yv, col = HAL_COL[[key]], lty = HAL_LTY[[key]], lwd = HAL_LWD[[key]])
}

## standard legend for the 4 elements (pass the achieved norms for lo/cv/hi).
## formatC right-aligns the norms in a fixed-width field so the parentheses line up.
hal_legend <- function(pos = "topright", norms) {
  n <- formatC(norms, format = "f", digits = 1, width = 4)
  legend(pos, bty = "n", cex = 0.85, seg.len = 2.6,
         legend = c("Truth",
                    sprintf("HAL, low norm  (%s)", n[1]),
                    sprintf("HAL, CV-chosen (%s)", n[2]),
                    sprintf("HAL, high norm (%s)", n[3])),
         col = HAL_COL[c("truth", "lo", "cv", "hi")],
         lty = HAL_LTY[c("truth", "lo", "cv", "hi")],
         lwd = HAL_LWD[c("truth", "lo", "cv", "hi")])
}

## ggplot counterpart of the base-R look, so Fig 2 matches Fig 1
hal_ggtheme <- function() {
  ggplot2::theme_classic(base_size = FIG_PT) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
}
