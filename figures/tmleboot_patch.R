## ============================================================
## Patch: restore observation weights in TMLEbootstrap's density HAL fit.
##
## WHY THIS EXISTS
##   hal9001::fit_hal() (0.4.6) takes `weights` as a TOP-LEVEL argument, has no
##   `...`, and its body runs, unconditionally:
##
##       fit_control$weights <- weights          # `weights` defaults to NULL
##
##   So any weights placed *inside* fit_control are overwritten with NULL and
##   never reach glmnet.
##
##   Neither published version of TMLEbootstrap works against hal9001 0.4.6:
##
##     * upstream (wilsoncai1992) passes `weights` at top level -- correct -- but
##       also passes fit_type / use_min / cv_select there. Those formals no longer
##       exist, and with no `...` to absorb them fit_hal errors outright:
##           unused arguments (fit_type = "glmnet", use_min = TRUE, cv_select = FALSE)
##
##     * the fork (zbutzin, branch fix-hal9001-weights-and-cleanup) fixes that by
##       moving those arguments into fit_control -- but sweeps `weights` in with
##       them. It runs, and silently drops the weights.
##
##   The density is fit on frequency-compressed bin rows, so dropping the Freq
##   weights makes every bin count once regardless of how many observations fall
##   in it. The resulting density is nearly constant, does not resolve the modes,
##   and is insensitive to lambda. Because the Wald interval is influence-curve
##   based -- EIC = 2 * (p_hat - Psi), SE = sd(EIC) / sqrt(n) -- a flat p_hat
##   collapses the CI width to near zero. That is the artifact this patch removes.
##
##   Sanity check on simulate_density_data(n_sim = 2000, n_mode = 3), whose truth
##   is an equal-weight Gaussian mixture with modes at -4/0/4, so the target has a
##   closed form: Psi = Integral p(x)^2 dx = 0.16926.
##
##       weights dropped (fork as installed) : Psi = 0.086,  plateau width 3.6e-4
##       weights restored (this patch)       : Psi = 0.165,  plateau width 1.5e-2
##
## WHAT IT DOES
##   Keeps the fork's fit_control migration (cv_select belongs there) and moves
##   `weights` back to the top level, where fit_hal expects it. That is the only
##   change. Everything else in densityHAL$fit is byte-for-byte the fork's.
##
##   Remove this file once the fork carries the fix; `verify_weights_patch()`
##   below will tell you loudly if it ever stops being needed or stops working.
## ============================================================

stopifnot(requireNamespace("TMLEbootstrap", quietly = TRUE),
          requireNamespace("hal9001", quietly = TRUE))

local({
  ns  <- asNamespace("TMLEbootstrap")
  gen <- get("densityHAL", envir = ns)

  ## Guard: this patch is only meaningful while fit_hal still has the top-level
  ## `weights` formal that overwrites fit_control$weights. If hal9001 ever changes
  ## that contract, stop rather than silently fitting the wrong thing.
  fm <- names(formals(hal9001::fit_hal))
  if (!("weights" %in% fm)) {
    stop("hal9001::fit_hal has no top-level `weights` argument. ",
         "This patch assumes it does -- re-check before trusting Figure 2.")
  }

  gen$public_methods$fit <- function(lambda = 2e-05, ...) {
    df_compressed <- self$longiData$generate_df_compress(x = self$x)
    fit_control_dens <- list(cv_select = FALSE)
    if (is.null(lambda)) fit_control_dens$cv_select <- TRUE
    self$hal_fit <- hal9001::fit_hal(
      X            = df_compressed[, "box"],
      Y            = df_compressed$Y,
      family       = "binomial",
      lambda       = lambda,
      weights      = df_compressed$Freq,   # <-- TOP LEVEL, not inside fit_control
      fit_control  = fit_control_dens,
      return_lasso = TRUE,
      return_x_basis = FALSE,
      yolo         = FALSE
    )
  }
})

## Confirm the weights actually reach glmnet, using a fit small enough to be cheap
## but weighted strongly enough that ignoring the weights cannot coincidentally
## produce the same coefficients.
verify_weights_patch <- function() {
  set.seed(7)
  X <- as.matrix(runif(300))
  Y <- rbinom(300, 1, 0.5)
  w <- rep(c(1, 50), 150)
  cf <- function(...) {
    as.numeric(hal9001::fit_hal(X = X, Y = Y, family = "binomial", lambda = 1e-4,
                                return_lasso = TRUE, yolo = FALSE, ...)$lasso_fit$beta)
  }
  unweighted <- cf(fit_control = list(cv_select = FALSE))
  top_level  <- cf(fit_control = list(cv_select = FALSE), weights = w)
  if (isTRUE(all.equal(unweighted, top_level))) {
    stop("Top-level `weights` had no effect on hal9001::fit_hal -- the density fit ",
         "would be unweighted and Figure 2's CI widths would be meaningless.")
  }
  invisible(TRUE)
}
verify_weights_patch()

message("TMLEbootstrap: densityHAL$fit patched -- observation weights restored.")
