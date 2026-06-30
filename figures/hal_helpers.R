# Dependencies ------------------------------------------------------------
library(glmnet)


# Data Generation ---------------------------------------------------------

#' Generate non-linear synthetic data
#' 
#' Creates univariate data following y = 2 * sin(π/2 * |x|) + ε
#' @param n Sample size
#' @param x_min Minimum x value (default -4)
#' @param x_max Maximum x value (default 4)
#' @param seed Random seed for reproducibility
#' @return Data frame with x1 and y
generate_nonlinear_data <- function(n, x_min = -4, x_max = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  x1 <- runif(n, x_min, x_max)
  y <- 2 * sin(pi / 2 * abs(x1)) + rnorm(n)
  
  return(data.frame(x1 = x1, y = y))
}


# HAL Design Matrix Construction -----------------------------------------

#' Build HAL design matrix with indicator basis functions
#' 
#' Creates design matrix where entry [i, j] = I(x_i >= x_j)
#' This forms the basis for univariate HAL regression
#' @param x_train Training covariate vector
#' @param x_new Optional new data points for prediction
#' @return Design matrix (n x n for training, or m x n for new data)
build_hal_design_matrix <- function(x_train, x_new = NULL) {
  if (is.null(x_new)) {
    x_new <- x_train
  }
  
  design <- outer(x_new, x_train, FUN = function(a, b) as.numeric(a >= b))
  return(design)
}


# Model Fitting -----------------------------------------------------------

#' Fit HAL regression using cross-validated LASSO
#' 
#' Fits HAL model with automatic lambda selection via CV
#' @param design_matrix HAL design matrix
#' @param y Response vector
#' @param nlambda Number of lambda values to try (default 100)
#' @param nfolds Number of CV folds (default 10)
#' @return cv.glmnet fitted object
fit_hal_cv <- function(design_matrix, y, nlambda = 100, nfolds = 10) {
  cv_fit <- cv.glmnet(
    x = design_matrix,
    y = y,
    family = "gaussian",
    alpha = 1,
    nlambda = nlambda,
    nfolds = nfolds,
    lambda.min.ratio = 1e-3
  )
  
  return(cv_fit)
}


#' Calculate variation norm for all lambda values
#' 
#' Variation norm is the L1 norm of coefficients (excluding intercept)
#' @param cv_fit Fitted cv.glmnet object
#' @return Vector of variation norms for each lambda
calculate_variation_norms <- function(cv_fit) {
  glmnet_fit <- cv_fit$glmnet.fit
  all_lambdas <- cv_fit$lambda
  
  # Extract all coefficients
  coef_matrix <- as.matrix(coef(glmnet_fit, s = all_lambdas))
  
  # Calculate L1 norm excluding intercept
  variation_norms <- colSums(abs(coef_matrix[-1, , drop = FALSE]))
  
  return(variation_norms)
}


#' Find lambda values corresponding to target variation norms
#' 
#' Identifies lambda values that produce variation norms closest to targets
#' @param cv_fit Fitted cv.glmnet object
#' @param target_norms Vector of target variation norms
#' @return List with lambda values and achieved norms
find_lambdas_by_norm <- function(cv_fit, target_norms) {
  variation_norms <- calculate_variation_norms(cv_fit)
  all_lambdas <- cv_fit$lambda
  
  results <- list()
  
  for (target in target_norms) {
    idx <- which.min(abs(variation_norms - target))
    results[[paste0("norm_", target)]] <- list(
      lambda = all_lambdas[idx],
      achieved_norm = variation_norms[idx],
      target_norm = target
    )
  }
  
  return(results)
}


# Prediction --------------------------------------------------------------

#' Make predictions for multiple lambda values
#' 
#' @param glmnet_fit Fitted glmnet object
#' @param design_new Design matrix for new data
#' @param lambdas Vector of lambda values for prediction
#' @return Matrix of predictions (rows = observations, cols = lambdas)
predict_hal_multiple <- function(glmnet_fit, design_new, lambdas) {
  predictions <- sapply(lambdas, function(lam) {
    predict(glmnet_fit, newx = design_new, s = lam, type = "response")
  })
  
  return(predictions)
}


# Visualization -----------------------------------------------------------

#' Plot HAL fits with different regularization levels
#' 
#' @param x_train Training x values
#' @param y_train Training y values
#' @param x_grid Grid of x values for plotting
#' @param predictions Matrix of predictions for each regularization level
#' @param labels Character vector of labels for each fit
#' @param main Plot title
#' @param line_types Line types for each fit
#' @param line_widths Line widths for each fit
plot_hal_fits <- function(x_train, y_train, x_grid, predictions, labels,
                         main = "HAL Regression with Different Regularization",
                         line_types = c(1, 2, 3),
                         line_widths = c(3, 2, 2)) {
  
  # Plot training data
  plot(x_train, y_train,
       pch = 16, col = "grey70",
       xlab = "x", ylab = "y / fitted f(x)",
       main = main)
  
  # Add fitted lines
  for (i in seq_len(ncol(predictions))) {
    lines(x_grid, predictions[, i], 
          lwd = line_widths[i], 
          lty = line_types[i])
  }
  
  # Add legend
  legend("topleft",
         legend = c("Observations", labels),
         col = c("grey70", rep("black", length(labels))),
         pch = c(16, rep(NA, length(labels))),
         lwd = c(NA, line_widths),
         lty = c(NA, line_types),
         bty = "n")
}


# Main Execution ----------------------------------------------------------

#' Run complete HAL simulation example
#' 
#' Demonstrates HAL regression with different regularization levels
#' @param n_obs Number of observations (default 500)
#' @param target_norms Target variation norms to compare (default c(4.8, 35.2))
#' @param seed Random seed (default 123)
run_hal_simulation <- function(n_obs = 500, target_norms = c(4.8, 35.2), seed = 123) {
  
  # 1. Generate data
  cat("Generating data...\n")
  data <- generate_nonlinear_data(n_obs, seed = seed)
  x_train <- data$x1
  y_train <- data$y
  
  # 2. Build HAL design matrix
  cat("Building HAL design matrix...\n")
  design_train <- build_hal_design_matrix(x_train)
  
  # 3. Fit HAL with cross-validation
  cat("Fitting HAL model with CV...\n")
  cv_fit <- fit_hal_cv(design_train, y_train)
  
  # 4. Calculate variation norms and find target lambdas
  cat("Finding lambdas for target variation norms...\n")
  lambda_cv <- cv_fit$lambda.min
  lambda_results <- find_lambdas_by_norm(cv_fit, target_norms)
  
  # Print results
  cat("\n=== Lambda Selection Results ===\n")
  cat("CV-selected lambda:", lambda_cv, "\n\n")
  
  for (result_name in names(lambda_results)) {
    result <- lambda_results[[result_name]]
    cat(sprintf("Target norm %.1f: lambda = %.6f, achieved norm = %.2f\n",
                result$target_norm, result$lambda, result$achieved_norm))
  }
  
  # 5. Make predictions on grid
  cat("\nGenerating predictions...\n")
  x_grid <- seq(-4, 4, length.out = 400)
  design_new <- build_hal_design_matrix(x_train, x_grid)
  
  lambdas_to_plot <- c(lambda_cv, sapply(lambda_results, function(x) x$lambda))
  predictions <- predict_hal_multiple(cv_fit$glmnet.fit, design_new, lambdas_to_plot)
  
  # 6. Create visualization
  cat("Creating plot...\n")
  labels <- c("CV-chosen", 
              paste0("≈ ", target_norms[1]), 
              paste0("≈ ", target_norms[2]))
  
  plot_hal_fits(
    x_train, y_train, x_grid, predictions, labels,
    main = sprintf("HAL (univariate indicators): CV, ~%.1f, ~%.1f", 
                  target_norms[1], target_norms[2])
  )
  
  # Return results
  invisible(list(
    data = data,
    cv_fit = cv_fit,
    lambda_results = lambda_results,
    predictions = predictions,
    x_grid = x_grid
  ))
}


# Execute simulation ------------------------------------------------------
if (sys.nframe() == 0) {
  # Only run if script is executed directly (not sourced)
  results <- run_hal_simulation(n_obs = 500, target_norms = c(4.8, 35.2), seed = 123)
}
