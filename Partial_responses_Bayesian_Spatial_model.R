# Select best iteration

best_iter <- readRDS("best_iter_conductance_surface.rds")

best_Conn <- best_iter$Conn
best_Paths <- best_iter$Lineas
best_Values <- best_iter$datos
best_stan_fit <- best_iter$ajuste

# Model settings

predictors <- c("changes", "clay", "elevation", "meanNDVI", "sand", "sdNDVI")
smooth_k <- 6L

n_predictors <- length(predictors)
n_curvature <- smooth_k - 2L

# Reconstruct GAM bases

predictor_centers <- vapply(best_Values[predictors], mean, na.rm=TRUE, numeric(1))
predictor_scales <- vapply(best_Values[predictors], stats::sd, numeric(1))

smooth_bases <- vector("list", n_predictors)
linear_index <- vector("list", n_predictors)
curvature_index <- vector("list", n_predictors)

for (j in seq_along(predictors)) {
  
  x <- (best_Values[[predictors[j]]] - predictor_centers[j]) / predictor_scales[j]
  
  smooth_basis <- mgcv::smoothCon(mgcv::s(x, k = smooth_k, bs = "cr"), data = data.frame(x = x), absorb.cons = TRUE, diagonal.penalty = TRUE)[[1L]]
  
  penalty <- diag(smooth_basis$S[[1L]])
  
  linear_index[[j]] <- which(abs(penalty) < 1e-6)
  curvature_index[[j]] <- which(penalty > 1e-6)
  smooth_bases[[j]] <- smooth_basis
}

names(smooth_bases) <- predictors
names(linear_index) <- predictors
names(curvature_index) <- predictors

# Extract posterior samples

posterior_samples <- rstan::extract(best_stan_fit, pars = c("beta_linear", "u_smooth"), permuted = TRUE)

# Calculate partial responses

partial_curves <- vector("list", n_predictors)
partial_importance <- numeric(n_predictors)

names(partial_curves) <- predictors
names(partial_importance) <- predictors

for (j in seq_along(predictors)) {
  
  smooth_columns <- seq.int((j - 1L) * n_curvature + 1L, j * n_curvature)
  
  x_grid <- seq(min(best_Values[[predictors[j]]]), max(best_Values[[predictors[j]]]), length.out = 100L)
  
  x_standardized <- (x_grid - predictor_centers[j]) / predictor_scales[j]
  
  prediction_basis <- mgcv::PredictMat(smooth_bases[[j]], data.frame(x = x_standardized))
  
  partial_effect <- outer(posterior_samples$beta_linear[, j], prediction_basis[, linear_index[[j]]]) +
    posterior_samples$u_smooth[, smooth_columns, drop = FALSE] %*% t(prediction_basis[, curvature_index[[j]], drop = FALSE])
  
  observed_basis <- smooth_bases[[j]]$X
  
  mean_effect <- as.numeric(posterior_samples$beta_linear[, j] * colMeans(observed_basis)[linear_index[[j]]] +
                              posterior_samples$u_smooth[, smooth_columns, drop = FALSE] %*% colMeans(observed_basis)[curvature_index[[j]]])
  
  partial_effect <- sweep(partial_effect, 1L, mean_effect, "-")
  
  credible_bands <- apply(partial_effect, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  
  partial_importance[j] <- stats::sd(credible_bands[2L, ])
  
  partial_curves[[j]] <- list(x = x_grid, lower = credible_bands[1L, ], median = credible_bands[2L, ], upper = credible_bands[3L, ])
}

# Order predictors by partial-response variation

plot_order <- order(partial_importance, decreasing = TRUE)

importance_table <- data.frame(predictor = predictors[plot_order], partial_response_sd = partial_importance[plot_order])

print(importance_table, row.names = FALSE, digits = 3)

# Plot ordered partial responses

y_limits <- c(-4, 4)

par(mfrow = c(2, 3), mar = c(2, 4.4, 2.4, 1), oma = c(4, 3, 0, 0), mgp = c(2, 0.75, 0), las = 1)

for (j in plot_order) {
  
  curve <- partial_curves[[j]]
  
  plot(curve$x, curve$median, type = "n", xlab = "", ylab = "", 
       main = predictors[j], ylim = y_limits, cex.axis = 1.5, cex.main = 2)
  
  polygon(c(curve$x, rev(curve$x)), c(curve$lower, rev(curve$upper)), col = "grey80", border = NA)
  
  lines(curve$x, curve$median, col = "black", lwd = 2.5)
  
  abline(h = 0, lty = 2, col = "grey45")
}

mtext("Partial effect", side = 2, outer = TRUE, line = -1, cex = 3, las = 3)

mtext("Predictor value", side = 1, outer = TRUE, line = 2, cex = 3)