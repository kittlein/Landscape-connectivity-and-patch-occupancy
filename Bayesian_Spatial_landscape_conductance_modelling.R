# Supplementary Script
# Iterative Bayesian spatial GAM for genetic connectivity

library(rstan)
library(mgcv)
library(raster)
library(gdistance)
library(sf)
library(sp)

# Input files

layers_file <- "layers.tif"
layer_names_file <- "layers_names.txt"
coordinates_file <- "coordinates.txt"
connectivity_file <- "connectivity.txt"

# Model configuration

n_iterations <- 3L # change to 500...takes long
smooth_k <- 6L
epsilon <- 1e-4
n_chains <- 4L
mcmc_iterations <- 2000L
warmup_iterations <- 1000L
max_treedepth <- 12L
random_seed <- 20260928L
block_size <- 20000L

predictors <- c("changes", "sand", "elevation", "meanNDVI", "clay", "sdNDVI")

crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# Read environmental layers

layers <- raster::stack(layers_file)

layer_names <- trimws(readLines(layer_names_file, warn = FALSE))
layer_names <- layer_names[nzchar(layer_names)]

if (length(layer_names) != raster::nlayers(layers)) stop("layers_names.txt must contain one name for each raster layer.")

names(layers) <- layer_names

if (!all(c("water", "environments", predictors) %in% names(layers))) stop("Required environmental layers are missing.")

raster_crs <- raster::projection(layers)

if (is.na(raster_crs) || !nzchar(raster_crs)) {
  raster::projection(layers) <- crs_sinu
  raster_crs <- raster::projection(layers)
}

# Read locality coordinates

locality_table <- read.table(coordinates_file, header = TRUE, stringsAsFactors = FALSE)

if (!all(c("Locality", "Long", "Lat") %in% names(locality_table))) stop("coordinates.txt must contain Locality, Long and Lat.")

localities <- as.character(locality_table$Locality)

localities_sf <- sf::st_as_sf(locality_table, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

raster_coordinates <- sf::st_coordinates(sf::st_transform(localities_sf, crs = sf::st_crs(raster_crs)))

locality_xy_km <- sf::st_coordinates(sf::st_transform(localities_sf, crs = crs_sinu)) / 1000

# Read directional connectivity

connectivity <- as.matrix(read.table(connectivity_file, header = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", "")))

n_localities <- length(localities)

if (!is.numeric(connectivity) || nrow(connectivity) != n_localities || ncol(connectivity) != n_localities) stop("connectivity.txt must contain a square numeric matrix.")

if (any(!is.na(connectivity) & (connectivity < 0 | connectivity > 1))) stop("Connectivity values must be between 0 and 1.")

rownames(connectivity) <- localities
colnames(connectivity) <- localities
diag(connectivity) <- NA_real_

mean_connectivity <- mean(connectivity, na.rm = TRUE)

source_order <- rep(seq_len(n_localities), each = n_localities - 1L)
destination_order <- unlist(lapply(seq_len(n_localities), function(i) seq_len(n_localities)[-i]), use.names = FALSE)

connectivity_vector <- unlist(lapply(seq_len(n_localities), function(i) connectivity[i, -i]), use.names = FALSE)

pair_xy_km <- cbind(locality_xy_km[source_order, , drop = FALSE], locality_xy_km[destination_order, , drop = FALSE])

colnames(pair_xy_km) <- c("x_source", "y_source", "x_destination", "y_destination")

# Prepare raster values

environmental_raster <- as.data.frame(raster::getValues(layers[[predictors]]))

water <- raster::getValues(layers[["water"]])
water[is.na(water)] <- 0

environment_codes <- raster::getValues(layers[["environments"]])

environmental_raster$environments <- as.character(environment_codes)

valid_cells <- which(stats::complete.cases(environmental_raster) & rowSums(!is.finite(as.matrix(environmental_raster[predictors]))) == 0L & water != 1)

environmental_values <- environmental_raster[valid_cells, c(predictors, "environments"), drop = FALSE]

environment_levels <- as.character(sort(unique(environment_codes[!is.na(environment_codes) & water != 1])))

# Initial connectivity surface

Conn <- raster::raster(nrows = raster::nrow(layers), ncols = raster::ncol(layers), ext = raster::extent(layers), crs = raster_crs)

Conn[] <- mean_connectivity
Conn[which(water == 1)] <- 0
names(Conn) <- "Connectivity"

# Stan model

stan_code <- '
data {
  int<lower=2> N;
  int<lower=1> J;
  int<lower=1> P;
  int<lower=1> Q;
  int<lower=1> F;
  vector[N] y;
  matrix[N, J] X_linear;
  matrix[N, Q] Z_curvature;
  matrix[N, F] X_environments;
  matrix[N, N] distances_km;
  real mean_y;
  real<lower=0> sd_y;
  vector<lower=0>[J] beta_prior_sd;
  vector<lower=0>[J] tau_prior_sd;
  real<lower=0> ell_prior_km;
  real<lower=0> jitter;
}
parameters {
  real alpha;
  vector[J] beta_linear;
  vector[F] gamma_environments;
  vector[Q] z_smooth;
  vector<lower=0>[J] tau_smooth;
  real<lower=0> eta;
  real<lower=0> ell_km;
  real<lower=0> sigma;
}
transformed parameters {
  vector[Q] u_smooth;
  for (j in 1:J) {
    for (p in 1:P) u_smooth[(j - 1) * P + p] = tau_smooth[j] * z_smooth[(j - 1) * P + p];
  }
}
model {
  vector[N] mu;
  matrix[N, N] Sigma;

  mu = rep_vector(alpha, N) + X_linear * beta_linear + Z_curvature * u_smooth + X_environments * gamma_environments;

  for (i in 1:N) Sigma[i, i] = square(eta) + square(sigma) + jitter;

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      real h;
      h = sqrt(3.0) * distances_km[i, j] / ell_km;
      Sigma[i, j] = square(eta) * (1 + h) * exp(-h);
      Sigma[j, i] = Sigma[i, j];
    }
  }

  alpha ~ normal(mean_y, 2 * sd_y);
  beta_linear ~ normal(rep_vector(0, J), beta_prior_sd);
  gamma_environments ~ normal(0, 2 * sd_y);
  tau_smooth ~ normal(rep_vector(0, J), tau_prior_sd);
  z_smooth ~ std_normal();
  eta ~ normal(0, 2 * sd_y);
  sigma ~ normal(0, 2 * sd_y);
  ell_km ~ lognormal(log(ell_prior_km), 1);

  y ~ multi_normal_cholesky(mu, cholesky_decompose(Sigma));
}
'

rstan::rstan_options(auto_write = TRUE)

stan_model <- rstan::stan_model(model_code = stan_code, model_name = "iterative_connectivity_matern32_environments")

n_cores <- parallel::detectCores(logical = FALSE)

if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

# Conductance rule

mean_conductance <- function(z) {
  if (anyNA(z)) return(0)
  mean(z)
}

# Iterative fitting

metrics_table <- data.frame()

cat("\nStarting iterative analysis\n")
cat("Iterations:", n_iterations, "| k =", smooth_k, "| chains =", n_chains, "\n\n")

list.Conn <- list()
list.Paths <- list()
list.Values <- list()
list.stan.fits <- list()

for (iter in 0:(n_iterations - 1L)) {
  
  cat("iter =", iter, "| calculating paths ... ")
  
  conductance <- gdistance::transition(Conn, transitionFunction = mean_conductance, directions = 16)
  
  conductance <- gdistance::geoCorrection(conductance, type = "c")
  
  paths <- gdistance::shortestPath(conductance, origin = raster_coordinates[1L, ], goal = raster_coordinates[-1L, , drop = FALSE], output = "SpatialLines")
  
  for (i in 2:n_localities) {
    paths_i <- gdistance::shortestPath(conductance, origin = raster_coordinates[i, ], goal = raster_coordinates[-i, , drop = FALSE], output = "SpatialLines")
    paths <- rbind(paths, paths_i)
  }
  
  sp::proj4string(paths) <- raster_crs
  
  cat("environment ... ")
  
  path_data <- as.data.frame(raster::extract(layers[[predictors]], paths, df = TRUE, fun = mean, na.rm = TRUE))
  
  path_environment_values <- raster::extract(layers[["environments"]], paths, factors = FALSE)
  
  modal_environment <- rep(NA_character_, length(path_environment_values))
  
  for (i in seq_along(path_environment_values)) {
    environment_i <- as.character(path_environment_values[[i]][!is.na(path_environment_values[[i]])])
    if (length(environment_i) > 0L) {
      frequencies_i <- table(environment_i)
      modal_environment[i] <- names(frequencies_i)[which.max(frequencies_i)]
    }
  }
  
  path_data$environments <- factor(modal_environment, levels = environment_levels)
  path_data$ID <- NULL
  path_data$Connectivity <- as.numeric(connectivity_vector)
  path_data$source <- localities[source_order]
  path_data$destination <- localities[destination_order]
  
  complete_paths <- stats::complete.cases(path_data[, c("Connectivity", predictors, "environments"), drop = FALSE])
  
  path_data <- path_data[complete_paths, , drop = FALSE]
  pair_coordinates_km <- pair_xy_km[complete_paths, , drop = FALSE]
  
  if (nrow(path_data) < 20L) stop("At least 20 complete paths are required.")
  
  # GAM bases
  
  N <- nrow(path_data)
  J <- length(predictors)
  P <- smooth_k - 2L
  Q <- J * P
  F <- length(environment_levels) - 1L
  
  predictor_centers <- vapply(path_data[predictors], mean, numeric(1))
  predictor_scales <- vapply(path_data[predictors], stats::sd, numeric(1))
  
  if (any(!is.finite(predictor_scales) | predictor_scales <= 0)) stop("One or more predictors have no variation.")
  
  y <- stats::qlogis(pmin(pmax(path_data$Connectivity, epsilon), 1 - epsilon))
  sd_y <- stats::sd(y)
  
  X <- matrix(0, N, J)
  Z <- matrix(0, N, Q)
  
  smooth_bases <- vector("list", J)
  linear_index <- vector("list", J)
  curvature_index <- vector("list", J)
  
  beta_prior_sd <- numeric(J)
  tau_prior_sd <- numeric(J)
  
  for (j in seq_len(J)) {
    
    x <- (path_data[[predictors[j]]] - predictor_centers[j]) / predictor_scales[j]
    
    smooth_basis <- mgcv::smoothCon(mgcv::s(x, k = smooth_k, bs = "cr"), data = data.frame(x = x), absorb.cons = TRUE, diagonal.penalty = TRUE)[[1L]]
    
    penalty <- diag(smooth_basis$S[[1L]])
    
    linear_index[[j]] <- which(abs(penalty) < 1e-6)
    curvature_index[[j]] <- which(penalty > 1e-6)
    
    smooth_columns <- seq.int((j - 1L) * P + 1L, j * P)
    
    X[, j] <- smooth_basis$X[, linear_index[[j]]]
    Z[, smooth_columns] <- smooth_basis$X[, curvature_index[[j]], drop = FALSE]
    
    beta_prior_sd[j] <- 2 * sd_y / stats::sd(X[, j])
    tau_prior_sd[j] <- 1.5 * sd_y / sqrt(mean(rowSums(Z[, smooth_columns, drop = FALSE]^2)))
    
    smooth_bases[[j]] <- smooth_basis
  }
  
  names(smooth_bases) <- predictors
  names(linear_index) <- predictors
  names(curvature_index) <- predictors
  
  # Environments factor
  
  environment_frequencies <- table(path_data$environments)
  reference_environment <- names(environment_frequencies)[which.max(environment_frequencies)]
  
  model_environment_levels <- c(reference_environment, setdiff(environment_levels, reference_environment))
  
  X_environments <- matrix(0, nrow = N, ncol = F)
  
  for (j in seq_len(F)) {
    X_environments[, j] <- as.integer(as.character(path_data$environments) == model_environment_levels[j + 1L])
  }
  
  environments_without_paths <- setdiff(model_environment_levels, unique(as.character(path_data$environments)))
  
  # Spatial distance matrix
  
  distance_matrix <- as.matrix(stats::dist(pair_coordinates_km))
  
  positive_distances <- distance_matrix[upper.tri(distance_matrix) & distance_matrix > 0]
  
  if (!length(positive_distances)) stop("No positive distances were found among source-destination pairs.")
  
  stan_data <- list(N = N, J = J, P = P, Q = Q, F = F, y = y, X_linear = X, Z_curvature = Z, X_environments = X_environments, distances_km = distance_matrix, mean_y = mean(y), sd_y = sd_y, beta_prior_sd = beta_prior_sd, tau_prior_sd = tau_prior_sd, ell_prior_km = median(positive_distances), jitter = 1e-8 * sd_y^2)
  
  cat("Stan ... ")
  
  fit <- rstan::sampling(stan_model, data = stan_data, chains = n_chains, cores = min(n_chains, n_cores), iter = mcmc_iterations, warmup = warmup_iterations, seed = random_seed + iter, control = list(adapt_delta = 0.95, max_treedepth = max_treedepth), refresh = 0L)
  
  # MCMC diagnostics
  
  parameter_names <- c("alpha", "beta_linear", "gamma_environments", "tau_smooth", "eta", "ell_km", "sigma")
  
  parameter_summary <- summary(fit, pars = parameter_names)$summary
  
  rhats <- parameter_summary[, "Rhat"]
  max_rhat <- if (all(is.finite(rhats))) max(rhats) else NA_real_
  
  sampler_parameters <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  divergences <- sum(vapply(sampler_parameters, function(z) sum(z[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(sampler_parameters, function(z) sum(z[, "treedepth__"] >= max_treedepth), numeric(1)))
  
  diagnostics_ok <- is.finite(max_rhat) && max_rhat <= 1.01 && divergences == 0L && treedepth_hits == 0L
  
  # Posterior means
  
  posterior_samples <- rstan::extract(fit, pars = c("alpha", "beta_linear", "gamma_environments", "u_smooth", "eta", "ell_km", "sigma"), permuted = TRUE)
  
  alpha_mean <- mean(posterior_samples$alpha)
  beta_mean <- colMeans(posterior_samples$beta_linear)
  gamma_mean <- colMeans(posterior_samples$gamma_environments)
  smooth_mean <- colMeans(posterior_samples$u_smooth)
  
  # Environmental predictions
  
  environmental_linear_predictor <- as.vector(alpha_mean + X %*% beta_mean + Z %*% smooth_mean + X_environments %*% gamma_mean)
  
  fitted_environmental <- stats::plogis(environmental_linear_predictor)
  
  # Conditional spatial predictions
  
  eta_mean <- mean(posterior_samples$eta)
  ell_mean <- mean(posterior_samples$ell_km)
  sigma_mean <- mean(posterior_samples$sigma)
  
  matern_h <- sqrt(3) * distance_matrix / ell_mean
  
  spatial_covariance <- eta_mean^2 * (1 + matern_h) * exp(-matern_h)
  
  total_covariance <- spatial_covariance
  
  diag(total_covariance) <- diag(total_covariance) + sigma_mean^2 + 1e-8 * sd_y^2
  
  chol_covariance <- chol(total_covariance)
  
  spatial_effect <- as.vector(spatial_covariance %*% backsolve(chol_covariance, forwardsolve(t(chol_covariance), y - environmental_linear_predictor)))
  
  fitted_conditional <- stats::plogis(environmental_linear_predictor + spatial_effect)
  
  # Fit statistics
  
  observed_connectivity <- path_data$Connectivity
  
  total_sum_squares <- sum((observed_connectivity - mean(observed_connectivity))^2)
  
  rmse_environmental <- sqrt(mean((observed_connectivity - fitted_environmental)^2))
  
  r2_environmental <- if (total_sum_squares > 0) 1 - sum((observed_connectivity - fitted_environmental)^2) / total_sum_squares else NA_real_
  
  r2_conditional <- if (total_sum_squares > 0) 1 - sum((observed_connectivity - fitted_conditional)^2) / total_sum_squares else NA_real_
  
  # Project connectivity surface
  
  new_Conn <- raster::raster(nrows = raster::nrow(layers), ncols = raster::ncol(layers), ext = raster::extent(layers), crs = raster_crs)
  
  new_Conn[] <- NA_real_
  
  for (block_start in seq.int(1L, length(valid_cells), by = block_size)) {
    
    cell_index <- seq.int(block_start, min(block_start + block_size - 1L, length(valid_cells)))
    
    new_data <- environmental_values[cell_index, , drop = FALSE]
    
    projected_linear_predictor <- rep(alpha_mean, nrow(new_data))
    
    environment_index <- match(as.character(new_data$environments), model_environment_levels)
    
    if (anyNA(environment_index)) stop("Some environments categories were not coded in the fitted model.")
    
    non_reference_index <- which(environment_index > 1L)
    
    projected_linear_predictor[non_reference_index] <- projected_linear_predictor[non_reference_index] + gamma_mean[environment_index[non_reference_index] - 1L]
    
    for (j in seq_along(predictors)) {
      
      x_new <- (new_data[[predictors[j]]] - predictor_centers[j]) / predictor_scales[j]
      
      prediction_basis <- mgcv::PredictMat(smooth_bases[[j]], data.frame(x = x_new))
      
      smooth_columns <- seq.int((j - 1L) * P + 1L, j * P)
      
      projected_linear_predictor <- projected_linear_predictor + prediction_basis[, linear_index[[j]]] * beta_mean[j]
      
      projected_linear_predictor <- projected_linear_predictor + as.vector(prediction_basis[, curvature_index[[j]], drop = FALSE] %*% smooth_mean[smooth_columns])
    }
    
    projected_connectivity <- stats::plogis(projected_linear_predictor)
    
    new_Conn[valid_cells[cell_index]] <- pmin(pmax(projected_connectivity, epsilon), 1 - epsilon)
  }
  
  new_Conn[which(water == 1)] <- 0
  
  names(new_Conn) <- "Connectivity"
  
  # Store iteration metrics
  
  iteration_metrics <- data.frame(iter = iter, n_paths = nrow(path_data), rmse_env_train = rmse_environmental, r2_env_train = r2_environmental, r2_cond_train = r2_conditional, rhat_max = max_rhat, divergences = divergences, treedepth_hits = treedepth_hits, diagnostics_ok = diagnostics_ok, n_environments_without_paths = length(environments_without_paths))
  
  metrics_table <- rbind(metrics_table, iteration_metrics)
  
  Conn <- new_Conn

  # keep rasters, paths, values, and fits
  list.Conn[[iter+1]] <- Conn
  list.Paths[[iter+1]] <- paths
  list.Values[[iter+1]] <- path_data
  list.stan.fits[[iter+1]] <- fit
  
  cat("RMSE =", signif(rmse_environmental, 5), "| R2 env =", signif(r2_environmental, 4), "| R2 cond =", signif(r2_conditional, 4), "| Rhat =", signif(max_rhat, 4), "| divergences =", divergences, "| treedepth =", treedepth_hits, "| environments without paths =", length(environments_without_paths), "\n")
}

# Final summary

candidate_rows <- which(is.finite(metrics_table$rmse_env_train))

valid_rows <- candidate_rows[metrics_table$diagnostics_ok[candidate_rows] %in% TRUE]

if (length(valid_rows) > 0L) {
  best_row <- valid_rows[which.min(metrics_table$rmse_env_train[valid_rows])]
} else {
  best_row <- candidate_rows[which.min(metrics_table$rmse_env_train[candidate_rows])]
}

best_iteration <- metrics_table$iter[best_row]

cat("\nBest iteration =", best_iteration, "| RMSE =", 
    signif(metrics_table$rmse_env_train[best_row], 5), "| R2 env =",
    signif(metrics_table$r2_env_train[best_row], 4), "| R2 cond =", 
    signif(metrics_table$r2_cond_train[best_row], 4), "| diagnostics =", 
    metrics_table$diagnostics_ok[best_row], "\n")
