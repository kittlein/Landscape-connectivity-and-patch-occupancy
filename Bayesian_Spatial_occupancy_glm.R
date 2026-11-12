# Bayesian spatial logistic model for patch occupancy

library(rstan)
library(openxlsx)
library(sf)
library(sp)
library(raster)

# Input

occupancy_file <- "patch_occupancy.xlsx"

n_chains <- 4L
n_iterations <- 4000L
n_warmup <- 2000L
random_seed <- 20261002L
max_treedepth <- 12L
ell_prior_median <- 50

crs_ll <- "+proj=longlat +datum=WGS84 +no_defs"
crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# Read patch data

patches <- openxlsx::read.xlsx(occupancy_file)

required_columns <- c("lon", "lat", "area", "Occupancy")

if (!all(required_columns %in% names(patches))) stop("patch_occupancy.xlsx must contain lon, lat, area and Occupancy.")

best_Conn <- raster::raster("best_iter_conductance.tif")

# Prepare patch coordinates

patches_sf <- sf::st_as_sf(patches, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

patches_sf <- sf::st_transform(patches_sf, crs = crs_sinu)

# Create patch buffers

patch_buffers <- sf::st_buffer(patches_sf, dist = 565 * sqrt(patches_sf$area))

# Extract connectivity

connectivity_extraction <- raster::extract(best_Conn, as(patch_buffers, "Spatial"), fun = mean, df = TRUE, na.rm = TRUE)

patches$Connectivity <- as.numeric(connectivity_extraction[[2L]])

# Prepare occupancy

if (is.factor(patches$Occupancy)) patches$Occupancy <- as.character(patches$Occupancy)

if (is.character(patches$Occupancy)) {
  
  occupancy_labels <- tolower(trimws(patches$Occupancy))
  
  patches$Occupancy <- ifelse(occupancy_labels %in% c("1", "true", "yes", "occupied", "present"), 1L,
                              ifelse(occupancy_labels %in% c("0", "false", "no", "unoccupied", "absent", "empty"), 0L, NA_integer_))
}

# Patch coordinates in km

patch_xy <- sf::st_coordinates(patches_sf)

patches$x_km <- patch_xy[, 1L] / 1000
patches$y_km <- patch_xy[, 2L] / 1000

# Remove incomplete observations

patches <- patches[complete.cases(patches[, c("Connectivity", "Occupancy", "x_km", "y_km")]), ]

if (!all(patches$Occupancy %in% 0:1)) stop("Occupancy must be coded as 0 or 1.")

if (length(unique(patches$Occupancy)) != 2L) stop("Both occupancy states must be present.")

if (sd(patches$Connectivity) == 0) stop("Connectivity has no variation among patches.")

# Scale connectivity

mean_connectivity <- mean(patches$Connectivity)

patches$connectivity_01 <- (patches$Connectivity - mean_connectivity) / 0.1

# Spatial distances

distance_matrix <- as.matrix(dist(patches[, c("x_km", "y_km")]))

# Stan data

stan_data <- list(N = nrow(patches),
                  y = as.integer(patches$Occupancy),
                  x = as.numeric(patches$connectivity_01),
                  D = distance_matrix,
                  ell_prior_median = ell_prior_median)

# Stan model

stan_code <- '
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> y;
  vector[N] x;
  matrix<lower=0>[N,N] D;
  real<lower=0> ell_prior_median;
}

parameters {
  real alpha;
  real beta_01;
  real<lower=0> eta;
  real<lower=0> ell_km;
  vector[N] z;
}

transformed parameters {
  vector[N] u;
  {
    matrix[N,N] K;

    for (i in 1:N) {

      K[i,i] = square(eta) + 1e-8;

      if (i < N) {

        for (j in (i + 1):N) {

          real h;

          h = sqrt(3.0) * D[i,j] / ell_km;

          K[i,j] = square(eta) * (1.0 + h) * exp(-h);

          K[j,i] = K[i,j];
        }
      }
    }

    u = cholesky_decompose(K) * z;
  }
}

model {
  alpha ~ normal(0,2);
  beta_01 ~ normal(0,1);
  eta ~ normal(0,1);
  ell_km ~ lognormal(log(ell_prior_median),1);
  z ~ std_normal();

  y ~ bernoulli_logit(alpha + beta_01 * x + u);
}

generated quantities {
  real OR_01;
  vector[N] p_conditional;
  vector[N] p_spatial_zero;

  OR_01 = exp(beta_01);

  p_conditional = inv_logit(alpha + beta_01 * x + u);

  p_spatial_zero = inv_logit(alpha + beta_01 * x);
}
'

# Compile model

rstan::rstan_options(auto_write = TRUE)

n_cores <- parallel::detectCores(logical = FALSE)

if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

options(mc.cores = min(n_chains, n_cores))

stan_model <- rstan::stan_model(model_code = stan_code, model_name = "patch_occupancy_matern32")

# Fit model

stan_fit <- rstan::sampling(stan_model,
                            data = stan_data,
                            chains = n_chains,
                            iter = n_iterations,
                            warmup = n_warmup,
                            seed = random_seed,
                            refresh = 100,
                            control = list(adapt_delta = 0.95, max_treedepth = max_treedepth))

# Extract posterior samples

posterior_samples <- rstan::extract(stan_fit,
                                    pars = c("alpha", "beta_01", "eta", "ell_km", "OR_01", "p_conditional", "p_spatial_zero"),
                                    permuted = TRUE)

# Posterior summary

beta_ci <- quantile(posterior_samples$beta_01, c(0.025, 0.5, 0.975))
or_ci <- quantile(posterior_samples$OR_01, c(0.025, 0.5, 0.975))
eta_ci <- quantile(posterior_samples$eta, c(0.025, 0.5, 0.975))
ell_ci <- quantile(posterior_samples$ell_km, c(0.025, 0.5, 0.975))

posterior_summary <- data.frame(parameter = c("beta_per_0.1", "OR_per_0.1", "spatial_eta", "ell_km"),
                                lower_95 = c(beta_ci[1], or_ci[1], eta_ci[1], ell_ci[1]),
                                median = c(beta_ci[2], or_ci[2], eta_ci[2], ell_ci[2]),
                                upper_95 = c(beta_ci[3], or_ci[3], eta_ci[3], ell_ci[3]))

# MCMC diagnostics

diagnostics <- summary(stan_fit, pars = c("alpha", "beta_01", "eta", "ell_km"))$summary

sampler_parameters <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)

divergences <- sum(vapply(sampler_parameters, function(x) sum(x[, "divergent__"]), numeric(1L)))

treedepth_hits <- sum(vapply(sampler_parameters, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1L)))

# Patch posterior probabilities

patches$prob_conditional <- apply(posterior_samples$p_conditional, 2L, median)

patches$prob_conditional_lower <- apply(posterior_samples$p_conditional, 2L, quantile, probs = 0.025)

patches$prob_conditional_upper <- apply(posterior_samples$p_conditional, 2L, quantile, probs = 0.975)

patches$prob_spatial_zero <- apply(posterior_samples$p_spatial_zero, 2L, median)

# Console output

cat("\nBayesian spatial logistic model\n\n")

print(posterior_summary, row.names = FALSE, digits = 4)

cat("\nP(beta > 0 | data) =", round(mean(posterior_samples$beta_01 > 0), 4), "\n")

cat("\nMCMC diagnostics\n\n")

print(diagnostics[, c("n_eff", "Rhat"), drop = FALSE], digits = 4)

cat("\nDivergences =", divergences, "\n")

cat("Maximum treedepth hits =", treedepth_hits, "\n")

par(mar = c(4.5, 4.8, 2.5, 1), mgp=c(3,0.8,0), mfrow=c(1,1))
hist(posterior_samples$beta_01, n=100, las=1, col="transparent", main=NULL, 
     xlab="", cex.lab=2, cex.axis=1.2)
mtext("Connectivity", side=1, cex=2, line=2)

ql=quantile(posterior_samples$beta_01, 0.05)
segments(ql,0, min(ql), 500, lty=3, lwd=3)
segments(max(ql),0, max(ql), 500, lty=3, lwd=3)

sum(posterior_samples$beta_01<0)/length(posterior_samples$beta_01)
