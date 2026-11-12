# Landscape connectivity and patch occupancy

Supplementary data and R code for *Linking genetic connectivity and habitat occupancy through environmental resistance in fragmented landscapes*.

## Overview

The study examines genetic connectivity among sampled localities of *Ctenomys mendocinus* and asks whether a landscape connectivity surface derived from genetic and environmental information is associated with patch occupancy. The analytical sequence begins with microsatellite genotypes, estimates directional connectivity among localities, relates those values to geographic distance and environmental conditions along least-cost paths, and finally evaluates occupancy using field observations obtained independently of the genetic sampling.

The spatial landscape analysis uses an iterative Bayesian Gaussian additive model for logit-transformed genetic connectivity. The subsequent patch-occupancy analysis uses a separate Bayesian spatial binomial model. These analyses answer different questions: the first describes variation in connectivity among directed locality pairs; the second estimates the association between predicted patch connectivity and the probability of occupancy.

## Files and data requirements

This repository contains **six R scripts**, this README, and `README.html`. The input files `microsatellites.xlsx`, `coordinates.txt`, `layers.tif`, [`layers_names.txt`](layers_names.txt), and `patch_occupancy.xlsx` are included. The repository also includes the intermediate files `microsatellites_genepop.gen`, `population_codes.csv`, and `connectivity.txt`, and the selected model outputs `best_iter_conductance.tif` and `best_iter_conductance_surface.rds`.

### Genetic data and localities

`microsatellites.xlsx` is expected to contain a sheet named `microsatellites` with multilocus genotypes. The first script converts those records to GENEPOP format. The next script uses `diveRsity::divMigrate()` to derive a matrix of relative, directional genetic connectivity. The landscape model reads its response values from the supplied `connectivity.txt`. These values quantify relative directional connectivity.

`coordinates.txt` supplies the sampled locality names and longitude/latitude coordinates. The landscape script requires columns named `Locality`, `Long`, and `Lat`; the matrix in `connectivity.txt` must use the same locality ordering. Coordinates link each directional source-destination pair to its route through the environmental rasters and to the spatial covariance model.

A geographic-distance analysis calculates straight-line distances between localities after transformation to the sinusoidal projection used in the scripts (`+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs`). Because connectivity is directional, the Mantel script considers both triangles of its matrix, calculates an observed correlation across both directions, and uses permutation outputs to calculate an empirical probability. The Mantel result concerns geographic distance; it is separate from the environmental and occupancy models.

### Environmental layers

`layers.tif` is expected to contain eight raster layers, in the order specified by `layers_names.txt`:

| Layer | Role in these scripts |
| --- | --- |
| `water` | Mask for cells assigned zero conductance; not a fitted predictor. |
| `environments` | Categorical environmental or habitat classes in the landscape model. |
| `changes` | Continuous land-cover-change predictor. |
| `clay` | Continuous soil-clay predictor. |
| `elevation` | Continuous elevation predictor. |
| `meanNDVI` | Continuous mean vegetation-index predictor. |
| `sand` | Continuous soil-sand predictor. |
| `sdNDVI` | Continuous vegetation-index variability predictor. |

The landscape script extracts summaries of these environmental conditions along least-cost routes. It transforms locality coordinates to the raster coordinate system for route calculations, and uses source and destination coordinates in a sinusoidal metric projection for the spatial covariance. The raster's coordinate reference system must correspond to its actual spatial coordinates.

### Patch occupancy

`patch_occupancy.xlsx` is expected to contain `lon`, `lat`, `area`, and `Occupancy`. The occupancy script creates a buffer for each surveyed patch, extracts the mean value of a selected connectivity raster within that buffer, and converts recognized occupancy labels such as `occupied` and `empty` to 1 and 0. Field occupancy observations are independent of the genetic data used to build the landscape surface. The fitted coefficient therefore describes an association between predicted connectivity and observed occupancy; it does not, by itself, establish a causal effect.

## Analytical workflow and scripts

### 1. Convert microsatellite genotypes

The conversion script reads the `microsatellites` sheet of `microsatellites.xlsx`, formats the locus and population records, and writes `microsatellites_genepop.gen`. This file is required by the next step. [Open the R script](convert_microsatellites_to_genepop.R).

### 2. Estimate directional genetic connectivity

The connectivity script applies `divMigrate()` to the GENEPOP file, displays a network using `qgraph`, and writes the `gRelMig` matrix to `connectivity.txt`. Running this step overwrites the supplied matrix. The order of localities in the coordinate file and matrix must agree for subsequent spatial analyses. [Open the R script](estimate_asymmetric_connectivity.R).

### 3. Evaluate connectivity against geographic distance

The Mantel script reads `coordinates.txt` and `connectivity.txt`, calculates pairwise straight-line distances in kilometres, and evaluates the two directional triangles of the connectivity matrix using 19,999 permutations for each Mantel call. It uses `rstudioapi` to set its working directory from the active RStudio document; adjust that line if running outside RStudio. [Open the R script](mantel_connectivity_distance.R).

### 4. Model the landscape connectivity surface

The landscape script starts from a constant conductance surface based on mean observed connectivity and sets water cells to zero. At each iteration it obtains least-cost paths between directed locality pairs, summarizes environmental values along those paths, and fits a Gaussian Bayesian additive model to logit-transformed connectivity. The six continuous predictors have penalized cubic regression splines (`k = 6`); `environments` enters as a factor. A Matérn 3/2 covariance represents residual spatial correlation using the coordinates of both the source and destination of each pair.

Posterior environmental effects are projected back to raster cells to update the surface for the next iteration. The raster projection uses the environmental component of the fitted model, without adding pair-specific spatial residuals. The script reports training `RMSE_env`, `R2_env`, and `R2_cond` and selects an iteration by the lowest environmental training RMSE among iterations with acceptable diagnostics when available. These are training metrics, not independent predictive validation. The supplied setting is `n_iterations <- 3L`; change it deliberately if a longer iterative run is needed. [Open the R script](Bayesian_Spatial_landscape_conductance_modelling.R).

### 5. Plot environmental partial responses

The partial-response script reconstructs spline bases, extracts posterior samples from a selected landscape fit, and plots centered curves with pointwise 95% credible intervals for the six continuous predictors. Its ordering statistic is the variation of each posterior median curve across the observed predictor range; it is descriptive and is not a test of global significance. It reads the included `best_iter_conductance_surface.rds` and uses its `Conn`, `Lineas`, `datos`, and `ajuste` components. Before using the curves, ensure that the predictor order used to reconstruct spline bases matches the order used in the fitted landscape model. [Open the R script](Partial_responses_Bayesian_Spatial_model.R).

### 6. Relate predicted connectivity to patch occupancy

The occupancy script reads the independent patch survey, extracts mean predicted connectivity around each patch from a selected raster, and fits a Bernoulli logit model. A Matérn 3/2 spatial effect based on patch coordinates accounts for spatial correlation not explained by connectivity. The script loads the included `best_iter_conductance.tif` as `best_Conn`. Its posterior coefficient and fitted occupancy curve refer to the particular surface supplied. [Open the R script](Bayesian_Spatial_occupancy_glm.R).

## Reproducibility and linked files

Install the required R packages and run the scripts from the repository folder. Use the supplied `connectivity.txt` for the distance and landscape analyses. The partial-response and occupancy scripts read the included RDS and GeoTIFF, respectively. To regenerate the genetic intermediates, run steps 1 and 2 in a separate working copy. After refitting the landscape model, the fit and iteration objects remain in the R session. To use a newly selected iteration in steps 5 and 6, save the fit and surface to the RDS and GeoTIFF read by those scripts.
