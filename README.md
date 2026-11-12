# Landscape connectivity and patch occupancy

Supplementary R code for *Linking genetic connectivity and habitat occupancy through environmental resistance in fragmented landscapes*.

## Overview

The study examines genetic connectivity among sampled localities of *Ctenomys mendocinus* and asks whether a landscape connectivity surface derived from genetic and environmental information is associated with patch occupancy. The analytical sequence begins with microsatellite genotypes, estimates directional connectivity among localities, relates those values to geographic distance and environmental conditions along least-cost paths, and finally evaluates occupancy using field observations obtained independently of the genetic sampling.

The spatial landscape analysis uses an iterative Bayesian Gaussian additive model for logit-transformed genetic connectivity. The subsequent patch-occupancy analysis uses a separate Bayesian spatial binomial model. These analyses answer different questions: the first describes variation in connectivity among directed locality pairs; the second estimates the association between predicted patch connectivity and the probability of occupancy.

## Files and data requirements

This ZIP contains **six R scripts**, six identical plain-text copies (`.R.txt`) for reading them without opening an R editor, [`layers_names.txt`](layers_names.txt), this README, and `README.pdf`. It does **not** contain the input datasets `microsatellites.xlsx`, `coordinates.txt`, `layers.tif`, or `patch_occupancy.xlsx`. These inputs must be obtained separately and placed where the scripts expect them, or their paths must be edited. The scripts generate `microsatellites_genepop.gen` and `connectivity.txt` in the first two steps.

### Genetic data and localities

`microsatellites.xlsx` is expected to contain a sheet named `microsatellites` with multilocus genotypes. The first script converts those records to GENEPOP format. The next script uses `diveRsity::divMigrate()` to derive a matrix of relative, directional genetic connectivity. The matrix supplies the response values for the landscape model; these estimates should be interpreted as relative connectivity rather than absolute migration rates.

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

The connectivity script applies `divMigrate()` to the GENEPOP file, displays a network using `qgraph`, and writes the directional matrix to `connectivity.txt`. The order of localities in the coordinate file and matrix must agree for subsequent spatial analyses. [Open the R script](estimate_asymmetric_connectivity.R).

### 3. Evaluate connectivity against geographic distance

The Mantel script reads `coordinates.txt` and `connectivity.txt`, calculates pairwise straight-line distances in kilometres, and evaluates the two directional triangles of the connectivity matrix using 19,999 permutations for each Mantel call. It uses `rstudioapi` to set its working directory from the active RStudio document; adjust that line if running outside RStudio. [Open the R script](mantel_connectivity_distance.R).

### 4. Model the landscape connectivity surface

The landscape script starts from a constant conductance surface based on mean observed connectivity and sets water cells to zero. At each iteration it obtains least-cost paths between directed locality pairs, summarizes environmental values along those paths, and fits a Gaussian Bayesian additive model to logit-transformed connectivity. The six continuous predictors have penalized cubic regression splines (`k = 6`); `environments` enters as a factor. A Matérn 3/2 covariance represents residual spatial correlation using the coordinates of both the source and destination of each pair.

Posterior environmental effects are projected back to raster cells to update the surface for the next iteration. The raster projection uses the environmental component of the fitted model, without adding pair-specific spatial residuals. The script reports training `RMSE_env`, `R2_env`, and `R2_cond` and selects an iteration by the lowest environmental training RMSE among iterations with acceptable diagnostics when available. These are training metrics, not independent predictive validation. The supplied setting is `n_iterations <- 3L`; change it deliberately if a longer iterative run is needed. [Open the R script](Bayesian_Spatial_landscape_conductance_modelling.R).

### 5. Plot environmental partial responses

The partial-response script reconstructs spline bases, extracts posterior samples from a selected landscape fit, and plots centered curves with pointwise 95% credible intervals for the six continuous predictors. Its ordering statistic is the variation of each posterior median curve across the observed predictor range; it is descriptive and is not a test of global significance. It requires a separately created `best_iter_conductance_surface.rds` with components `Conn`, `Lineas`, `datos`, and `ajuste`. That RDS is not included in this ZIP and the landscape script does not automatically write it. Before using the curves, ensure that the predictor order used to reconstruct spline bases matches the order used in the fitted landscape model. [Open the R script](Partial_responses_Bayesian_Spatial_model.R).

### 6. Relate predicted connectivity to patch occupancy

The occupancy script reads the independent patch survey, extracts mean predicted connectivity around each patch from a selected raster, and fits a Bernoulli logit model. A Matérn 3/2 spatial effect based on patch coordinates accounts for spatial correlation not explained by connectivity. The script expects the selected raster to exist as `best_Conn` in the current R session; the ZIP does not include a saved selected surface. Its posterior coefficient and fitted occupancy curve refer to the particular surface supplied. [Open the R script](Bayesian_Spatial_occupancy_glm.R).

## Reproducibility and linked files

Run the scripts in the sequence above, supplying the external datasets and the required R packages. The distance analysis can be run after `connectivity.txt` has been created. After the iterative landscape model, explicitly select the desired iteration and provide `best_Conn` to the occupancy script; separately prepare the RDS required for partial responses. In the supplied scripts, the landscape fit and iteration objects remain in the R session, so they are not automatically available after restarting R.

The PDF in this ZIP describes the same workflow. Its six script links open the corresponding `.R.txt` copies as plain text according to the reader's file associations; the Markdown links above point to the original `.R` files. Extract the entire ZIP and keep `README.pdf`, the `.R` files, and their `.R.txt` copies in the same folder. The PDF cannot force a particular text viewer, and some PDF viewers restrict links to local files; the relative links in this Markdown README and the printed filenames provide alternative access. All six script links above refer to files in this ZIP, but the missing input data and intermediate RDS mean that this archive alone cannot reproduce every analysis.
