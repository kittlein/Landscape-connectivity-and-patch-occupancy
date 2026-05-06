# Landscape-connectivity-and-patch-occupancy
This repository provides the data and code used in the manuscript *Linking genetic connectivity and habitat occupancy through environmental resistance in fragmented landscapes*. It includes scripts for processing genetic and environmental data, building resistance-based connectivity models, generating spatial predictions, and evaluating the association between predicted connectivity and habitat patch occupancy.

### Overview

The study evaluates genetic connectivity in *Ctenomys mendocinus* across the eastern portion of its range and examines whether connectivity inferred from environmental resistance predicts habitat patch occupancy. The workflow integrates microsatellite genetic data, environmental raster layers, and iterative modeling in `R` to generate a spatially explicit connectivity surface and test its relationship with field-based occupancy data.

### Data included in this repository

#### Microsatellite genotypes

Microsatellite genotype data are provided in `micros_sat.csv`. This file contains the multilocus genotypes used to estimate genetic differentiation and directional connectivity among sampled localities.

These data were used for:

- population-level genetic characterization
- estimation of pairwise differentiation among localities
- inference of relative directional connectivity
- construction of the response matrix used in landscape resistance modeling

#### Environmental layers

Environmental predictors are provided in `environment.tif` as a raster brick including all layers described in the manuscript. These layers represent the environmental conditions used to model connectivity across the landscape.

The raster brick includes the variables referenced in the manuscript, including:

- habitat/environment classes
- land-cover change
- mean NDVI
- NDVI variability
- elevation (DEM)
- sand content
- clay content
- surface water

These layers were used to summarize environmental conditions along least-cost routes and to predict connectivity spatially across the study region.

### Analytical workflow

The analyses were conducted in `R` and reproduce the main steps described in the manuscript and Supporting Information. In broad terms, the workflow includes:

1. importing microsatellite genotype data from `micros_sat.csv`
2. estimating genetic differentiation and directional connectivity among localities
3. loading `environment.tif` as a raster brick of environmental predictors
4. generating least-cost paths among localities
5. extracting environmental summaries along routes
6. fitting iterative Random Forest models to predict connectivity
7. projecting predicted connectivity across the landscape
8. evaluating the association between predicted connectivity and habitat patch occupancy

### Code availability

All code used to run the analyses is provided in `R` format in this repository. The scripts implement the full workflow for:

- data import and preprocessing
- microsatellite-based connectivity estimation
- raster preparation and environmental extraction
- least-cost path analysis
- Random Forest modeling
- iterative connectivity-surface updating
- occupancy validation and figure generation

### Reproducibility

Together, `micros_sat.csv`, `environment.tif`, and the accompanying `R` scripts provide the core inputs required to reproduce the main analytical pipeline presented in the manuscript. These materials are intended to support transparency, reproducibility, and reuse in comparable landscape genetics studies.

### Notes

Field-based patch occupancy data and genetic sampling were obtained independently, allowing connectivity inferred from genetic data to be evaluated against an external occupancy dataset. The repository is therefore designed to document both the analytical workflow and the data structure underlying the reported results.
