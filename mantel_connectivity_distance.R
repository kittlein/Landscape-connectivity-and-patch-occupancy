# Mantel test between genetic connectivity and linear geographic distance

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(sf)
library(vegan)

# Input files

coordinates_file  <- "coordinates.txt"
connectivity_file <- "connectivity.txt"

# Coordinate reference system

crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# Read locality coordinates

coords <- read.table(coordinates_file, header = TRUE, stringsAsFactors = FALSE)

stopifnot(all(c("Locality", "Long", "Lat") %in% names(coords)))

localities_sf <- st_as_sf(coords, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

localities_sf <- st_transform(localities_sf, crs = crs_sinu)

# Read directional genetic connectivity matrix

connectivity <- as.matrix(read.table(connectivity_file, header = FALSE))

stopifnot(nrow(connectivity) == nrow(coords), ncol(connectivity) == nrow(coords))

# Calculate linear geographic distances

xy <- st_coordinates(localities_sf)

linear_distance_km <- as.matrix(dist(xy)) / 1000

diag(linear_distance_km) <- NA

diag(connectivity) <- NA

# Mantel tests for both triangles of the directional connectivity matrix

set.seed(123)

mantel_lower <- mantel(as.dist(linear_distance_km), as.dist(connectivity), permutations = 19999, method = "pearson")

mantel_upper <- mantel(as.dist(linear_distance_km), as.dist(t(connectivity)), permutations = 19999, method = "pearson")

# Join permutation distributions from both Mantel tests

perms <- c(mantel_lower$perm, mantel_upper$perm)

# Observed correlation using both triangles

valid <- !is.na(linear_distance_km) & !is.na(connectivity)

r_obs <- cor(linear_distance_km[valid], connectivity[valid], method = "pearson")

# Calculate empirical probability

p_mantel <- (sum(perms <= r_obs) + 1) / (length(perms) + 1)

print(r_obs)
print(p_mantel)
