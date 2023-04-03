# Load libraries
library(spdep)
library(dplyr)
library(sf)

# Load data
data("nyc_bor", package = "ptools")

# Monte Carlo neighborhood construction
set.seed(1234) # Set seed for reproducibility
n_subset <- 50 # Number of spatial locations in the subset
k <- 5 # Number of nearest neighbors for neighborhood construction
subset_coords <- nyc_bor[sample(nrow(nyc_bor), n_subset),]$geometry %>% st_coordinates()
nn_subset <- knn2nb(knearneigh(subset_coords, k = k))
ind <- knn2list(knearneigh(subset_coords, k = k))
W_irreg <- nb2listw(nn_subset, style = "B")

# Display summary of neighborhood matrix
summary(W_irreg)


## knn sul -1 del cluster
