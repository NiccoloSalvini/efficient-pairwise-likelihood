# Load libraries
library(spdep)
library(dplyr)

# Load data
data("nyc_bor", package = "ptools")
# nyc_boroughs is a dataset included in the spdep package, consisting of polygons of New York City boroughs.

# Construct neighborhood matrix for irregular grid dataset
polyCoords <- st_coordinates(nyc_bor) # non funziona (non ci sono coordinate)
W_irreg <- knn2nb(knearneigh(polyCoords, k = 4))

# Estimate spatial regression model using pairwise likelihood
se_pl <- errorsar(plm(Crime ~ Income, data=nyc_boroughs), listW=W_irreg, SDEM=TRUE)

# Display summary results
summary(se_pl)
