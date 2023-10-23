# Load data
data("nyc_boroughs")
# nyc_boroughs is a dataset included in the spdep package, consisting of polygons of New York City boroughs.

# Construct kernel-based neighborhood matrix
nyc_boroughs_coords <- st_coordinates(nyc_boroughs)
W_kernel <- spkernelnb(coordinates = nyc_boroughs_coords, alpha = 1.5, knn = 5)

# Display summary of neighborhood matrix
summary(W_kernel)


## fare
