# Load data
data("nyc_boroughs")
# nyc_boroughs is a dataset included in the spdep package, consisting of polygons of New York City boroughs.

# Construct Voronoi tessellation neighborhood matrix
nyc_boroughs_coords <- st_coordinates(nyc_boroughs)
n_neighbors <- 5 # number of neighbors
dists <- nbdists(knearneigh(nyc_boroughs_coords, k = n_neighbors)) x
W_voronoi <- nb2mat(dists, glist = NULL, style = "B")

# Display summary of neighborhood matrix
summary(W_voronoi)
