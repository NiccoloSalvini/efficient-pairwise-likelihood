# First, create a sample dataset of geo-referenced points
set.seed(1)
n_points <- 100

lat <- runif(n_points, 40, 45)
lon <- runif(n_points, -80, -75)

points_df <- data.frame(lat = lat, lon = lon)

# Calculate the pairwise distances between all the points
library(geosphere)

dist_matrix <- distm(points_df[, c("lon", "lat")], fun = distHaversine)

# Replace the diagonal with an arbitrarily large number
diag(dist_matrix) <- Inf

# Find the indices of the two closest points
# Keep track of which points have already been matched
matches <- data.frame(from = numeric(0), to = numeric(0))

while (nrow(matches) < n_points / 2) {

  # Find the indices of the minimum value
  min_index <- which(dist_matrix == min(dist_matrix), arr.ind = TRUE)

  # Get the indices of the two closest points
  closest_indices <- min_index[, 2:1]

  # Check if either point has already been matched
  if ((closest_indices[1] %in% matches$from) || (closest_indices[2] %in% matches$to)) {
    # If one has, remove the minimum distance from the matrix and continue
    dist_matrix[closest_indices[1], closest_indices[2]] <- Inf
  } else {
    # Otherwise, add the pair to the matches dataframe
    matches[nrow(matches) + 1,] <- closest_indices

    # Remove the matched points from the distance matrix
    dist_matrix[closest_indices[1], ] <- Inf
    dist_matrix[, closest_indices[1]] <- Inf
    dist_matrix[closest_indices[2], ] <- Inf
    dist_matrix[, closest_indices[2]] <- Inf
  }
}



library(ggplot2)
# Join the matches dataframe with the original points dataframe
matches_df <- cbind(from = points_df[matches$from, ], to = points_df[matches$to, ])

# Plot the points and the matches
ggplot() +
  geom_point(data = points_df, aes(x = lon, y = lat)) +
  geom_segment(data = matches_df, aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), color = "red", size = 1)

##############################################






# First, create a sample dataset of geo-referenced points
set.seed(1)
n_points <- 100

lat <- runif(n_points, 40, 45)
lon <- runif(n_points, -80, -75)

points_df <- data.frame(lat = lat, lon = lon)


# Calculate mean distance and filter matches dataframe
# recalculate distance matrix since you have overwritten the previous
library(geosphere)
dist_matrix <- distm(points_df[, c("lon", "lat")], fun = distHaversine)
mean_dist <- mean(dist_matrix[!is.infinite(dist_matrix)])
matches <- matches[dist_matrix[as.matrix(matches)] < 0.2 * mean_dist,]

# Randomly select subset of shortest distanced points within threshold distance
threshold <- 233845.3
subset_size <- 2

# Keep track of shortest distances and pair indices
shortest_dists <- numeric()
shortest_pairs <- list()

# Loop through distances in matches dataframe
for (i in 1:nrow(matches)) {
  dist <- distm(points_df[matches[i, "from"], c("lon", "lat")], points_df[matches[i, "to"], c("lon", "lat")], fun = distHaversine)

  # If within threshold and shorter than current shortest distance,
  # update shortest_dists and shortest_pairs
  if (is.finite(dist) && dist <= threshold) {
    if (length(shortest_dists) < subset_size) {
      shortest_dists <- c(shortest_dists, dist)
      shortest_pairs[[length(shortest_pairs) + 1]] <- matches[i,]
    } else {
      if (dist < max(shortest_dists)) {
        max_index <- which.max(shortest_dists)
        shortest_dists[max_index] <- dist
        shortest_pairs[[max_index]] <- matches[i,]
      }
    }
  }
}

library(ggplot2)


## manca random_pairs_df, dove è?


# Plot all geo-referenced points
# First, create a new dataframe that includes the coordinates of the from and to points
random_pairs_coords <- data.frame(lon_from = points_df$lon[random_pairs_df$from], lat_from = points_df$lat[random_pairs_df$from], lon_to = points_df$lon[random_pairs_df$to], lat_to = points_df$lat[random_pairs_df$to])

# Plot all geo-referenced points
ggplot(points_df, aes(x = lon, y = lat)) +
  geom_point(alpha = 0.3, size = 0.5) +
  theme_bw() +
  labs(title = "Spatial Dataset of Points")

# Plot selected random couplets within threshold distance
ggplot() +
  geom_point(data = points_df, aes(x = lon, y = lat), alpha = 0.3, size = 0.5) +
  geom_point(data = random_pairs_coords, aes(x = lon_from, y = lat_from), alpha = 0.7, size = 2, color = "red") +
  geom_point(data = random_pairs_coords, aes(x = lon_to, y = lat_to), alpha = 0.7, size = 2, color = "red") +
  geom_segment(data = random_pairs_coords, aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to), color = "red", alpha = 0.3) +
  theme_bw() +
  labs(title = "Random Subset of Couplets within Threshold Distance")



## retry
# Original code to create dataset and find pairwise distances
# Original code to create dataset and find pairwise distances
set.seed(1)
n_points <- 100
lat <- runif(n_points, 40, 45)
lon <- runif(n_points, -80, -75)
points_df <- data.frame(lat = lat, lon = lon)
library(geosphere)
dist_matrix <- distm(points_df[, c("lon", "lat")], fun = distHaversine)
diag(dist_matrix) <- Inf

# Algorithm to find nearest couplets and filter by mean distance
matches <- data.frame(from = numeric(0), to = numeric(0))
while (nrow(matches) < n_points / 2) {
  min_index <- which(dist_matrix == min(dist_matrix), arr.ind = TRUE)
  closest_indices <- min_index[, 2:1]
  if ((closest_indices[1] %in% matches$from) || (closest_indices[2] %in% matches$to)) {
    dist_matrix[closest_indices[1], closest_indices[2]] <- Inf
  } else {
    matches[nrow(matches) + 1,] <- closest_indices
    dist_matrix[closest_indices[1], ] <- Inf
    dist_matrix[, closest_indices[1]] <- Inf
    dist_matrix[closest_indices[2], ] <- Inf
    dist_matrix[, closest_indices[2]] <- Inf
  }
}
dist_matrix <- distm(points_df[, c("lon", "lat")], fun = distHaversine)
mean_dist <- mean(dist_matrix[!is.infinite(dist_matrix)])
matches <- matches[dist_matrix[as.matrix(matches)] < 0.2 * mean_dist,]

# Randomly selecting couplets and excluding neighboring ones within a given distance threshold
n <- 10 # number of randomly selected couplets
d <- 10000 # threshold distance in meters
random_couplets <- matches[sample(nrow(matches), n),]

midpoints <- apply(random_couplets, 1, function(x) {
  with(points_df, (data.frame(lat = (lat[x[1]] + lat[x[2]]) / 2,
                              lon = (lon[x[1]] + lon[x[2]]) / 2)))
})

circles <- lapply(midpoints, function(midpoint) {
  circle <- tryCatch(gcIntermediate(c(midpoint[[1]]$lon, midpoint[[1]]$lat), c(midpoint[[1]]$lon,
                                    midpoint[[1]]$lat), d, addStartEnd = TRUE),
                     error=function(e) NULL)
  if (!is.null(circle)) {
    circle_df <- data.frame(
      lon = circle[,"lon"],
      lat = circle[,"lat"]
    )
  } else {
    circle_df <- NULL
  }
  circle_df
})

circles <- circles[!sapply(circles, is.null)]

excluded_couplets <- c()
for (i in seq_along(random_couplets[,1])) {
  candidates <- setdiff(matches, random_couplets[i,])
  circle_center <- midpoints[i]
  if (!is.null(circles[[i]])) {
    within_circle <- apply(candidates, 1, function(x) {
      with(points_df, (distm(data.frame(lon = c(lon[x[1]], lon[x[2]]),
                                        lat = c(lat[x[1]], lat[x[2]])),
                             fun = distHaversine) <= d))
    })
    excluded_couplets <- c(excluded_couplets, candidates[rowSums(within_circle) > 0,])
  }
}
final_couplets <- setdiff(matches, excluded_couplets)
