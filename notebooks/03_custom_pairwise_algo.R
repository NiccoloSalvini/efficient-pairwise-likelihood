library(tictoc)

# First, create a sample dataset of geo-referenced points
set.seed(1)
n_points <- 100
treshold_dist <- 2


lat <- runif(n_points, 40, 45)
lon <- runif(n_points, -80, -75)

points_df <- data.frame(lat = lat, lon = lon)

# Calculate the pairwise distances between all the points
library(geosphere)

mean_dist <- mean(dist_matrix[lower.tri(dist_matrix)])

# Set distances greater than 2 times the mean to Inf
dist_matrix[dist_matrix > discard_treshold * mean_dist] <- Inf


# Replace the diagonal with an arbitrarily large number
diag(dist_matrix) <- Inf

# Find the indices of the two closest points
# Keep track of which points have already been matched
matches <- data.frame(from = numeric(0), to = numeric(0))

tic()
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
toc()

# Plot the matched pairs with ggplot2
library(ggplot2)

# Join the matches dataframe with the original points dataframe
matches_df <- cbind(from = points_df[matches$from, ], to = points_df[matches$to, ])

# Plot the points and the matches
ggplot() +
  geom_point(data = points_df, aes(x = lon, y = lat)) +
  geom_segment(data = matches_df, aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), color = "red", size = 1)





## function wrapper
library(geosphere)
library(ggplot2)

match_nearest_points <- function(lat, lon, plot_output = FALSE, dist_type = "distHaversine") {
  # Check that the lat and lon vectors are of equal length
  if (!length(lat) == length(lon)) {
    stop("Error: length of lat and lon vectors must be equal.")
  }

  # Check that lat and lon vectors are numeric
  if(!is.numeric(lat) || !is.numeric(lon)) {
    stop("Error: lat and lon vectors must be numeric.")
  }

  # Create a data frame of the input points
  points_df <- data.frame(lat = lat, lon = lon)

  # Calculate the pairwise distances between all points
  dist_matrix <- distm(points_df[, c("lon", "lat")], fun = distHaversine)

  # Replace the diagonal with an arbitrarily large number
  diag(dist_matrix) <- Inf

  # Find the indices of the two closest points
  # Keep track of which points have already been matched
  matches <- data.frame(from = numeric(0), to = numeric(0))

  while (nrow(matches) < length(lat) / 2) {
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

  # Plot the matched pairs with ggplot2 if requested
  if (plot_output) {
    # Join the matches dataframe with the original points dataframe
    matches_df <- cbind(from = points_df[matches$from, ], to = points_df[matches$to, ])

    # Plot the points and the matches
    ggplot() +
      geom_point(data = points_df, aes(x = lon, y = lat)) +
      geom_segment(data = matches_df, aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), color = "red", size = 1)
  }

  # Return the matched pairs as a data frame
  return(data.frame(from = points_df[matches$from, ], to = points_df[matches$to, ]))
}



## calcolo distanzza media, se eccede 2/3 volte viua
## prendere i massimi e toglierli
## prendo punto da coppia (medio) e poi faccio buffer circolare per escludere coppie vicine
## su griglia regolare
