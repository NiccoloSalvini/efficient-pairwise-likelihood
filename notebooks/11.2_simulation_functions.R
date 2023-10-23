library(RANN)
library(dplyr)
library(tidyr)


viz_couples <- function(data = data, coppiette = couplets_df) {
  ggplot(data,  aes(x = coord_1, y = coord_2)) +
    geom_point() +
    geom_segment(data = coppiette,
                 aes(x = el_1_coord_1, y = el_1_coord_2, xend = el_2_coord_1, yend = el_2_coord_2),
                 color = "red",  size = 1.1)+
    geom_text(aes(label = rownames(data), alpha = .3), hjust=0, vjust=1.3, check_overlap = TRUE, size = 3) +
    theme(legend.position='none')

}

GenerateSpatialCouplets<- function(data) {

  # compute nearest neighbours using kdtree algorithm
  nn <- nn2(data, data, treetype = "kd", searchtype = "standard", k = 2)

  # Extract the index and distance information from nn output
  nn_idx <- nn$nn.idx
  nn_dists <- nn$nn.dists

  n = nrow(data)

  # Create an empty distance matrix
  distance_matrix <- matrix(0, nrow = n, ncol = n)

  # initialize an empty list for paired points
  couplets <- list()

  # vector to keep track of the paired points
  paired_points <- c()

  # pair up points such that each point is only paired once
  # Fill the distance matrix with the pairwise distances
  for (i in seq_len(nrow(nn$nn.idx))) {
    point1_index <- nn$nn.idx[i, 1]
    point2_index <- nn$nn.idx[i, 2]
    dist <- nn_dists[i, 2]

    if (point1_index %in% paired_points || point2_index %in% paired_points) {
      next
    } else {
      # pair the two points and add them to the list and the vector of paired points
      point1 <- data[point1_index, ]
      point2 <- data[point2_index, ]
      distance_matrix[point1_index, point2_index] <- dist
      distance_matrix[point2_index, point1_index] <- dist

      couplets[[length(couplets) + 1]] <- list(point1, point2)
      paired_points <- c(paired_points, point1_index, point2_index)
    }
  }

  # convert couplets list to a data frame
  if(couplets[[1]][[1]] %>% colnames() %>% length() == 2){

    couplets_df <- data.frame(do.call(rbind, couplets)) %>%
      unnest_wider(c("X1", "X2"), names_sep = "_", names_repair = "check_unique")

    names(couplets_df) = c("el_1_coord_1",
                           "el_1_coord_2",
                           "el_2_coord_1",
                           "el_2_coord_2")


  } else {
    couplets_df <- data.frame(do.call(rbind, couplets)) %>%
      unnest_wider(c("X1", "X2"), names_sep = "_", names_repair = "check_unique")
    names(couplets_df) = c("el_1_coord_1",
                           "el_1_coord_2",
                           "el_1_x",
                           "el_1_epsilon",
                           "el_1_y",
                           "el_2_coord_1",
                           "el_2_coord_2",
                           "el_2_x",
                           "el_2_epsilon",
                           "el_2_y")
  }

  return(list(couplets_df, paired_points, distance_matrix))
}




GenerateSpatialCoupletsBuffered<- function(data, radius_filter = "min", buffer = 50) {

  d <- as.numeric(as.matrix(dist(data)))
  d[d<=0] <- NA

  switch(radius_filter,
         "min" = {
           initial_radius = min(d, na.rm=TRUE)
         },
         "mean" = {
           initial_radius = mean(d, na.rm=TRUE)
         },
         "max" = {
           initial_radius = max(d, na.rm=TRUE)
         },
         "median" = {
           initial_radius = median(d, na.rm=TRUE)
         },
         "custom" = {
           initial_radius = min(d, na.rm=TRUE) + buffer
         },
         stop(paste0("No method for ", radius_filter))
  )

  cat(initial_radius)

  # compute nearest neighbours using kdtree algorithm
  nn <- nn2(data, data, treetype = "kd", searchtype = "radius", k = 2, radius = initial_radius)

  # Extract the index and distance information from nn output
  nn_idx <- nn$nn.idx
  nn_dists <- nn$nn.dists

  n = nrow(data)

  # Create an empty distance matrix
  distance_matrix <- matrix(0, nrow = n, ncol = n)

  # initialize an empty list for paired points
  couplets <- list()

  # vector to keep track of the paired points
  paired_points <- c()

  # pair up points such that each point is only paired once
  # Fill the distance matrix with the pairwise distances
  for (i in seq_len(nrow(nn$nn.idx))) {
    point1_index <- nn$nn.idx[i, 1]
    point2_index <- nn$nn.idx[i, 2]
    dist <- nn_dists[i, 2]

    if (point1_index %in% paired_points || point2_index %in% paired_points || point2_index == 0) {
      next
    } else {
      # pair the two points and add them to the list and the vector of paired points
      point1 <- data[point1_index, ]
      point2 <- data[point2_index, ]
      distance_matrix[point1_index, point2_index] <- dist
      distance_matrix[point2_index, point1_index] <- dist

      couplets[[length(couplets) + 1]] <- list(point1, point2)
      paired_points <- c(paired_points, point1_index, point2_index)
    }
  }

  # convert couplets list to a data frame
  if(couplets[[1]][[1]] %>% colnames() %>% length() == 2){

    couplets_df <- data.frame(do.call(rbind, couplets)) %>%
      unnest_wider(c("X1", "X2"), names_sep = "_", names_repair = "check_unique")

    names(couplets_df) = c("el_1_coord_1",
                           "el_1_coord_2",
                           "el_2_coord_1",
                           "el_2_coord_2")


  } else {
    couplets_df <- data.frame(do.call(rbind, couplets)) %>%
      unnest_wider(c("X1", "X2"), names_sep = "_", names_repair = "check_unique")
    names(couplets_df) = c("el_1_coord_1",
                           "el_1_coord_2",
                           "el_1_x","Chiara Perricone (ICONSULTING)" <c.perricone@iconsulting.biz>
                           "el_1_epsilon",
                           "el_1_y",
                           "el_2_coord_1",
                           "el_2_coord_2",
                           "el_2_x",
                           "el_2_epsilon",
                           "el_2_y")
  }

  return(list(couplets_df, paired_points, distance_matrix))
}


# Define functions to calculate correlation and generate spatially dependent observations
InverseExponentialCorrelation <- function(d, phi) {
  exp(-d^2 / phi)
}

InverseExponentialCorrelationDubin <- function(d, phi) {
  exp(-d / phi)

}

# custom Inverse
InverseExponentialCorrelationCustomDecline <- function(d, phi, esp) {
  exp(-d^esp / phi)
}


# This needs to be edited (it needs to generate lat and long and relative x and y values for each couplet)
GenerateSpatialObservations <- function(n, L, sd = 1) {
  z <- rnorm(n, sd = sd)
  L %*% z
}

# Calculate BML estimates using closed form expressions
# questa va rivista
CalculateBMLEstimates <- function(x_i, x_l, y_i, y_l, psi_hat) {
  q <- length(x_i)
  alpha1 <- sum(x_i^2) + sum(x_l^2)
  alpha2 <- sum(y_i^2) + sum(y_l^2)
  alpha3 <- sum(x_i*y_i) + sum(x_l*y_l)
  alpha4 <- sum(x_i*y_l) + sum(x_l*y_i)
  alpha5 <- sum(x_i*x_l)
  alpha6 <- sum(y_i*y_l)

  beta_hat <- (alpha3 - psi_hat * alpha4) / (alpha1 - 2 * psi_hat * alpha5)
  sigma_hat_sq <- (alpha2 + beta_hat^2 * alpha1 - 2 * beta_hat * alpha3 -
                     2 * psi_hat * alpha6 - 2 * psi_hat * beta_hat^2 * alpha5 +
                     2 * psi_hat * beta_hat * alpha4) / (2 * q * (1 - psi_hat^2))
  psi_hat <- (alpha6 - beta_hat * alpha4 + beta_hat^2 * alpha5) / (q * sigma_hat_sq)

  return(
    list(
      beta_hat_pl = beta_hat,
      sigma_hat_sq_pl = sigma_hat_sq,
      psi_hat_pl = psi_hat))
}


rescale_min <- function(x, verbose = FALSE) {

  if(verbose){
    cat("Diagnostic for distance matrix\n")
  }
  min_x <- min(x[x != 0])  # Exclude diagonal elements

  if(verbose){
    cat("min for matrix d couple distances", min_x, "\n")
  }

  max_x <- max(x)

  if(verbose){
    cat("max for matrix d couple distances", max_x, "\n")
  }

  # Scale the values to the desired range
  scaled_x <- (x - min_x) / (max_x - min_x) * 4 + 1

  min_scaled_x <- min(scaled_x[scaled_x != 0])

  if(verbose){
    cat("most values now are ", min_scaled_x, "\n")
  }

  # Replace values which are below recaled mean with 0
  scaled_x[scaled_x<=min_scaled_x] <- 0
  scaled_x[scaled_x>min_scaled_x] <- 1

  scaled_x
}
