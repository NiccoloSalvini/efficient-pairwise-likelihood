library(RANN)
library(dplyr)
library(tidyr)

source(here("notebooks", "11.2_simulation_functions.R"))


viz_couples <- function(data = data, coppiette = couplets_df) {
  ggplot(data,  aes(x = coord_1, y = coord_2)) +
    geom_point() +
    geom_segment(data = coppiette,
                 aes(x = el_1_coord_1, y = el_1_coord_2, xend = el_2_coord_1, yend = el_2_coord_2),
                 color = "red")+
    geom_text(aes(label = rownames(data)), hjust=0, vjust=1.3)

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

  return(list(beta_hat = beta_hat, sigma_hat_sq = sigma_hat_sq, psi_hat = psi_hat))
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

# Set simulation parameters
phi_values <- c(1, 0.8)
psi_est_values =  c(0.367, 0.286)
sample_sizes <- c(100, 400, 900, 2500)
sample_sizes_big_data <- c(1000, 4000, 9000, 25000) # per 10 rispetto alla dimensiobalità del paper di Arbia
num_replications <- 3
num_replications <- 10000
beta <- 1
sigma_sq <- 1

# Initialize result matrices
bias_matrix <- matrix(0, nrow = length(phi_values), ncol = length(sample_sizes))
relative_bias_matrix <- matrix(0, nrow = length(phi_values), ncol = length(sample_sizes))
mse_matrix <- matrix(0, nrow = length(phi_values), ncol = length(sample_sizes))


# Initialize list to store BML estimates
bml_estimates_list <- list()

# Initialize result tibble
simulation_results <- tibble()

# Loop through different phi values and sample sizes
set.seed(1234)

for (i in 1:length(phi_values)) {
  phi <- phi_values[i]

  for (j in 1:length(sample_sizes)) {
    n <- sample_sizes[j]

    set.seed(1234)
    coord_1 <- runif(n, min = 400, max = 5000)
    coord_2 <- runif(n, min = 400, max = 5000)

    data <- data.frame(
      coord_1 = coord_1,
      coord_2 = coord_2
    )


    couplets  <- GenerateSpatialCouplets(data)

    # Rescale the distances
    scaled_d <- rescale_min(couplets[[3]])

    # custom declined Inverse Correlation
    # Phi <- InverseExponentialCorrelationCustomDecline(scaled_d, phi, esp = 3)
    Phi <- InverseExponentialCorrelation(scaled_d, phi)
    Phi[Phi==1] <- 0
    diag(Phi) = 1

    # Decompose Phi using Cholesky decomposition
    # quando iterazione n = 400 https://stackoverflow.com/questions/60000554/make-a-random-correlation-matrix-semi-definite-positive
    L <- chol(Phi)

    # Generate epsilon values
    epsilon <- GenerateSpatialObservations(n, L)

    # Generate x values
    # questi vanno generati una sola volta con la loro struttura spaziale , quindi con la loro L
    # ma va fatto prima perchp cambiano per la n ma la phi rimane sempre la stessa, cioè 1
    x <- GenerateSpatialObservations(n, L, sd = 2)

    # Generate y values
    y <- beta*x + epsilon

    data$x_val  <- x
    data$epsilon_val  <- epsilon
    data$y_val  <- y

    # regerate couplets (this may be optimised)
    couplets = GenerateSpatialCouplets(data)

    # TODO Couple spatial data, introduce params for buffer


    couplets_values <- couplets[[1]]
    couplets_indexes <- couplets[[2]]

    # Run BML estimation
    for (l in 1:length(psi_est_values)) {
      psi_est <- psi_est_values[l]

      cat("psi hat value", psi_est, "\n")

      bml_estimates <- CalculateBMLEstimates(
        x_i = couplets_values$el_1_x,
        x_l = couplets_values$el_2_x,
        y_i = couplets_values$el_1_y,
        y_l =  couplets_values$el_2_y,
        psi_est
        )

      # Store BML estimates in the list
      bml_estimates_list <- append(bml_estimates_list, bml_estimates)

      # Calculate bias, relative bias, and mean squared error
      bias <- abs(1 - bml_estimates$beta_hat)
      relative_bias <- bias / 1
      mse <- (bml_estimates$beta_hat - mean(bml_estimates$beta_hat))^2 + bias^2

      # Store results in matrices
      bias_matrix[i, j] <- bias
      relative_bias_matrix[i, j] <- relative_bias
      mse_matrix[i, j] <- mse
      print(mse_matrix)


      # Store BML estimates for sigma hat, psi hat, and beta hat
      sigma_hat <- bml_estimates$sigma_hat_sq
      psi_hat_val <- bml_estimates$psi_hat
      beta_hat_val <- bml_estimates$beta_hat

      # Print BML estimates
      cat(paste0("Simulation Parameter: phi =", phi, ", n =", n, "\n"))
      cat(paste0("sigma_hat =", sigma_hat, ", psi_hat =", psi_hat_val, ", beta_hat =", beta_hat_val, "\n\n"))


      result <- tibble(
        phi = phi,
        sample_size = n,
        bias = bias,
        relative_bias = relative_bias,
        mse = mse,
        sigma_hat = sigma_hat,
        psi_hat = psi_hat_val,
        beta_hat = beta_hat_val
      )

      simulation_results <- bind_rows(simulation_results, result)

    }
  }
}


# Print results
print(bias_matrix)
print(relative_bias_matrix)
print(mse_matrix)

# Print simulation results
print(simulation_results)


## measuring algorithm time to couple observations ----
library(purrr)
library(tictoc)
library(ggplot2)
library(microbenchmark)


generateSpatialCoupltesWithTime <- function(n) {
  coord_1 <- runif(n , min = 1000, max = 5000)
  coord_2 <- runif(n, min = 1000, max = 5000)

  data <- data.frame(
    coord_1 = coord_1,
    coord_2 = coord_2
  )

  GenerateSpatialCouplets(data)

}

# test
generateSpatialCoupltesWithTime(100)


## with microbenchmark
## potresti esplorare con profvis i colli di bottiglia
coupling_times = microbenchmark(
  n_100 = generateSpatialCoupltesWithTime(n = 100),
  n_500 = generateSpatialCoupltesWithTime(n = 500),
  n_1000 = generateSpatialCoupltesWithTime(n = 1000),
  n_5000 = generateSpatialCoupltesWithTime(n = 5000),
  n_10000 = generateSpatialCoupltesWithTime(n = 10000),
  n_50000 = generateSpatialCoupltesWithTime(n = 50000),
  n_100000 = generateSpatialCoupltesWithTime(n = 100000),
  times = 10,
  unit = "s")

theme_set(theme_light())

microbenchmark:::autoplot.microbenchmark(coupling_times)
microbenchmark:::boxplot.microbenchmark(coupling_times)

library(dplyr)
as_tibble(coupling_times) %>%
  mutate(expr = as.factor(expr)) %>%
  group_by(expr) %>%
  summarise(mean_time = mean(time)/1000000000) %>%
  ggplot(aes(x = mean_time, y = expr)) +
  geom_line( group = 1) +
  geom_point() +
  labs(
  title = "Coupling Algorithm time (mean time, experiment repetead 10 times)",
    x = "Time (milliseconds)",
    y = ""
    )





