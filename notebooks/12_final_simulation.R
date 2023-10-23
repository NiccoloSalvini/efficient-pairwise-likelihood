library(RANN)
library(dplyr)
library(tidyr)
library(here)

# source(here("notebooks", "13_compare_with_full.R"))
source(here("notebooks", "11.2_simulation_functions.R"))


# Loop through different phi values and sample sizes
# Set simulation parameters
phi_values <- c(1, 0.8)
psi_est_values =  c(0.367, 0.286)
sample_sizes <- c(100, 400, 900, 2500)
not_so_big_sample_sizes <- c(200, 800, 1800, 5000)
sample_sizes_big_data <- c(1000, 4000, 9000, 25000) # per 10 rispetto alla dimensiobalità del paper di Arbia
num_replications <- 3
beta <- 1
sigma_sq <- 1
bml_estimates_list <- list()
simulation_results = tibble()

# Start simulation
for (i in 1:length(phi_values)) {

  phi <- phi_values[i]
  cat("phi hat value:", phi, "\n")

  psi_est <- psi_est_values[i]
  cat("psi hat value:", psi_est, "\n")

  for (j in 1:length(not_so_big_sample_sizes)) {
    n <- not_so_big_sample_sizes[j]

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

    # Inverse Correlation when phi = 1
    Phi_x <- InverseExponentialCorrelation(scaled_d, phi = 1)
    Phi_x[Phi_x==1] <- 0
    diag(Phi_x) <- 1

    # Decompose Phi using Cholesky decomposition
    # quando iterazione n = 400 https://stackoverflow.com/questions/60000554/make-a-random-correlation-matrix-semi-definite-positive
    L_x <- chol(Phi_x)

    # Generate x valuesn(this needs to be non stochastic)
    x <- GenerateSpatialObservations(n, L_x, sd = 2)

    for (k in 1:num_replications) {

      # regenerate spatial structure
      Phi <- InverseExponentialCorrelation(scaled_d, phi)
      Phi[Phi==1] <- 0
      diag(Phi) <- 1

      # Decompose Phi using Cholesky decomposition
      # quando iterazione n = 400 https://stackoverflow.com/questions/60000554/make-a-random-correlation-matrix-semi-definite-positive
      L <- chol(Phi)

      # Generate epsilon values
      epsilon <- GenerateSpatialObservations(n, L)

      # Generate y values
      y <- beta*x + epsilon

      data$x_val  <- x
      data$epsilon_val  <- epsilon
      data$y_val  <- y

      # # regerate couplets (this may be optimised)
      couplets = GenerateSpatialCouplets(data)

      # TODO Couple spatial data, introduce params for buffer

      couplets_values <- couplets[[1]]
      couplets_indexes <- couplets[[2]]

      # Calculate BML estimates
      bml_estimates <- CalculateBMLEstimates(
          x_i = couplets_values$el_1_x,
          x_l = couplets_values$el_2_x,
          y_i = couplets_values$el_1_y,
          y_l =  couplets_values$el_2_y,
          psi_est
      )

      # Print estimates
      cat("--- iteration num:", k, "\n")
      cat(" < estimated params > ", "\n")
      cat("beta hat value:", bml_estimates$beta_hat, "\n")
      cat("sigma hat sq value:", bml_estimates$sigma_hat_sq, "\n")
      cat("psi hat value:", bml_estimates$psi_hat, "\n")

      # Store BML estimates
      bml_estimates_list[[k]] <- bml_estimates

    }

    bml_estimates_res = tibble(bml_estimates_list) %>%
      unnest_wider(bml_estimates_list, names_sep = "_")

    # Calculate bias to retrieve relative bias
    bias_beta_hat <- abs(1 - mean(bml_estimates_res$bml_estimates_list_beta_hat))
    bias_sigma_hat_sq <- abs(1 - mean(bml_estimates_res$bml_estimates_list_sigma_hat_sq))
    bias_psi_hat <- abs(1 - mean(bml_estimates_res$bml_estimates_list_psi_hat))

    # Store relative bias
    relative_bias_beta_hat <- bias_beta_hat / beta
    relative_bias_sigma_hat_sq <- bias_sigma_hat_sq / sigma_sq
    relative_bias_bias_psi_hat <- bias_psi_hat / phi

    # Store Mean Squared Error
    mse_beta_hat <- (mean(bml_estimates_res$bml_estimates_list_beta_hat-mean(c(bml_estimates_res$bml_estimates_list_beta_hat))))^2 + bias_beta_hat^2
    mse_sigma_hat_sq <- (mean(bml_estimates_res$bml_estimates_list_sigma_hat_sq-mean(c(bml_estimates_res$bml_estimates_list_sigma_hat_sq))))^2 + bias_sigma_hat_sq^2
    mse_psi_hat <- (mean(bml_estimates_res$bml_estimates_list_psi_hat-mean(c(bml_estimates_res$bml_estimates_list_psi_hat))))^2 + bias_psi_hat^2


    cat(" < metrics >", "\n")
    cat("relative bias value for beta hat:", mean(relative_bias_beta_hat), "\n")
    cat("relative bias value for sigma hat sq:", mean(relative_bias_sigma_hat_sq), "\n")
    cat("relative bias value for psi hat", mean(relative_bias_bias_psi_hat), "\n")
    cat("mse value for beta hat:", mean(mse_beta_hat), "\n")
    cat("mse value for sigma hat sq:", mean(mse_sigma_hat_sq), "\n")
    cat("mse value for psi hat", mean(mse_psi_hat), "\n")

    result <- tibble(
      phi = phi,
      sample_size = n,
      true_sigma = sigma_sq,
      sigma_hat = mean(bml_estimates_res$bml_estimates_list_sigma_hat_sq),
      true_psi = psi_est,
      psi_hat = mean(bml_estimates_res$bml_estimates_list_psi_hat),
      true_beta = beta,
      beta_hat = mean(bml_estimates_res$bml_estimates_list_beta_hat),
      relative_bias_beta_hat = relative_bias_beta_hat,
      relative_bias_sigma_hat_sq = relative_bias_sigma_hat_sq,
      relative_bias_bias_psi_hat = relative_bias_bias_psi_hat,
      mse_beta_hat = mse_beta_hat,
      mse_sigma_hat_sq = mse_sigma_hat_sq,
      mse_psi_hat = mse_psi_hat
    )

    simulation_results <- bind_rows(simulation_results, result)

  }

}




