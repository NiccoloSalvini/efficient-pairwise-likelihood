library(spData)
library(spdep)
library(sp)
library(spatialreg)
library(RANN)
library(dplyr)
library(tidyr)
library(broom)
library(tictoc)

## compare couples algo with full likelihood
source(here("notebooks", "11.2_simulation_functions.R"))

# this needs already a scaled matrix
generateSpatialData <- function(scaled_matrix = scaled_d, phi = 0.8) {

  # Inverse Correlation when phi = 1
  Phi_x <- InverseExponentialCorrelation(scaled_d, phi = 1)
  Phi_x[Phi_x==1] <- 0
  diag(Phi_x) <- 1

  # Decompose Phi using Cholesky decomposition
  # quando iterazione n = 400 https://stackoverflow.com/questions/60000554/make-a-random-correlation-matrix-semi-definite-positive
  L_x <- chol(Phi_x)

  # Generate x valuesn(this needs to be non stochastic)
  x <- GenerateSpatialObservations(n, L_x, sd = 2)

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

  return(data)

}


produceSemEstimates <- function(data) {
  coords <- data[ , c("coord_1", "coord_2")]   # coordinates
  data   <- data[ , c("x_val", "y_val")]          # data
  #crs    <- CRS("+proj=utm +zone=18 +datum=NAD83")# proj4string of coords

  # make the SpatialPointsDataFrame object
  spdf <- SpatialPointsDataFrame(coords      = coords,
                                 data        = data)
  #proj4string = crs)


  # Create spatial weights matrix using inverse distance
  Wd <- dnearneigh(as.matrix(coordinates(spdf)), d1=0, d2=3000)

  W.list <- nb2listw(Wd, style="B", zero.policy=TRUE) # style='B' for binary weights

  # Define the response variable and the covariates:
  y <- data$y_val
  x <- data$x_val

  # Fit the spatial error model using maximum likelihood estimation
  sem <- errorsarlm(y ~ 0 + x, data=data, listw=W.list, Durbin = FALSE)

  return(list(
    beta_hat_fl = sem$coefficients,
    sigma_hat_sq_fl = sem$s2,
    psi_hat_fl = sem$lambda
  )
  )


}

# Loop through different phi values and sample sizes
# Set simulation parameters
phi_values <- c(1, 0.8)
psi_est_values =  c(0.367, 0.286)
sample_sizes <- c(100, 400, 900, 2500)
not_so_big_sample_sizes <- c(200, 800, 1800, 5000)
sample_sizes_big_data <- c(1000, 4000, 9000, 25000) # per 10 rispetto alla dimensiobalità del paper di Arbia
num_replications <- 3
buff_vec = seq(50, 800, by = 150)
beta <- 1
sigma_sq <- 1
bml_estimates_list <- list()
fl_estimates_list <- list()
simulation_results = tibble()

# Start simulation
for (b in 1:length(buff_vec)) {

  for (i in 1:length(phi_values)) {

    phi <- phi_values[i]
    cat("phi hat value:", phi, "\n")

    psi_est <- psi_est_values[i]
    cat("psi hat value:", psi_est, "\n")

    for (j in 1:length(not_so_big_sample_sizes)) {
      n <- not_so_big_sample_sizes[j]

      set.seed(12345)

      coord_1 <- runif(n, min = 400, max = 5000)
      coord_2 <- runif(n, min = 400, max = 5000)

      data <- data.frame(
        coord_1 = coord_1,
        coord_2 = coord_2
      )

      buffer = buff_vec[b]
      couplets  <- GenerateSpatialCoupletsBuffered(data, radius_filter = "custom", buffer = buffer)

      # Rescale the distances
      scaled_d <- rescale_min(couplets[[3]])

      # Generate spatial data for full likelihood estimates
      fl_data = generateSpatialData(scaled_d)

      # Inverse Correlation when phi = 1
      Phi_x <- InverseExponentialCorrelation(scaled_d, phi = 1)
      Phi_x[Phi_x==1] <- 0
      diag(Phi_x) <- 1

      # Decompose Phi using Cholesky decomposition
      # quando iterazione n = 400 https://stackoverflow.com/questions/60000554/make-a-random-correlation-matrix-semi-definite-positive
      L_x <- chol(Phi_x)

      # Generate x valuesn(this needs to be non stochastic)
      x <- GenerateSpatialObservations(n, L_x, sd = 2)

      start = tic(quiet = F)

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
        couplets = GenerateSpatialCoupletsBuffered(data, radius_filter = "custom", buffer = buffer)

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


        fl_estimates = produceSemEstimates(
          data = fl_data
        )

        # Print estimates for pairwise likelihood
        cat("--- iteration num:", k, "\n")
        cat(" < estimated params > ", "\n")

        cat("beta hat value for pairwise:", bml_estimates$beta_hat_pl, "\n")
        cat("sigma hat sq value for pairwise:", bml_estimates$sigma_hat_sq_pl, "\n")
        cat("psi hat value for pairwise:", bml_estimates$psi_hat_pl, "\n")

        cat("beta hat value for full likelihood:", fl_estimates$beta_hat_fl, "\n")
        cat("sigma hat sq value for full likelihood:", fl_estimates$sigma_hat_sq_fl, "\n")
        cat("psi hat value for full likelihood:", fl_estimates$psi_hat_fl, "\n")

        # Store BML estimates
        bml_estimates_list[[k]] <- bml_estimates
        fl_estimates_list[[k]] <- fl_estimates

      }

      end = toc(quiet = F)

      estimates_res = bind_cols(bind_rows(bml_estimates_list), bind_rows(fl_estimates_list))

      # Calculate bias to retrieve relative bias for pl and fl
      bias_beta_hat_pl <- abs(1 - mean(estimates_res$beta_hat_pl))
      bias_sigma_hat_sq_pl <- abs(1 - mean(estimates_res$sigma_hat_sq_pl))
      bias_psi_hat_pl <- abs(1 - mean(estimates_res$psi_hat_pl))

      bias_beta_hat_fl <- abs(1 - mean(estimates_res$beta_hat_fl))
      bias_sigma_hat_sq_fl <- abs(1 - mean(estimates_res$sigma_hat_sq_fl))
      bias_psi_hat_fl <- abs(1 - mean(estimates_res$psi_hat_fl))

      # Store relative bias for pl and fl
      relative_bias_beta_hat_pl <- bias_beta_hat_pl / beta
      relative_bias_sigma_hat_sq_pl <- bias_sigma_hat_sq_pl / sigma_sq
      relative_bias_psi_hat_pl <- bias_psi_hat_pl / phi

      relative_bias_beta_hat_fl  <- bias_beta_hat_fl / beta
      relative_bias_sigma_hat_sq_fl <- bias_sigma_hat_sq_fl / sigma_sq
      relative_bias_psi_hat_fl <- bias_psi_hat_fl / phi

      # Store Mean Squared Error for pl and fl (as defined in Arbia 2014)
      mse_beta_hat_pl <- (mean(estimates_res$beta_hat_pl-mean(c(estimates_res$beta_hat_pl))))^2 + bias_beta_hat_pl^2
      mse_sigma_hat_sq_pl <- (mean(estimates_res$sigma_hat_sq_pl-mean(c(estimates_res$sigma_hat_sq_pl))))^2 + bias_sigma_hat_sq_pl^2
      mse_psi_hat_pl <- (mean(estimates_res$psi_hat_pl-mean(c(estimates_res$psi_hat_pl))))^2 + bias_psi_hat_pl^2

      mse_beta_hat_fl <- (mean(estimates_res$beta_hat_fl-mean(c(estimates_res$beta_hat_fl))))^2 + bias_beta_hat_fl^2
      mse_sigma_hat_sq_fl <- (mean(estimates_res$sigma_hat_sq_fl-mean(c(estimates_res$sigma_hat_sq_fl))))^2 + bias_sigma_hat_sq_fl^2
      mse_psi_hat_fl <- (mean(estimates_res$psi_hat_fl-mean(c(estimates_res$psi_hat_fl))))^2 + bias_psi_hat_fl^2


      cat(" < metrics for pl and fl>", "\n")
      cat("relative bias value for beta hat PL:", mean(relative_bias_beta_hat_pl), "\n")
      cat("relative bias value for sigma hat sq PL:", mean(relative_bias_sigma_hat_sq_pl), "\n")
      cat("relative bias value for psi hat PL", mean(relative_bias_psi_hat_pl), "\n")
      cat("mse value for beta hat PL:", mean(mse_beta_hat_pl), "\n")
      cat("mse value for sigma hat sq PL:", mean(mse_sigma_hat_sq_pl), "\n")
      cat("mse value for psi hat PL:", mean(mse_psi_hat_pl), "\n")

      cat("relative bias value for beta hat FL:", mean(relative_bias_beta_hat_fl), "\n")
      cat("relative bias value for sigma hat sq FL:", mean(relative_bias_sigma_hat_sq_fl), "\n")
      cat("relative bias value for psi hat FL:", mean(relative_bias_psi_hat_fl), "\n")
      cat("mse value for beta hat FL:", mean(mse_beta_hat_fl), "\n")
      cat("mse value for sigma hat sq FL:", mean(mse_sigma_hat_sq_fl), "\n")
      cat("mse value for psi hat FL:", mean(mse_psi_hat_fl), "\n")


      result <- tibble(
        phi = phi,
        buffer = buffer,
        sample_size = n,
        true_sigma = sigma_sq,
        sigma_hat_pairwise = mean(estimates_res$sigma_hat_sq_pl),
        sigma_hat_fullLik = mean(estimates_res$sigma_hat_sq_fl),
        true_psi = psi_est,
        psi_hat_pairwise = mean(estimates_res$psi_hat_pl),
        psi_hat_fullLik = mean(estimates_res$psi_hat_fl),
        true_beta = beta,
        beta_hat_pairwise = mean(estimates_res$beta_hat_pl),
        beta_hat_fullLik = mean(estimates_res$beta_hat_fl),
        relative_bias_beta_hat_pairwise = relative_bias_beta_hat_pl,
        relative_bias_beta_hat_fullLik = relative_bias_beta_hat_fl,
        relative_bias_sigma_hat_sq_pairwise = relative_bias_sigma_hat_sq_fl,
        relative_bias_sigma_hat_sq_fullLik = relative_bias_sigma_hat_sq_fl,
        relative_bias_psi_hat_pairwise = relative_bias_psi_hat_pl,
        relative_bias_psi_hat_fullLik = relative_bias_psi_hat_fl,
        mse_beta_hat_pairwise = mse_beta_hat_pl,
        mse_beta_hat_fullLik = mse_beta_hat_fl,
        mse_sigma_hat_sq_pairwise = mse_sigma_hat_sq_pl,
        mse_sigma_hat_sq_fuillLik = mse_sigma_hat_sq_fl,
        mse_psi_hat_pairwise = mse_psi_hat_pl,
        mse_psi_hat_fullLik = mse_psi_hat_fl,
        time_elapsed = end$toc - end$tic
      )

      simulation_results <- bind_rows(simulation_results, result)

    }

  }
}


saveRDS(simulation_results, file = here(paste0("data/simres_buff", Sys.Date(),".rds")))
