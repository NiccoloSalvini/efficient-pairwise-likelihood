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


# produceSemEstimates(data = generateSpatialData(n = 100))
# produceSemEstimates(data = generateSpatialData(n = 400))
# produceSemEstimates(data = generateSpatialData(n = 900))
# produceSemEstimates(data = generateSpatialData(n = 2500))


# Loop through different phi values and sample sizes
# Set simulation parameters
phi_values <- c(1, 0.8)
psi_est_values =  c(0.367, 0.286)
sample_sizes <- c(100, 400, 900, 2500)
not_so_big_sample_sizes <- c(200, 800, 1800, 5000)
sample_sizes_big_data <- c(1000, 4000, 9000, 25000) # per 10 rispetto alla dimensiobalità del paper di Arbia
num_replications <- 100
beta <- 1
sigma_sq <- 1
bml_estimates_list <- list()
fl_estimates_list <- list()
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



## coupling times ---- :
## measuring algorithm time to couple observations ----
library(purrr)
library(tictoc)
library(ggplot2)
library(dplyr)
library(microbenchmark)


generateSpatialCoupltesWithTime <- function(n, phi = 0.8) {
  coord_1 <- runif(n , min = 1000, max = 5000)
  coord_2 <- runif(n, min = 1000, max = 5000)

  data <- data.frame(
    coord_1 = coord_1,
    coord_2 = coord_2
  )

  couplets = GenerateSpatialCouplets(data)

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

}

# test
generateSpatialCoupltesWithTime(100)


## with microbenchmark
## potresti esplorare con profvis i colli di bottiglia
coupling_times = microbenchmark(
  n_100_08 = generateSpatialCoupltesWithTime(n = 100),
  n_100_1 = generateSpatialCoupltesWithTime(n = 100, phi = 1),
  n_500_08 = generateSpatialCoupltesWithTime(n = 500),
  n_500_1 = generateSpatialCoupltesWithTime(n = 500),
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


# as_tibble(coupling_times) %>%
#   mutate(expr = as.factor(expr)) %>%
#   group_by(expr) %>%
#   summarise(mean_time = mean(time)/1000000000) %>%
#   ggplot(aes(x = mean_time, y = expr)) +
#   geom_line( group = 1) +
#   geom_point() +
#   labs(
#     title = "Coupling Algorithm time (mean time, experiment repetead 10 times)",
#     x = "Time (milliseconds)",
#     y = ""
#   )
#


Coupling_Times <- coupling_times %>%  group_by(expr) %>%  summarise(median(time/1000))
Coupling_Times['ops'] <- c(100,500,1000,5000,10000, 50000, 100000)
names(Coupling_Times)[2] <- 'median' ;

ggplot(Coupling_Times,aes(ops,median),size = 1.08) +
  geom_point(aes(colour = 'data curve'),size=1.1) +
  geom_line(aes(colour = 'data curve')) +
  stat_function(fun = nlogn,size = 1.08,aes(colour = 'n*logn curve')) +
  labs(x = 'length of array',y = 'median of computing time [miliseconds]')


## TODO
## - rimane da comparare accuratezza parametri di Pairwise contro full-likelihood
## - rimane da comparae computation complexity across different data size
## - reoptimise algorithm, risottomettilo a ChatGPT
## - scrivere articolino (3/4) pagine con solo risultati e misurazione dei tempi
## - vedere fino a che puinto cambiano le cose se il raggio dell'algo cresce.
## - vedere fino a che punto cambiano le cose quando il processo generativo cambia
## - idea altro articolo, testare la metodologia su altro tipo di regressione spaziale, SARAR, whatever whatever




