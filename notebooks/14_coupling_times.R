## coupling times:
## measuring algorithm time to couple observations ----
library(purrr)
library(tictoc)
library(ggplot2)
library(dplyr)
library(microbenchmark)
library(here)

source(here("notebooks", "11.2_simulation_functions.R"))

## grouping couples how long it takes ----
justGenerateCouples <- function(n) {
  coord_1 <- runif(n , min = 1000, max = 5000)
  coord_2 <- runif(n, min = 1000, max = 5000)

  data <- data.frame(
    coord_1 = coord_1,
    coord_2 = coord_2
  )

  couplets = GenerateSpatialCouplets(data)
}

## I assume this takes approx KDtree times
coupling_times = microbenchmark(
  n_100 = justGenerateCouples(n = 100),
  n_200 = justGenerateCouples(n = 200),
  n_500 = justGenerateCouples(n = 500),
  n_700 = justGenerateCouples(n = 700),
  n_1000 = justGenerateCouples(n = 1000),
  n_2000 = justGenerateCouples(n = 2000),
  n_5000 = justGenerateCouples(n = 5000),
  n_10000 = justGenerateCouples(n = 10000),
  n_50000 = justGenerateCouples(n = 50000),
  # n_100000 = justGenerateCouples(n = 100000),
  times = 10,
  unit = "s"
  )

theme_set(theme_light())

microbenchmark:::autoplot.microbenchmark(coupling_times)
microbenchmark:::boxplot.microbenchmark(coupling_times)


## generate couples and spatial strucure (this needed for simulation) ----
## some

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

  couplets = GenerateSpatialCouplets(data)

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

}


# test generate with time!
generateSpatialCoupltesWithTime(100)


## with microbenchmark
## potresti esplorare con profvis i colli di bottiglia
coupling_time_pairwise = microbenchmark(
  n_100 = generateSpatialCoupltesWithTime(n = 100),
  n_200 = generateSpatialCoupltesWithTime(n = 200),
  n_500 = generateSpatialCoupltesWithTime(n = 500),
  n_700 = generateSpatialCoupltesWithTime(n = 700),
  n_1000 = generateSpatialCoupltesWithTime(n = 1000),
  n_2000 = generateSpatialCoupltesWithTime(n = 2000),
  n_5000 = generateSpatialCoupltesWithTime(n = 5000),
  n_10000 = generateSpatialCoupltesWithTime(n = 10000),
  #n_50000 = generateSpatialCoupltesWithTime(n = 50000),
  times = 10,
  unit = "s")

theme_set(theme_light())

# TODO plot smoothing line
microbenchmark:::autoplot.microbenchmark(coupling_time_pairwise)
microbenchmark:::boxplot.microbenchmark(coupling_time_pairwise)



# as_tibble(coupling_time_pairwise) %>%
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



Coupling_Times <- coupling_times %>%  group_by(expr) %>%  summarise(median(time/1000))
Coupling_Times['ops'] <- c(100,500,1000,5000,10000, 50000, 100000)
names(Coupling_Times)[2] <- 'median' ;

ggplot(Coupling_Times,aes(ops,median),size = 1.08) +
  geom_point(aes(colour = 'data curve'),size=1.1) +
  geom_line(aes(colour = 'data curve')) +
  stat_function(fun = nlogn,size = 1.08,aes(colour = 'n*logn curve')) +
  labs(x = 'length of array',y = 'median of computing time [miliseconds]')



## computing times full likelihood ----
library(spdep)
library(sp)
library(spatialreg)

generateSpatialData <- function(n, phi = 0.8) {

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

  return(data)

}


produceSemEstimates <- function(n) {

  data = generateSpatialData(n)

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


coupling_times_fulllikelihood = microbenchmark(
  n_100 = produceSemEstimates(n = 100),
  n_200 = produceSemEstimates(n = 200),
  n_500 = produceSemEstimates(n = 500),
  n_700 = produceSemEstimates(n = 700),
  n_1000 = produceSemEstimates(n = 1000),
  n_2000 = produceSemEstimates(n = 2000),
  n_5000 = produceSemEstimates(n = 5000),
  n_10000 = produceSemEstimates(n = 10000),
  #n_100000 = produceSemEstimates(simdata = 50000),
  times = 10,
  unit = "s")

theme_set(theme_light())


microbenchmark:::autoplot.microbenchmark(coupling_times_fulllikelihood)


library(scales)
library(ggetho)
library(patchwork)
library(ggplot2)

coupling_times_fulllikelihood = read_rds(here("data", "coupling_times_fulllikelihood.rds"))


ct_likelihood = ggplot(tibble(coupling_times_fulllikelihood) %>%
         mutate(
           expr = as.numeric(str_remove(expr,"n_" )),
           time = lubridate::dmilliseconds(time/100000)
           ), aes(x = expr, y =  time)) +
  geom_point(alpha = 3) +
  geom_smooth(method = "glm", se = TRUE,  formula = y ~ x, method.args = list(family = gaussian(link = 'log')), alpha = 3) +
  scale_y_seconds(limits = c(0, 30000)) +
  labs(
    title = "Full Likelihood",
    x = "")

coupling_time_pairwise = read_rds(here("data", "coupling_time_pairwise.rds"))

ct_pairwise = ggplot(tibble(coupling_time_pairwise) %>%
         mutate(
           expr = as.numeric(str_remove(expr,"n_" )),
           time = lubridate::dmilliseconds(time)
         ), aes(x = expr, y =  time/100000)) +
  geom_point() +
  geom_smooth(method = "glm", se = TRUE,  formula = y ~ x, method.args = list(family = gaussian(link = 'log')), alpha = 3) +
  scale_y_seconds(limits = c(0,30000)) +
  labs(
    title = "Pairwise",
    x = ""
       )


tot = ct_likelihood + ct_pairwise

ggsave(tot, filename =here("img", "plots", "ct_loglig_vs_pairs.png"))


ggplot(smooth_line_data) +
  geom_point(aes(x = mean_time/100000, y = expr)) +
  geom_smooth(stat = "smooth", se = TRUE) +
  theme_light()





microbenchmark:::boxplot.microbenchmark(coupling_times_fulllikelihood)





## pathchwork of computational times


comput_complexity = bind_rows( as_tibble(coupling_times_fulllikelihood) %>%
  mutate(type = "fullLik"),
as_tibble(coupling_time_pairwise) %>%
  mutate(type  = "pairwise")
)

ct_on_same_scale = ggplot(tibble(comput_complexity) %>%
         mutate(
           expr = as.numeric(str_remove(expr,"n_" )),
           time = lubridate::dmilliseconds(time)
         ), aes(x = expr, y =  time/100000, col =  type)) +
  geom_point() +
  geom_smooth(method = "glm", se = TRUE,  formula = y ~ x, method.args = list(family = gaussian(link = 'log')), alpha = 3) +
  scale_y_seconds(limits = c(0,30000)) +
  labs(
    title = "PL vs. FL",
    x = ""
  )


ggsave(plot = ct_on_same_scaele,filename = here("img","plots" , "ct_on_same_scale.pmg"))
