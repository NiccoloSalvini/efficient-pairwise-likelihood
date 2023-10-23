library(spData)
library(spdep) # for spatial regression functions
library(sp)
library(spatialreg)

data(baltimore)

# prepare coordinates, data, and proj4string
coords <- baltimore[ , c("X", "Y")]   # coordinates
data   <- baltimore[ , 1:15]          # data
crs    <- CRS("+proj=utm +zone=18 +datum=NAD83")# proj4string of coords

# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data,
                               proj4string = crs)

# Create spatial weights matrix using inverse distance
Wd <- dnearneigh(as.matrix(coordinates(spdf)), d1=0, d2=3000)

W.list <- nb2listw(Wd, style="B", zero.policy=TRUE) # style='B' for binary weights

# Define the response variable and the covariates:
y <- baltimore$PRICE
x <- cbind(baltimore$NROOM, baltimore$DWELL, baltimore$AGE, baltimore$NSTOR, baltimore$SQFT)

# Fit the spatial error model using maximum likelihood estimation
sem <- errorsarlm(y ~ x, data=baltimore, listw=W.list)
summary(sem)

# Retrieve the likelihood value and the estimated parameters
likelihood <- summary(sem)$logLik_lm.model; likelihood
beta <- coef(sem); beta
lambda <- sem$lambda; lambda
sig2 <- sem$s2; sig2



## fit SEM on simulated data (1000 obs singular covariare):
# function to generate spatial dataset and find minimum radius
find_mean_radius <- function(n = 1000, initial_radius = 0.01881959) {
  # create a random spatial dataset
  set.seed(1234)
  x_i <- runif(n)
  x_l <- runif(n)
  y <- runif(n)
  x <- runif(n)
  data <- data.frame(x_i = x_i, x_l = x_l, x = x , y = y)

  d <- as.numeric(as.matrix(dist(data))[,1])
  d[d<=0] <- NA

  # return the minimum distance between any two points
  mean_distance <- mean(d, na.rm = TRUE)
  radius <- ifelse(initial_radius > mean_distance, initial_radius, mean_distance)

  return(list("data" = data, "radius" = radius))
}

simul_data = find_mean_radius()$data


coords <- simul_data[ , c("x_i", "x_l")]   # coordinates
data   <- simul_data[ , 3:4]          # data
# crs    <- CRS("+proj=utm +zone=18 +datum=NAD83")# proj4string of coords

# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data)
                               #proj4string = crs)

# Create spatial weights matrix using inverse distance
Wd <- dnearneigh(as.matrix(coordinates(spdf)), d1=0, d2=3000)

W.list <- nb2listw(Wd, style="B", zero.policy=TRUE) # style='B' for binary weights

# Define the response variable and the covariates:
y <- simul_data$y
x <- cbind(simul_data$x)

# Fit the spatial error model using maximum likelihood estimation
sem <- errorsarlm(y ~ x, data=simul_data, listw=W.list)
summary(sem)

# Retrieve the likelihood value and the estimated parameters
likelihood <- summary(sem)$logLik_lm.model; likelihood
beta <- coef(sem); beta
rho <- sem$rho; rho
sig2 <- sem$s2; sig2



