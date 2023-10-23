library(dplyr)
library(nnls)
library(RANN)

# function to calculate pairwise likelihood between two points
calc_pairwise_likelihood <- function(point1, point2, model_type, model_parameters){
  distance <- sqrt(sum((point1 - point2)^2))
  if(model_type == "gaussian"){
    sigma <- model_parameters
    likelihood <- exp(-(distance^2)/(2*sigma^2))
  } else if(model_type == "matern"){
    nu <- model_parameters[1]
    rho <- model_parameters[2]
    likelihood <- (2^(1-nu)/gamma(nu))*(distance/rho)^nu*besselK(distance/rho, nu)
  } else if(model_type == "exponential"){
    rho <- model_parameters
    likelihood <- exp(-distance/rho)
  } else {
    stop("Invalid model type specified.")
  }
  return(likelihood)
}

# function to implement the Coupled KD-Tree algorithm
ckdt_algorithm <- function(data, k, model_type, model_parameters){
  nn_list <- nn2(data, data, k=k, searchtype="standard", treetype ="kd")
  pair_list <- lapply(seq_along(nn_list$nn.idx), function(i){
    data[i, ] %>%
      replicate(length(nn_list$nn.idx[i, ]), ., simplify = FALSE) %>%
      bind_rows(data[nn_list$nn.idx[i, ], ]) %>%
      filter(row_number() %% 2 == 1) %>%
      pull() %>%
      as.matrix() %>%
      matrix(ncol=2, byrow=TRUE) %>%
      apply(1, as.vector) %>%
      as.matrix()
  }) %>%
    do.call(rbind, .)
  return(sum(apply(pair_list, 1, function(pair){
    calc_pairwise_likelihood(point1=pair[1, ], point2=pair[2, ], model_type=model_type, model_parameters=model_parameters)
  })))
}

# specify data and model parameters
data <- matrix(rnorm(1000), ncol=2)
k <- 5
model_type <- "gaussian"
model_parameters <- 1

# run the algorithm and output results
ckdt_likelihood <- ckdt_algorithm(data=data, k=k, model_type=model_type, model_parameters=model_parameters)

full_likelihood <- 0 # implementation of full likelihood estimation method here

cat(paste0("Coupled KD-Tree likelihood: ", ckdt_likelihood))
cat("\n")
cat(paste0("Full likelihood: ", full_likelihood))


## prov aunwrapped
nn_list <- nn2(data, data, k=k, searchtype="standard", treetype ="kd")
pair_list <- lapply(seq_along(nn_list$nn.idx), function(i){
  data[i, ] %>%
    replicate(length(nn_list$nn.idx[i, ]), ., simplify = FALSE) %>%
    bind_rows(data[nn_list$nn.idx[i, ], ]) %>%
    filter(row_number() %% 2 == 1) %>%
    pull() %>%
    as.matrix() %>%
    matrix(ncol=2, byrow=TRUE) %>%
    apply(1, as.vector) %>%
    as.matrix()
}) %>%
  do.call(rbind, .)
return(sum(apply(pair_list, 1, function(pair){
  calc_pairwise_likelihood(point1=pair[1, ], point2=pair[2, ], model_type=model_type, model_parameters=model_parameters)
})))




# load required packages
library(sp)
library(FNN)

# define function to calculate pairwise likelihood
ckdt <- function(coords, data, K) {
  # build kd-tree from coordinates
  kd <- kd.tree(coords)
  # find k-nearest-neighbors pairs
  nn <- knn(kd, K + 1)$nn[, -1]
  # construct all non-redundant pairwise couples
  pairs <- t(combn(1:nrow(nn), 2))
  # compute pairwise distance and likelihood score
  dist <- apply(pairs, 1, function(x) spDistsN1(data[x, ], coords[nn[x, ], ]))
  score <- apply(dist, 2, function(x) sum(log(dnorm(x, 0, 1))))
  # sum all pairwise scores to obtain likelihood estimate
  return(sum(score))
}




## version2
library(RANN)
library(sp)

ckdt <- function(coords, data, K) {
  n <- nrow(coords)
  # compute nearest neighbors using RANN
  nn <- nn2(coords, k = K + 1, searchtype = "standard", treetype = "kd")$nn[, -1]
  # construct all non-redundant pairwise couples
  pairs <- t(combn(1:nrow(nn), 2))
  # compute pairwise distance and likelihood score
  dist <- apply(pairs, 1, function(x) spDistsN1(data[x, ], coords[nn[x, ], ]))
  score <- apply(dist, 2, function(x) sum(log(dnorm(x, 0, 1))))
  # sum all pairwise scores to obtain likelihood estimate
  sum(score)
}

# simulate data from Gaussian model
set.seed(123)
n <- 1000
coords <- matrix(runif(n*2), n, 2)
data <- rnorm(n)
ckdt(coords, data, 10) # test function with K = 10
# compare with full likelihood estimator
library(gstat)
g <- gstat(NULL, id="data", formula=data~1, data=data.frame(x=coords[,1], y=coords[,2]))
v <- vgm(psill=1, model="Gau", range=0.1)
fit <- gstat:::fit.variogram(g, v, fit.method=4)
fit$lik # full likelihood estimator result

# simulate data from Matérn model
set.seed(234)
n <- 5000
coords <- matrix(runif(n*2), n, 2)
data <- gstat::rmatern(n, model=list(nu=1.5), cov.pars=c(0.1, 0.2))
ckdt(coords, data, 20) # test function with K = 20
# compare with full likelihood estimator
g <- gstat(NULL, id="data", formula=data~1, data=data.frame(x=coords[,1], y=coords[,2]))
v <- vgm(psill=1, model="Mtn", range=0.1, nugget=0.01)
fit <- gstat:::fit.variogram(g, v, fit.method=4)
fit$lik # full likelihood estimator result

# simulate data from exponential model
set.seed(345)
n <- 10000
coords <- matrix(runif(n*2), n, 2)
data <- gstat::rexponential(coords, model=list(rho=0.01))
ckdt(coords, data, 30) # test function with K = 30
# compare with full likelihood estimator
g <- gstat(NULL, id="data", formula=data~1, data=data.frame(x=coords[,1], y=coords[,2]))
v <- vgm(psill=1, model="Exp", range=0.01, nugget=0.1)
fit <- gstat:::fit.variogram(g, v, fit.method=4)
fit$lik # full likelihood estimator result




x1 <- runif(100, 0, 2*pi)
x2 <- runif(100, 0,3)
DATA <- data.frame(x1, x2)
nearest <- nn2(DATA,DATA)


n <- nrow(coords)
# compute nearest neighbors using RANN
nn <- nn2(coords, k = 10 + 1, searchtype = "standard", treetype = "kd")
# construct all non-redundant pairwise couples
pairs <- t(combn(1:nrow(nn$nn.idx[, -1]), 2))
# compute pairwise distance and likelihood score
dist <- apply(pairs, 1, function(x) spDistsN1(data[x, ], coords[nn$nn.idx[x, ], ]))
score <- apply(dist, 2, function(x) sum(log(dnorm(x, 0, 1))))
# sum all pairwise scores to obtain likelihood estimate
sum(score)


pairs <- t(combn(1:nrow(nn$nn.idx[, 1]), 2))


library(ggplot2)
library(RANN)
# create a random spatial dataset
set.seed(123)
x <- runif(10000)
y <- runif(10000)
data <- data.frame(x = x, y = y)

# compute nearest neighbours using kdtree algorithm
nn <- nn2(data, data, treetype = "kd", searchtype = "standard", k = 2)

# extract the indexes of the nearest couplets

# plot the nearest couplets using ggplot2
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = data.frame(data[nn$nn.idx[, 1], ], data[nn$nn.idx[, 2], ]),
               aes(x = x, y = y, xend = x.1, yend = y.1),
               color = "red")

library(ggplot2)
library(RANN)

## find min radius
initial_radius = 0.08287174
( d <- as.numeric(as.matrix(dist(data))[,1]) )
d[d<=0] <- NA
which(d <= initial_radius)
min(d,na.rm=TRUE)


nn_rad <- nn2(data, data, treetype = "kd", searchtype = "radius", k= 2, radius = 0.08287174) #0.04287174

nn_rad$nn.idx[, 1]
data[nn_rad$nn.idx[, 1],]

nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1]
data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1],]

## filtrare per quelle diverse da 0
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = data.frame(data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1],], data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 2],]),
               aes(x = x, y = y, xend = x.1, yend = y.1),
               color = "red")




library(ggplot2)
library(RANN)

## find min radius
initial_radius = 0.08287174
( d <- as.numeric(as.matrix(dist(data))[,1]) )
d[d<=0] <- NA
which(d <= initial_radius)
min(d,na.rm=TRUE)


nn_rad <- nn2(data, data, treetype = "kd", k= 3)

nn_rad$nn.idx[, 1]
data[nn_rad$nn.idx[, 1],]

nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1]
data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1],]

## filtrare per quelle diverse da 0
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = data.frame(data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1],], data[nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 2],]),
               aes(x = x, y = y, xend = x.1, yend = y.1),
               color = "red")




## unique sulle colonne.
## sottomissione lo faccio e poi metto mio account
## spatial deep learning slides



#########################################
#########################################
### trova coppie non terzetti ###########
#########################################
#########################################

## adesso le due strade:
## 1 quella dove iterativamente ricalcoli l'alber
## 2 quella dove prendi la matrice e calcoli solo gli unici per colonna

## 1. prima strada
## sulla falsa riga del greedy.


## iteratively exploit kdtree to find couplets
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


## alternatively
library(ggplot2)
library(RANN)
# create a random spatial dataset
set.seed(123)
x <- runif(100)
y <- runif(100)
data <- data.frame(x = x, y = y)
# create a vector to keep track of which points have been paired
paired <- rep(FALSE, nrow(data))

# initialize empty output
pairs <- data.frame(x1 = numeric(0), y1 = numeric(0), x2 = numeric(0), y2 = numeric(0))

# loop until all points are paired
while (!all(paired)) {
  # find the nearest neighbours of unpaired points
  unpaired <- which(!paired) # indices of unpaired points
  nn <- nn2(data[unpaired,], data, treetype = "kd", searchtype = "standard", k = 2)
  # extract the closest neighbour for each unpaired point
  nn1 <- nn$nn.idx[,2]
  # find the indices of the closest couplets
  couplets <- cbind(unpaired, nn1)
  # extract the closest couplet
  closest <- couplets[which.min(nn$nn.dists[,1]),]
  # add the couplet to the output
  pairs <- rbind(pairs, data.frame(x1 = data[closest[1], "x"], y1 = data[closest[1], "y"],
                                   x2 = data[closest[2], "x"], y2 = data[closest[2], "y"]))
  # mark the points as paired
  paired[closest] <- TRUE
}

# plot the resulting pairs
ggplot(pairs, aes(x = x1, y = y1)) +
  geom_point() +
  geom_point(aes(x = x2, y = y2), color = "red") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray", alpha = 0.5)







## alternatively
library(ggplot2)
library(RANN)
# create a random spatial dataset
set.seed(123)
x <- runif(100)
y <- runif(100)
data <- data.frame(x = x, y = y)

# set up a vector to keep track of which points were already paired
paired <- rep(FALSE, nrow(data))

nn <- nn2(data, data, treetype = "kd", searchtype = "standard", k = 3)

# create an empty list to store the couplets
couplets <- list()

# loop through all the points
i <- 1
while(i <= nrow(data)) {
  # if the point is already paired, move on to the next one
  if(paired[i]) {
    i <- i + 1
    next
  }

  # extract the indexes of the nearest couplet
  j <- nn$nn.idx[i, 2]

  # if the second nearest neighbor is the same as the first (i.e. point i is at a corner of a rectangle), use the third nearest neighbor instead
  if(nn$nn.idx[j, 2] == i) {
    j <- nn$nn.idx[i, 3]
  }

  # if the second nearest neighbor is already paired, use the third nearest neighbor instead
  if(paired[j]) {
    j <- nn$nn.idx[i, 3]
  }

  # add the couplet to the list of couplets
  couplet <- rbind(data[i, ], data[j, ])
  couplet$id <- c(i, j)
  couplets[[length(couplets) + 1]] <- couplet

  # mark the two points as paired
  paired[i] <- TRUE
  paired[j] <- TRUE

  i <- i + 1
}

# plot the couplets
couplets_df <- bind_rows(couplets, .id = "couplet_id")

# separate the x, y, x.1, y.1 columns
couplets_df= couplets_df %>%
  pivot_longer(cols = -couplet_id, names_to = "coordinate", values_to = "value")  %>%
  pivot_wider(names_from = coordinate , values_from = value) %>%
  unnest_wider(c("x","y","id"), names_sep = ".")


# plot the couplets
ggplot(couplets_df) +
  geom_segment(aes(x = x.1, y = y.1, xend = x.2, yend = y.2), color = "red") +
  geom_point(data = data, aes(x, y)) +
  theme_bw()




## now try to be a little bit more clear
## DEFINITIVO
## richiede anche un po' di ottimizzazione
####################################################################################################################
####################################################################################################################
####################################################################################################################
#########################     QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII #######################
####################################################################################################################
####################################################################################################################

##
## lavora coi punti disaccoppiati che sono l'inverso degli accoppaiti
##

library(ggplot2)
library(RANN)
library(dplyr)
library(tidyr)
# create a random spatial dataset
set.seed(1234)
x <- runif(1000)
y <- runif(1000)
data <- data.frame(x = x, y = y)

# compute nearest neighbours using kdtree algorithm
nn <- nn2(data, data, treetype = "kd", searchtype = "standard", k = 2)

# initialize an empty list for paired points
couplets <- list()

# vector to keep track of the paired points
paired_points <- c()

# pair up points such that each point is only paired once
for (i in seq_len(nrow(nn$nn.idx))) {
  point1_index <- nn$nn.idx[i,1]
  point2_index <- nn$nn.idx[i,2]

  if (point1_index %in% paired_points || point2_index %in% paired_points) {
    next
  } else {
    # pair the two points and add them to the list and the vector of paired points
    point1 <- data[point1_index, ]
    point2 <- data[point2_index, ]

    couplets[[length(couplets) + 1]] <- list(point1, point2)
    paired_points <- c(paired_points,point1_index,point2_index)
  }
}
couplets_df = data.frame(do.call(rbind, couplets))  %>%
  unnest_wider(c("X1", "X2"), names_sep = "_")

# plot the coupled points using ggplot2
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = couplets_df,
               aes(x = X1_x, y = X1_y, xend = X2_x, yend = X2_y),
               color = "red")



## double test
data %>%
  rownames_to_column("id")  %>%
  ggplot() +
  geom_point(aes(coord_1, coord_2)) +
  geom_segment(data = poc[[1]],
               aes(x = el_1_coord_1, y = el_1_coord_2, xend = el_2_coord_1, yend = el_2_coord_2),
               color = "red")




## test
data_to_plot = data %>%
  as_tibble() %>%
  rownames_to_column("ID")

ggplot(data_to_plot, aes(coord_1, coord_2))+
  geom_point() +
  geom_segment(data = couplets,
               aes(x = el_1_coord_1, y = el_1_coord_2, xend = el_2_coord_1, yend = el_2_coord_2),
               color = "red") +
  geom_text(aes(label = ID), hjust=0, vjust=1.3)





###################################################################
###################################################################
##########       <<<<<< ADESSO PROVA COL BUFFER >>>>>     #########
###################################################################
###################################################################



## (funziona ma deve essere testato)
##

## adesso ci vuole il buffer, speriamo di risolvere col raggio!
library(ggplot2)
library(RANN)
library(dplyr)
library(tidyr)
# create a random spatial dataset
set.seed(1234)
x <- runif(1000)
y <- runif(1000)
data <- data.frame(x = x, y = y, )


## find min radius
initial_radius = 0.01881959
d <- as.numeric(as.matrix(dist(data)))
d[d<=0] <- NA
min(d,na.rm=TRUE)



# compute nearest neighbours using kdtree algorithm
nn <- nn2(data, data, treetype = "kd", searchtype = "radius", k= 2, radius = initial_radius)
# initialize an empty list for paired points
couplets <- list()

# vector to keep track of the paired points
paired_points <- c()

# pair up points such that each point is only paired once
for (i in seq_len(nrow(nn$nn.idx))) {
  point1_index <- nn$nn.idx[i,1]
  point2_index <- nn$nn.idx[i,2]

  if (point1_index %in% paired_points || point2_index %in% paired_points) {
    next
  } else {
    # pair the two points and add them to the list and the vector of paired points
    point1 <- data[point1_index, ]
    point2 <- data[point2_index, ]

    couplets[[length(couplets) + 1]] <- list(point1, point2)
    paired_points <- c(paired_points,point1_index,point2_index)
  }
}
couplets_df = data.frame(do.call(rbind, couplets))  %>%
  unnest_wider(c("X1", "X2"), names_sep = "_")

# plot the coupled points using ggplot2
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = couplets_df,
               aes(x = X1_x, y = X1_y, xend = X2_x, yend = X2_y),
               color = "red")


## apply algo to Baltimore dataset
## House sales price and characteristics for a spatial hedonic regression,
##  Baltimore, MD 1978. X,Y on Maryland grid, projection type unknown.
library(spData)
library(dplyr)
library(janitor)
library(sf)
library(ggplot2)
library(RANN)
baltimore_data =tibble(baltimore) %>%
  clean_names()

baltimore_sf <- baltimore_data %>% st_as_sf(., coords = c("x","y"))


## visualize baltimore data
ggplot(baltimore_sf) +
  geom_sf(aes(size = price, colour = -price))

baltimore_subset = select(baltimore_data, price, x,y)

## find min radius
# initial_radius = 0.01881959
( d <- as.numeric(as.matrix(dist(baltimore_subset))[,1]) )
min(d,na.rm=TRUE) # 14.84082
initial_radius = 39.68643
d[d<=0] <- NA
which(d <= initial_radius)

## forse vale la pena selezionare la media del raggio
mean(d,na.rm=TRUE) # 14.84082

# compute nearest neighbours using kdtree algorithm
nn <- nn2(baltimore_subset, baltimore_subset, treetype = "kd", searchtype = "radius", k= 2, radius = initial_radius)

# initialize an empty list for paired points
couplets <- list()

# vector to keep track of the paired points
paired_points <- c()

# pair up points such that each point is only paired once
for (i in seq_len(nrow(nn$nn.idx))) {
  point1_index <- nn$nn.idx[i,1]
  point2_index <- nn$nn.idx[i,2]

  if (point1_index %in% paired_points || point2_index %in% paired_points) {
    next
  } else {
    # pair the two points and add them to the list and the vector of paired points
    point1 <- baltimore_subset[point1_index, ]
    point2 <- baltimore_subset[point2_index, ]

    couplets[[length(couplets) + 1]] <- list(point1, point2)
    paired_points <- c(paired_points,point1_index,point2_index)
  }
}
couplets_df = data.frame(do.call(rbind, couplets))  %>%
  unnest_wider(c("X1", "X2"), names_sep = "_")

# plot the coupled points using ggplot2
ggplot(baltimore_subset ) +
  geom_point(aes(x, y, size  = price)) +
  geom_segment(data = couplets_df,
               aes(x = X1_x, y = X1_y, xend = X2_x, yend = X2_y),
               color = "red")





###################################################################
###################################################################
#####################      written a function    ##################
###################################################################
###################################################################


library(ggplot2)
library(RANN)
library(dplyr)
library(tidyr)

# function to generate spatial dataset and find minimum radius
find_min_radius <- function(n = 1000, initial_radius = 0.01881959) {
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
  min_distance <- min(d, na.rm = TRUE)
  radius <- ifelse(initial_radius > min_distance, initial_radius, min_distance)

  return(list("data" = data, "radius" = radius))
}


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


prova = find_mean_radius()
nn <- nn2(prova$data, prova$data, treetype = "kd", searchtype = "radius", k= 2, radius = prova$radius)
# initialize an empty list for paired points
couplets <- list()

# vector to keep track of the paired points
paired_points <- c()

# pair up points such that each point is only paired once
for (i in seq_len(nrow(nn$nn.idx))) {
  point1_index <- nn$nn.idx[i,1]
  point2_index <- nn$nn.idx[i,2]

  if (point1_index %in% paired_points || point2_index %in% paired_points) {
    next
  } else {
    # pair the two points and add them to the list and the vector of paired points
    point1 <- data[point1_index, ]
    point2 <- data[point2_index, ]

    couplets[[length(couplets) + 1]] <- list(point1, point2)
    paired_points <- c(paired_points,point1_index,point2_index)

  }
}
# function to pair up points within given radius
pair_points_within_radius <- function(data, radius) {

  # compute nearest neighbours using kdtree algorithm
  # need to look for the first two columns which are coordinates
  # (always place coordinates in the first two place)
  nn <- nn2(data, data, treetype = "kd", searchtype = "radius", k= 2, radius = radius)

  # initialize an empty list for paired points
  couplets <- list()

  # vector to keep track of the paired points
  paired_points <- c()

  # pair up points such that each point is only paired once
  for (i in seq_len(nrow(nn$nn.idx))) {
    point1_index <- nn$nn.idx[i,1]
    point2_index <- nn$nn.idx[i,2]

    if (point1_index %in% paired_points || point2_index %in% paired_points) {
      next
    } else {
      # pair the two points and add them to the list and the vector of paired points
      point1 <- data[point1_index, ]
      point2 <- data[point2_index, ]

      couplets[[length(couplets) + 1]] <- list(point1, point2)
      paired_points <- c(paired_points,point1_index,point2_index)
    }
  }

  # return a data frame with paired points
  return(data.frame(do.call(rbind, couplets)) %>% unnest_wider(c("X1", "X2"), names_sep = "_"))
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

pair_points <- function(x_i, x_l, y, x, initial_radius) {
  library(RANN)
  library(dplyr)
  library(tidyr)

  # create a data frame
  data <- data.frame(x_i, x_l, y, x)

  # find min radius
  ( d <- as.numeric(as.matrix(dist(data))[,1]) )
  d[d<=0] <- NA
  which(d <= initial_radius)
  mean(d, na.rm=TRUE)

  # compute nearest neighbours using kdtree algorithm
  nn <- nn2(data, data, treetype = "kd", searchtype = "radius", k= 2, radius = initial_radius)

  # initialize an empty list for paired points
  couplets <- list()

  # vector to keep track of the paired points
  paired_points <- c()

  # pair up points such that each point is only paired once
  for (i in seq_len(nrow(nn$nn.idx))) {
    point1_index <- nn$nn.idx[i,1]
    point2_index <- nn$nn.idx[i,2]

    if (point1_index %in% paired_points || point2_index %in% paired_points) {
      next
    } else {
      # pair the two points and add them to the list and the vector of paired points
      point1 <- data[point1_index, ]
      point2 <- data[point2_index, ]
      couplets[[length(couplets) + 1]] <- list(point1, point2)
      paired_points <- c(paired_points,point1_index,point2_index)
    }
  }

  # create a data frame from the list of paired points
  couplets_df <- data.frame(do.call(rbind, couplets)) %>%
    # unnest the paired point columns
    unnest_wider(c("X1", "X2"), names_sep = "_")

  # return the data frame with paired point coordinates and corresponding response and covariate values
  return(couplets_df[, c("x_i_1", "x_l_1", "x_i_2", "x_l_2", "y_1", "y_2", "x_1", "x_2")])
}


pair_points(x_i = runif(1000), x_l = runif(1000), y = runif(1000), x = runif(1000), initial_radius = 0.01881959)


