# Load required libraries
library(sp) # for spatial objects
library(igraph) # for graph functions
library(coda) # for MCMC output processing


# Define functions for clustering and threshold-based coupling
kmeans_cluster <- function(coords, n){
  # Perform k-means clustering of spatial coordinates
  kmeans_obj <- kmeans(coords, n)
  # Assign cluster labels to each spatial point
  clusters <- kmeans_obj$cluster
  return(clusters)
}

couplets_threshold <- function(coords, values, dist_threshold = 500, corr_threshold = 0.5){
  # Identify candidate couplets based on distance and correlation thresholds
  n <- length(coords)
  candidate_pairs <- t(combn(n, 2))
  candidate_distances <- apply(candidate_pairs, 1, function(x) dist(rbind(coords[x[1],], coords[x[2],])))
  candidate_correlations <- apply(candidate_pairs, 1, function(x) cor(rbind(values[x[1]], values[x[2]])))
  selected_pairs <- candidate_pairs[candidate_distances <= dist_threshold & abs(candidate_correlations) >= corr_threshold,]
  return(selected_pairs)
}



# Define function for optimization of pairwise likelihood
optimize_likelihood <- function(coords, values, selected_pairs, spatial_model, niter = 1000, thinning = 10){
  # Set up the likelihood function to optimize
  likelihood_fun <- function(theta, coords, values, selected_pairs, spatial_model){
    beta <- theta[1]
    tau2 <- exp(theta[2])
    sigma2 <- exp(theta[3])
    A <- sparse.adjacency(graph(selected_pairs), weighted = TRUE)
    X <- spdep::spMtx(cbind(1, coords))
    V <- diag(1/(sigma2 * spatial_model))
    XtVX <- t(X) %*% V %*% X
    XtVy <- t(X) %*% V %*% values
    invcov <- XtVX + tau2 * A
    logdet <- log(det(invcov))
    coef <- solve(invcov, XtVy)
    quad <- t(values) %*% V %*% values - t(values) %*% V %*% X %*% coef -
      t(coef) %*% t(X) %*% V %*% values + t(coef) %*% XtVX %*% coef
    likelihood <- -0.5 * (logdet + quad + length(values) * log(2*pi))
    return(likelihood)
  }
  # Set up the prior distributions for the parameters
  prior <- function(theta){
    beta <- theta[1]
    tau2 <- theta[2]
    sigma2 <- theta[3]
    lprior_beta <- dnorm(beta, 0, 1, log = TRUE)
    lprior_tau2 <- dgamma(tau2, 0.01, 0.01, log = TRUE)
    lprior_sigma2 <- dunif(sigma2, 0, 1000, log = TRUE)
    lprior <- lprior_beta + lprior_tau2 + lprior_sigma2
    return(lprior)
  }
  # Set up the MCMC chain and run it
  theta_start <- c(0, 0, 1)
  theta_chain <- mcmc(rep(NA, length(theta_start)), niter = niter, thin = thinning)
  accept_count <- rep(0, length(theta_start))
  theta_chain[1,] <- theta_start
  for(i in 2:niter){
    # Propose a new set of parameter values
    theta_prop <- rnorm(length(theta_start), mean = theta_chain[i-1,], sd = c(0.1, 0.1, 0.1))
    lprior_prop <- prior(theta_prop)
    # Compute the likelihood for the proposed values
    ll_prop <- likelihood_fun(theta_prop, coords, values, selected_pairs, spatial_model)
    # Compute the log-acceptance probability
    log_accept_prob <- ll_prop + lprior_prop - ll_chain[i-1] - lprior_chain[i-1]
    # Accept or reject the proposal
    if(log(runif(1)) < log_accept_prob){
      theta_chain[i,] <- theta_prop
      ll_chain[i] <- ll_prop
      lprior_chain[i] <- lprior_prop
      accept_count <- accept_count + 1
    } else{
      theta_chain[i,] <- theta_chain[i-1,]
      ll_chain[i] <- ll_chain[i-1]
      lprior_chain[i] <- lprior_chain[i-1]
    }
  }
  # Discard burn-in and thin the chain
  burnin <- niter/10
  thinned_chain <- as.matrix(as.mcmc(theta_chain[(burnin+1):niter,]))
  return(list(thinned_chain, accept_count))
}

# Define function for overall algorithm
coupled_likelihood <- function(coords, values, n_clusters, dist_threshold = 500, corr_threshold = 0.5,
                               niter = 1000, thinning = 10, spatial_model = 1){
  # Cluster spatial coordinates
  clusters <- kmeans_cluster(coords, n_clusters)
  # Initialize objects for output
  selected_pairs <- matrix(NA, nrow = 0, ncol = 2)
  likelihood_chain <- matrix(NA, nrow = 0, ncol = 3)
  accept_counts <- matrix(NA, nrow = 0, ncol = 3)
  # Loop over clusters and select pair candidates, optimize likelihood, and collect results
  for(clust in 1:n_clusters){
    clust_bool <- clusters == clust
    clust_coords <- coords[clust_bool,]
    clust_values <- values[clust_bool]
    clust_pairs <- couplets_threshold(clust_coords, clust_values, dist_threshold, corr_threshold)
    if(nrow(clust_pairs) > 0){
      clust_likelihood <- optimize_likelihood(clust_coords, clust_values, clust_pairs, spatial_model, niter, thinning)
      selected_pairs <- rbind(selected_pairs, clust_pairs)
      likelihood_chain <- rbind(likelihood_chain, clust_likelihood[[1]])
      accept_counts <- rbind(accept_counts, clust_likelihood[[2]])
    }
  }
  # Combine results across all clusters
  spatial_likelihood <- likelihood_chain[,1]
  tau2_chain <- exp(likelihood_chain[,2])
  sigma2_chain <- exp(likelihood_chain[,3])
  mean_accept_count <- mean(accept_counts, na.rm = TRUE)
  # Return selected pair indices, posterior chains and mean acceptance rate
  return(list(selected_pairs, spatial_likelihood, tau2_chain, sigma2_chain, mean_accept_count))
}



## debug!
for(clust in 1:8){
  clust_bool <- clusters == clust
  clust_coords <- coords[clust_bool,]
  clust_values <- values[clust_bool]
  clust_pairs <- couplets_threshold(clust_coords, clust_values, dist_threshold, corr_threshold)
  print(clust_pairs)
  if(nrow(clust_pairs) > 0){
    clust_likelihood <- optimize_likelihood(clust_coords, clust_values, clust_pairs, spatial_model, niter, thinning)
    selected_pairs <- rbind(selected_pairs, clust_pairs)
    likelihood_chain <- rbind(likelihood_chain, clust_likelihood[[1]])
    accept_counts <- rbind(accept_counts, clust_likelihood[[2]])
  }
}


## test it
library(spatstat)

# Generate simulated data on irregular grid
set.seed(123)
n <- 1000
window <- owin(c(0, 1), c(0, 1))
X <- runifpoint(n, window)
Y <- rnorm(n, mean = 0.5 * X$x + 0.5 * X$y, sd = 0.1)
coords <- as.data.frame(X)
values <- Y


library(sp)
library(igraph)
library(coda)

# Run the coupled likelihood algorithm
set.seed(123)
results <- coupled_likelihood(coords = coords, values = values, n_clusters = 10, dist_threshold = 0.1, corr_threshold = 0.5,
                              niter = 2000, thinning = 10, spatial_model = 1)
selected_pairs <- results[[1]]
spatial_likelihood <- results[[2]]
tau2_chain <- results[[3]]
sigma2_chain <- results[[4]]
mean_accept_count <- results[[5]]

library(ggplot2)

# Plot the simulated data and the selected pair candidates
data <- data.frame(coords, values)
pairs <- data.frame(coords[c(selected_pairs[,1], selected_pairs[,2]),], pair_index = rep(1:nrow(selected_pairs), each = 2))
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes(color = values), size = 3) +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_segment(data = pairs, aes(x = X1, y = Y1, xend = X2, yend = Y2), color = "black")


# Plot the posterior distributions of tau^2 and sigma^2
tau2_df <- data.frame(tau2_chain)
sigma2_df <- data.frame(sigma2_chain)
ggplot() +
  geom_density(data = tau2_df, aes(x = V1), fill = "blue", alpha = 0.5) +
  geom_density(data = sigma2_df, aes(x = V1), fill = "red", alpha = 0.5) +
  theme_bw() +
  labs(x = "Parameter value", y = "Density", color = "Parameter") +
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~variable, scales = "free")









