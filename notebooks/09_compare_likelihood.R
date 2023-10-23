library(Matrix)
compute_params <- function(data, couplets, sigma2=NULL){

  q <- nrow(couplets)
  X <- data$x
  Y <- data$y
  X1 <- couplets$X1_x
  X2 <- couplets$X2_x
  Y1 <- couplets$X1_y
  Y2 <- couplets$X2_y

  ## Compute sufficient statistics
  alpha1 <- sum(X1*X1) + sum(X2*X2)
  alpha2 <- sum(Y1*Y1) + sum(Y2*Y2)
  alpha3 <- sum(X1*Y1) + sum(X2*Y2)
  alpha4 <- sum(X1*Y2) + sum(X2*Y1)
  alpha5 <- sum(X1*X2)
  alpha6 <- sum(Y1*Y2)

  ## Initialize parameters
  beta_hat <- NULL
  psi_hat <- NULL
  sigma2_hat <- NULL

  for (i in 1:length(couplets)){
    x1 <- couplets$X1_x[i]
    x2 <- couplets$X2_x[i]
    y1 <- couplets$X1_y[i]
    y2 <- couplets$X2_y[i]
    e1 <- Y[x1] - beta_hat*X[x1] - psi_hat*(Y[x2] - beta_hat*X[x2])
    e2 <- Y[x2] - beta_hat*X[x2] - psi_hat*(Y[x1] - beta_hat*X[x1])
    if (is.null(sigma2)){
      psi_hat <- (e1*e2)/(q*sum((e1^2 + e2^2)/(1-psi_hat^2)))
      beta_hat <- (alpha3 - psi_hat*alpha4)/(alpha1 - 2*psi_hat*alpha5)
      sigma2_hat <- (alpha2 + beta_hat^2*alpha1 - 2*beta_hat*alpha3 - 2*psi_hat*alpha6 - 2*psi_hat*beta_hat^2*alpha5 + 2*psi_hat*beta_hatalpha4)/(2*q*(1-psi_hat^2))
    } else {
        psi_hat <- (e1*e2)/(qsigma2)
        beta_hat <- (alpha3 - psi_hat*alpha4)/(alpha1 - 2*psi_hat*alpha5)
        sigma2_hat <- sigma2 }
  }

  ## Compute score
  score <- c(beta_hat, sigma2_hat, psi_hat)
  fisher_info <- c(alpha1-2*psi_hat*alpha5,
                   alpha2-2*psi_hat*alpha6,
                   alpha3-alpha4*psi_hat,
                   alpha1*(1-psi_hat^2)/(2*sigma2_hat),
                   alpha2*(1-psi_hat^2)/(2*sigma2_hat),
                   alpha4*psi_hat/(sigma2_hat))

  names(score) <- c("beta_hat", "sigma2_hat", "psi_hat")
  names(fisher_info) <- c("alpha1", "alpha2", "alpha3", "I11", "I22", "I12")

  return(list(score=score, fisher_info=fisher_info))
  }




params <- compute_params(data, couplets_df)
params$score
params$fisher_info




## fifth iteratin
# Define the conditional distribution of each observation given its neighbors (spatial error model)
se_model <- function(x, y, psi, sigma2, beta) {
  eps <- y - beta * x
  density <- dnorm(y - psi * eps, mean = 0, sd = sqrt(sigma2), log = TRUE)
  return(density)
  }

# Define the log-likelihood function for a pair of adjacent observations
pairwise_likelihood <- function(point1, point2, psi, sigma2, beta) {
  x1 <- point1[1]
  y1 <- point1[2]
  x2 <- point2[1]
  y2 <- point2[2]
  alpha1 <- sum(c(x1^2, x2^2))
  alpha2 <- sum(c(y1^2, y2^2))
  alpha3 <- sum(c(x1y1, x2y2))
  alpha4 <- sum(c(x1y2, x2y1))
  alpha5 <- sum(x1x2)
  alpha6 <- sum(y1y2)

  # Compute the pairwise estimator parameters
   beta_hat <- (alpha3 - psi * alpha4) / (alpha1 - 2 * psi * alpha5)
   psi_hat <- (alpha6 - beta_hat * alpha4 + beta_hat^2 * alpha5) / (2 * (1 - psi^2) * nrow(nn$nn.idx))
   sigma2_hat <- (alpha2 + beta_hat^2 * alpha1 - 2 * beta_hat * alpha3 - 2 * psi_hat * alpha6 - 2 * psi_hat * beta_hat^2 * alpha5 + 2 * psi_hat * beta_hat * alpha4) / (2 * nrow(nn$nn.idx) * (1 - psi_hat^2))

  # Return the log-likelihood contribution for this pair of points
   if (psi_hat >= 1) { # Check for unstable estimates
     return(-Inf)
   } else {
       loglik <- se_model(x1, y1, psi_hat, sigma2_hat, beta_hat) + se_model(x2, y2, psi_hat, sigma2_hat, beta_hat)
       return(loglik)
       }
   }

  # Use the pairwise likelihood to estimate the model parameters for each pair of adjacent observations
  couplets_params <- matrix(NA, ncol = 3, nrow = length(couplets))
  for (i in seq_along(couplets)) {
    point1 <- couplets[[i]][[1]]
    point2 <- couplets[[i]][[2]]

  #. Maximize the pairwise likelihood to estimate model parameters
    pairwise_optim <- optim(c(0, 0, 0), pairwise_likelihood, point1 = point1, point2 = point2, control = list(fnscale = -1))

  # Save the estimated model parameters
     couplets_params[i, ] <- c(pairwise_optim$par, -pairwise_optim$value)
   }

  # Compute the full likelihood estimator parameters using the estimated pairwise parameters
  beta_full <- sum(couplets_params[, 1] * couplets_params[, 3]) / sum(couplets_params[, 3])
  psi_full <- sum(couplets_params[, 2] * couplets_params[, 3]) / sum(couplets_params[, 3])
  sigma2_full <- sum((couplets_params[, 3] * (couplets_params[, 1]^2 + couplets_params[, 2]^2 - 2*couplets_params[, 1]*beta_full - 2*couplets_params[, 2]*psi_full + 2*couplets_params[, 1]*couplets_params[, 2]*psi_full))/ (2*nrow(nn$nn.idx)(1-psi_full^2))) / sum(couplets_params[, 3])




## further guess
##


# compute pairwise likelihood parameters for each pair of adjacent locations
psi_beta_sigma2 <- lapply(couplets, function(cpl) {
  # extract the pairs of adjacent locations
  x_i <- cpl$X1_x_i
  x_l <- cpl$X2_x_i
  y_i <- cpl$X1_y
  y_l <- cpl$X2_y

  # compute sufficient statistics
  alpha1 <- sum(x_i ^ 2) + sum(x_l ^ 2)
  alpha2 <- sum(y_i ^ 2) + sum(y_l ^ 2)
  alpha3 <- sum(x_i * y_i) + sum(x_l * y_l)
  alpha4 <- sum(x_i * y_l) + sum(x_l * y_i)
  alpha5 <- sum(x_i * x_l)
  alpha6 <- sum(y_i * y_l)

  # compute maximum pairwise likelihood estimator for beta
  beta_mpl <-
    (alpha3 - psi_mpl * alpha4) / (alpha1 - 2 * psi_mpl * alpha5)

  # compute maximum pairwise likelihood estimator for sigma2
  eps <- (y_i - beta_mpl * x_i) - psi_mpl * (y_l - beta_mpl * x_l)
  sigma2_mpl <- sum(eps ^ 2) / (2 * length(eps) * (1 - psi_mpl ^ 2))

  # compute maximum pairwise likelihood estimator for psi
  psi_mpl <-
    (sum(y_i * y_l) - beta_mpl * alpha4 + beta_mpl ^ 2 * alpha5) / (length(eps) * sigma2_mpl)

  # return results
  return(list(
    psi = psi_mpl,
    beta = beta_mpl,
    sigma2 = sigma2_mpl
  ))
})

#combine pairwise likelihood parameters for all pairs into a single data frame
psi_beta_sigma2_df <- data.frame(do.call(rbind, psi_beta_sigma2))



