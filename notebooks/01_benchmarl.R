# Load libraries
library(spdep)
library(dplyr)
library(broom)

# Load data
data("ChicagoN")
# ChicagoN is a dataset included in the spdep package, consisting of 90000 observations of crime incidents in Chicago

# Create spatial weights matrix
W <- dnearneigh(coordinates(ChicagoN), d1 = 0, d2 = 250)

# Estimate SAR model using pairwise likelihood
sar_pl <- sar(plm(Crime ~ Income + Housing + Transit, data = ChicagoN), listW=W, SDEM=TRUE)

# Display summary results
tidy(sar_pl)
