library(sp)
library(rgdal)
library(rgeos)
library(gstat)

# Load the meuse dataset
data(meuse)
coordinates(meuse) <- c("x", "y")
proj4string(meuse) <- CRS("+init=epsg:28992")

# Build the k-d tree
library(RANN)
kdtree <- nn2(meuse@coords)

## we need to map kdtrees to knn object
# Form the knn object
n = 10
knn_obj <- list(nn = kdtree$nn.idx[, 2:n],
                np = length(meuse)
                dimension = ncol(meuse@coords),
                k = n,
                x = t(as.matrix(meuse@coords)))
class(knn_obj)= "knn"

knear6 <- knearneigh(meuse@coords, k=6)

# Form the spatial couples
n <- 10 # choose 10 nearest neighbors
couples <- knn2nb(knn_obj)

# Estimate the pairwise likelihood
model <- gstat(id="zinc", formula=zinc~1, data=meuse, nmax=1000)
covariance_matrix <- krige.conv(model=model, coords=meuse@coords, couples, covfct=model$krige.basis$covariance)
chol_decomp <- chol(covariance_matrix)
likelihood <- apply(chol_decomp, 1, function(x) prod(diag(x)))

# Evaluate the performance
# Here, we will assume that the true likelihood is unknown for this dataset
# However, you can compare the estimated likelihood with the results from other methods for comparison
## plot likelihood
library(ggplot2)

# Combine the estimated likelihood with the spatial data
likelihood_df <- data.frame(coordinates(meuse), likelihood)
colnames(likelihood_df) <- c("x", "y", "likelihood")
likelihood_df$zinc <- meuse$zinc

# Plot the estimated likelihood as a heatmap
ggplot(likelihood_df, aes(x=x, y=y)) +
  geom_tile(aes(fill=likelihood)) +
  scale_fill_gradient(low="white", high="blue") +
  geom_point(aes(color=zinc)) +
  coord_equal() +
  theme_bw()



## plot couplets
library(ggplot2)
n = 10
# Form the spatial couples using knn2nb function
couples <- knn2nb(knn_obj, k=n, include.all=TRUE, useC=FALSE)

# Convert spatial couples to data frame
couple_df <- data.frame(t(sapply(couples, function(x) c(x$x, x$y))), stringsAsFactors = F)
colnames(couple_df) <- c("x1", "y1", "x2", "y2")

# Merge the spatial data with couple information
couple_data <- merge(meuse, couple_df, by.x=c("x", "y"), by.y=c("x1", "y1"), all.x=TRUE)
couple_data <- merge(couple_data, meuse, by.x=c("x2", "y2"), by.y=c("x", "y"), all.x=TRUE)
colnames(couple_data) <- c("x1", "y1", "zinc1", "x2", "y2", "zinc2")

# Plot the couplets as line segments
ggplot(couple_data, aes(x=x1, y=y1, xend=x2, yend=y2)) +
  geom_segment(aes(color=abs(zinc1-zinc2)),
               size=0.2, alpha=0.6) +
  coord_equal() +
  theme_bw() +
  scale_color_gradient(low="darkgreen", high="red")

