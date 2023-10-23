
## try to find couplets with triangulation
library(ggplot2)
library(RANN)
library(deldir)

## find min radius
initial_radius = 0.08287174
( d <- as.numeric(as.matrix(dist(data))[,1]) )
d[d<=0] <- NA
which(d <= initial_radius)
min(d,na.rm=TRUE)


nn_rad <- nn2(data, data, treetype = "kd", searchtype = "radius", k= 2, radius = 0.1408147)

nn.idx <- nn_rad$nn.idx[, 1]

## create triangulation
tri <- deldir(data$x, data$y)

tri_edges <- tri$dirsgs[tri$delsgs == 0,]

## select only edges that connect nearest neighbors
nn_edges <- data.frame(from = nn.idx[nn_rad$nn.idx[,2] != 0], to = nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 2])
nn_edges <- nn_edges[order(nn_edges$from),]

## select only edges that form a couplet
couplet_edges <- nn_edges[match(nn_edges$from, tri_edges[,1]) == match(nn_edges$to, tri_edges[,2]),]

## plot the couplets
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = data.frame(data[couplet_edges$from,], data[couplet_edges$to,]),
               aes(x = x, y = y, xend = x.1, yend = y.1),
               color = "red")



# create a triangulation of the data
library(deldir)
tri <- deldir(data$x, data$y)

# extract the edges from the triangulation
edges <- tri$dirsgs

# find all edges that connect two non-boundary points
non_boundary_edges <- edges[rowSums(is.na(edges)) == 0, ]
couplet_edges <- non_boundary_edges[rowSums(non_boundary_edges %in% nn_rad$nn.idx[nn_rad$nn.idx[,2] != 0, 1]) == 2, ]

# plot the couplets as line segments
ggplot(data) +
  geom_point(aes(x, y)) +
  geom_segment(data = data.frame(couplet_edges %>% unique()), aes(x = x, y = y, xend = x.1, yend = y.1), color = "red")


library(spdep)
library(ggplot2)

# create example spatial dataset with 10 points
set.seed(123)
coords <- matrix(rnorm(20), ncol = 2)
spatialdata <- SpatialPointsDataFrame(coords, data.frame(ID = 1:10))

# calculate spatial weights matrix
W <- knn2nb(knearneigh(coords, k = 2))
W <- nb2listw(W)

# identify non-redundant couplets
nonred_couplets <- dnearneigh(d1 = coords, d2 = 0, longlat = FALSE, x = W)

# visualize non-redundant couplets
couplets_df <- data.frame(lapply(nonred_couplets, unlist))
ggplot(couplets_df, aes(x = x, y = y)) +
  geom_point(size = 3) +
  geom_segment(aes(x = x.1, y = y.1, xend = x, yend = y), alpha = 0.5) +
  theme_bw()





library(RANN)
library(igraph)

# Generate sample data
set.seed(123)
n <- 100
x <- runif(n, 0, 100)
y <- runif(n, 0, 100)
data <- cbind(x, y)

# Build kd-tree
kdtree <- nn2(data, data, k = 1, treetype =  "kd")
edges <- data.frame(from = integer(), to = integer(), type = character())

# Recursive function to add edges to `edges` data frame
add_edges <- function(kdtree, parent) {
  if (!is.leaf(kdtree)) {
    left_child <- add_vertices(1)
    right_child <- add_vertices(1)
    edges <<- bind_rows(edges, data.frame(from = parent, to = left_child, type = "vertical"),
                        data.frame(from = parent, to = right_child, type = "vertical"))
    add_edges(kdtree$vl, left_child)
    add_edges(kdtree$vr, right_child)
  } else {
    leaf_node <- add_vertices(1)
    edges <<- bind_rows(edges, data.frame(from = parent, to = leaf_node))
  }
}

# Create igraph object
g <- make_empty_graph()
add_edges(kdtree, 1)
g <- add_edges(g, edges)

# Plot igraph tree
plot(g, layout = layout_as_tree(g, root = vcount(g), flip.y = TRUE), vertex.shape = "none", edge.arrow.mode = 0,
     vertex.color = "lightblue", vertex.size = 5, edge.color = ifelse(E(g)$type == "vertical", "blue", "gray"))



# Load required libraries
library(RANN)
library(igraph)
library(ggplot2)

# Generate random dataset
set.seed(123)
data <- matrix(runif(20), ncol = 2)

# Build kdtree
kdtree <- RANN::nn2(data, k = 2, treetype  = "kd")

# Build igraph object from kdtree
tree_graph <- igraph::graph_from_data_frame(kdtree[["nn.idx"]])

# Plot the tree using igraph
plot(tree_graph, layout = layout_as_tree(tree_graph,root = 1) , vertex.shape = "none", edge.width = 1.5, vertex.size = 5)

# Build a ggplot object from kdtree
tree_df <- data.frame(from = rep(1:length(kdtree[["nn.idx"]]), each = 2), to = as.vector(kdtree[["nn.idx"]]), x = data[, 1], y = data[, 2])
tree_plot <- ggplot(data = tree_df, aes(x = x, y = y)) +
  geom_segment(aes(x = x[from], y = y[from], xend = x[to], yend = y[to]), size = 1.5) +
  geom_point(size = 5)

# Plot the tree using ggplot2
tree_plot





library(tripack)

# create example spatial dataset
x <- c(0, 1, 2, 0, 1, 2)
y <- c(0, 0, 0, 1, 1, 1)
tri <- tri.mesh(x, y)

# extract edge information from triangulation
edges <- tri$edges
edge.count <- nrow(edges)

# find edges shared by adjacent triangles
couplets <- matrix(0, ncol = 2, nrow = edge.count)
counter <- 0

for(i in 1:edge.count) {
  # check if edge is shared by more than one triangle
  shared.triangles <- tri$fnd.T[edgecount[i,1],] == tri$fnd.T[edgecount[i,2],]
  if(sum(shared.triangles) > 1) {
    counter <- counter + 1
    couplets[counter,] <- edges[i,]
  }
}

# remove unused rows from couplets matrix
couplets <- couplets[1:counter,]
