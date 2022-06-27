
computeAdjacencySimple <- function(cntry.sf) {
  # extract the adjacency matrix
  adj_dat <- nb2graph(poly2nb(cntry.sf))
  
  N <- adj_dat$N;    # number of spatial units at the level of spatial interactions
  node1 <- adj_dat$node1; # "origin" node list
  node2 <- adj_dat$node2; # "destination" node
  N_edges <- adj_dat$N_edges; # number of edges
  
  
  nn_mat <- Matrix::sparseMatrix(i = node1, j = node2, x = 1, symmetric = T)
  isolated_vertices <- which(Matrix::rowSums(nn_mat) == 0)
  
  # Create graph to extract disconnected islands
  ng <- igraph::graph_from_adjacency_matrix(nn_mat)
  # Get the clusters
  ng_cl <- igraph::clusters(ng)
  
  if (ng_cl$no > 1) {
    cat("Found", ng_cl$no, "clusters of sizes {", paste(ng_cl$csize, collapse = ","), "}  adding edges from islands to mainland. \n")
    
    cluster_ids <- seq_len(ng_cl$no)
    mainland <- which(ng_cl$csize == max(ng_cl$csize))[1]
    mainland_ids <- which(ng_cl$membership == mainland)
    
    # Loop over island
    smooth_centroids <- sf::st_geometry(sf::st_centroid(cntry.sf))
    for (i in cluster_ids[-mainland]) {
      island_ids <- which(ng_cl$membership == i)
      
      # Workaround for MARCC since lwgeom is cannot be installed
      st_dist <- function(x, y) {
        cx <- st_coordinates(x)
        cy <- st_coordinates(y)
        distmat <- apply(cx, 1, function(j) apply(cy, 1, function(i)
          sqrt((i[2] - j[1])^2 + (i[1] - j[1])^2)))
        distmat
      }
      # find closest mainland pixels (n_pix_islands x n_pix_mainland matrix)
      dist_to_main <- st_dist(smooth_centroids[island_ids], smooth_centroids[mainland_ids])
      # get nearest ids for each island pixel
      nearest_main_ids <- apply(dist_to_main, 1, function(x) which(x == min(x))[1])
      nearest_dist <- unlist(mapply(x = 1:length(island_ids), y = nearest_main_ids, function(x, y) dist_to_main[x, y]))
      # get overall nearest mainland pixel
      nearest_isl_id <- which(nearest_dist == min(nearest_dist))[1]
      # connect the nearest island to the mainland (symetry)
      nn_mat[island_ids[nearest_isl_id], mainland_ids[nearest_main_ids[nearest_isl_id] ] ] <- 1
      nn_mat[mainland_ids[nearest_main_ids[nearest_isl_id] ], island_ids[nearest_isl_id] ] <- 1
    }
    
    # Check that everything is connected now
    if (igraph::clusters(igraph::graph_from_adjacency_matrix(nn_mat))$no > 1) {
      print(unique(igraph::clusters(igraph::graph_from_adjacency_matrix(nn_mat))$no > 1))
      stop("Something went wrong with island connection.")
    } else {
      cat("Done island connection\n")
    }
    
    # Re-extract adjacency list
    adj_dat2 <- mat2listw(nn_mat) %>% listw2sn() 
    
    unique_pairs <- adj_dat2[adj_dat2$from < adj_dat2$to, ]
    
    N <- length(unique(c(adj_dat2$from, adj_dat2$to)));    # number of spatial units at the level of spatial interactions
    if (adj_dat$N != N)
      stop("Something went wrong in re-connection")
    
    node1 <- unique_pairs$from # "origin" node list
    node2 <- unique_pairs$to; # "destination" node
    N_edges <- nrow(unique_pairs); # number of edges
  }
  
  res <- list(
    sf_object = cntry.sf,
    adj_list = list(N = N,
                    node1 = node1,
                    node2 = node2,
                    N_edges = N_edges)
  )
  
  return(res)
}