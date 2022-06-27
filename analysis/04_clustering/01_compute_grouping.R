# This script computes the clustering of seasonalities

# Preamble ----------------------------------------------------------------

library(rstan)
library(tidyverse)
library(sf)
options(mc.cores = 4)
library(spdep)

source("analysis/nbgraph_utils.R")
source("analysis/utils.R")

# Parse options
opt_list <- list(
  optparse::make_option(opt_str = c("-k", "--K"), type = "numeric",
                        default = 4, help = "n_groups"),
  optparse::make_option(opt_str = c("-s", "--sigma"), type = "numeric",
                        default = 1, help = "sigma"),
  optparse::make_option(opt_str = c("-t", "--tau_theta"), type = "numeric",
                        default = 1, help = "sigma")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

dir.create("interm/clustering_int")

# Data --------------------------------------------------------------------

# Get best seasonality estimates
best_seas <- readRDS("generated_data/model_outputs/best_seas_mean_annual_incidence.rds")

# Get latitudes and longitudes
sf_objects <- readRDS("generated_data/clustering_sf_objects.rds") %>% 
  mutate(lat = st_coordinates(st_centroid(geometry))[,2],
         long = st_coordinates(st_centroid(geometry))[,1]) %>% 
  arrange(gid)

# Complete missing values for null models
best_seas <- best_seas %>% 
  filter(is.na(mean)) %>% 
  rowwise() %>% 
  group_modify(function(x, y) {
    filled <- inner_join(x %>% select(-gid), 
                         sf_objects, by = c("country", "run_level" = "gadm_lev")) %>% 
      mutate(mean = 0) 
    
    map_df(1:12, ~ mutate(filled, month = .))
  }) %>% 
  ungroup() %>% 
  bind_rows(best_seas %>% 
              filter(!is.na(mean))) %>% 
  arrange(gid)

# Extract
best_seas_data <- as_tibble(best_seas) %>%
  dplyr::filter(!is.na(month)) %>% 
  dplyr::select(gid, month, mean) %>% 
  mutate(month = month.abb[month]) %>% 
  pivot_wider(values_from = "mean",
              names_from = "month") 

# Data into matrix for kmeans
mat <- best_seas_data %>%
  inner_join(sf_objects %>% 
               as_tibble() %>% 
               dplyr::select(gid, lat, long) %>% 
               mutate(lat = scale(lat),
                      long = scale(long))) %>% 
  as_tibble() %>% 
  arrange(gid) %>% 
  dplyr::select(lat, long, one_of(month.abb)) %>% 
  as.matrix() 

# Admin unit adjacency
adj <- computeAdjacencySimple(sf_objects)

# Number of clusters
K <- opt$K

# Kmeans for init
init_kmeans <- kmeans(mat, K, iter.max = 100, nstart = 100)
ind_latlong <- which(colnames(mat) %in% c("lat", "long"))
centers <-  map(1:K, ~init_kmeans$centers[., -ind_latlong])
sigma <- opt$sigma
tau_theta <- opt$tau_theta

# Setup data
data <- list(N = nrow(mat), 
             N_edge = adj$adj_list$N_edges,
             node1 = adj$adj_list$node1,
             node2 = adj$adj_list$node2,
             D = 12, 
             K = K, 
             y = map(1:nrow(mat), ~ mat[., -ind_latlong]),
             prior_mu_location = centers,
             sigma = sigma,
             prior_mu_scale = .5,
             tau_theta = tau_theta
             )

# Run stan model ----------------------------------------------------------

model <- stan_model("analysis/stan/space_kmeans_simple.stan",
                    auto_write = F)
n_chains <- 4
stan_fit <- rstan::sampling(model, 
                            data = data,
                            init = map(1:n_chains,
                                       ~list(mus = map(centers, ~ rnorm(length(.), ., .05)))),
                            iter = 1250,
                            warmup = 250,
                            control = list(adapt_delta = .95,
                                           max_treedepth = 15,
                                           metric = "unit_e"),
                            chains = n_chains,
                            sample_file = str_glue("interm/clustering_int/kmeans_{opt$K}_clusters_sigma{sigma}_tau{tau_theta}"),
                            pars = c("soft_z"),
                            include = F)

saveRDS(stan_fit, file = glue::glue("generated_data/output_clustering/kmeans_{opt$K}_clusters_sigma{sigma}_tau{tau_theta}.rds"))
