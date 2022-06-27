# This script extract the results of the clustering


# Preamble ----------------------------------------------------------------

library(rstan)
library(tidyverse)
library(sf)
library(Cairo)

# Parse options
opt_list <- list(
  optparse::make_option(opt_str = c("-k", "--K"), type = "numeric",
                        default = 4, help = "n_groups"),
  optparse::make_option(opt_str = c("-r", "--redo"), type = "logical",
                        default = F, help = "redo"),
  optparse::make_option(opt_str = c("-s", "--sigma"), type = "numeric",
                        default = 1, help = "sigma"),
  optparse::make_option(opt_str = c("-t", "--tau_theta"), type = "numeric",
                        default = 1, help = "sigma")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

dir.create("figures")

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

fit_file <- glue::glue("generated_data/output_clustering/kmeans_{opt$K}_clusters_sigma{opt$sigma}_tau{opt$tau_theta}.rds")

k <- as.numeric(opt$K)  
out_file <- str_replace(fit_file, "kmeans", "kmeans_output")

if (!file.exists(out_file) | opt$redo) {
  
  print(out_file)
  stan_fit <- try(readRDS(fit_file))
  
  if(!inherits(stan_fit, "try-error")) {
    
    # Plot centers
    mus <- summary(stan_fit, pars = "mus")
    
    mus_s <- map_df(1:4, function(x) {
      mus$c_summary[,,x] %>% 
        as_tibble() %>% 
        mutate(chain = x,
               param = rownames(mus$c_summary[,,x]),
               cluster = str_extract(param, "(?<=\\[)[0-9]+"),
               month = str_extract(param, "(?<=\\,)[0-9]+") %>% 
                 as.numeric())
    }) %>% 
      mutate(cluster = as.numeric(cluster))
    
    p_mus <- ggplot(mus_s, aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, 
                               y = month, color = factor(chain))) +
      geom_errorbarh(height = 0) +
      geom_point() +
      coord_flip() +
      facet_wrap(~cluster) +
      theme_bw()
    
    ggsave(plot = p_mus, 
           filename = glue::glue("figures/plot_mus_{opt$K}_sigma{opt$sigma}_tau{opt$tau_theta}.png"), 
           width = 6,
           height = 4, 
           dpi = 300,
           type = "cairo")
    
    
    # Get membership probs
    probs <- map(1:4, function(x) {
      summary(stan_fit, pars = "probs")$c_summary[,1,x] %>% 
        matrix(ncol = k, byrow = T)
    })
    
    clusters <- map(probs, ~ apply(., 1, which.max)) %>% 
      bind_cols() %>% 
      magrittr::set_colnames(str_c("chain", 1:4)) %>% 
      mutate(gid = best_seas_data$gid)
    
    max_probs <- map(probs, function(x) {
      best_ind <- apply(x, 1, which.max)
      lapply(1:length(best_ind), function(y) x[y, best_ind[y]]) %>% 
        do.call(rbind, .)
    }) %>% 
      bind_cols() %>% 
      mutate(gid = best_seas_data$gid)
    
    colnames(max_probs) <- c(str_c("chain", 1:4), "gid")
    
    sf_objects2 <- inner_join(sf_objects, 
                              clusters %>% 
                                pivot_longer(cols = contains("chain"),
                                             names_to = "chain",
                                             values_to = "cluster")) %>% 
      inner_join( max_probs %>% 
                    pivot_longer(cols = contains("chain"),
                                 names_to = "chain",
                                 values_to = "max_prob")) %>% 
      mutate(max_prob_cat = cut(max_prob, c(0, .25, .5, .75, 1)))
    
    p_clusters <- ggplot(sf_objects2, aes(fill = factor(cluster), alpha = max_prob_cat)) +
      geom_sf(aes(geometry = geometry), size = .1, color = "white") +
      facet_wrap(~chain) +
      scale_alpha_discrete(range = c(.5, 1)) +
      ggthemes::theme_map() 
    
    ggsave(p_clusters, filename = glue::glue("figures/plot_clusters_{opt$K}_sigma{opt$sigma}_tau{opt$tau_theta}.png"), 
           width = 7, 
           height = 6,
           dpi = 400,
           type = "cairo")
    
    # Find cluster matching based on ids
    dist_mat <- array(dim = c(k, k , 4))
    
    for (m in 1:4) {
      for (i in 1:k) {
        for (j in 1:k) {
          
          ref <- sf_objects2 %>% 
            dplyr::filter(chain == "chain1", cluster == i) %>% 
            pull(gid)
          
          comp <- sf_objects2 %>% 
            dplyr::filter(str_detect(chain, as.character(m)), cluster == j) %>% 
            pull(gid)
          
          dist_mat[i, j, m] <- length(setdiff(ref, comp))
        }
      }
    }
    
    # Use first chain as ref
    relabel_df <- map_df(1:4, function(x) {
      tibble(chain = x,
             cluster = 1:k,
             relabeled_cluster = map_dbl(1:k, ~ which.min(dist_mat[., , x])))
    }) 
    
    # Recompute membership probs
    probs <- map(1:4, function(x) {
      ps <- rstan::extract(stan_fit, pars = "probs", permuted = F)[, x, ] 
      ps_mat <- reticulate::array_reshape(ps, c(nrow(ps), k, ncol(ps)/k), order = c("C"))
      
      relabel_dict <- dplyr::filter(relabel_df, chain == x) %>% 
        arrange(cluster) %>% 
        pull(relabeled_cluster)
      
      ps_mat <- ps_mat[, relabel_dict, ]
      
      return(ps_mat)
    })
    
    # Compute mean prob
    new_probs <- map_df(1:dim(probs[[1]])[3], function(x) {
      prob_means <- map(1:4, function(y) {
        probs[[y]][, , x]
      }) %>% 
        do.call(rbind, .) %>% 
        apply(2, mean)
      
      tibble(gid = best_seas_data$gid[x],
             cluster = 1:length(prob_means),
             prob = prob_means)
    })
    
    clusters2 <- new_probs %>% 
      group_by(gid) %>% 
      summarise(cluster = cluster[which.max(prob)])
    
    max_prob2 <- new_probs %>% 
      group_by(gid) %>% 
      summarise(max_prob = max(prob))
    
    sf_objects3 <- inner_join(sf_objects, clusters2) %>% 
      inner_join(max_prob2) %>% 
      mutate(max_prob_cat = cut(max_prob, c(0, .25, .5, .75, 1)))
    
    p_clusters_single <- ggplot(sf_objects3, aes(fill = factor(cluster), alpha = max_prob_cat)) +
      geom_sf(aes(geometry = geometry), size = .1, color = "white") +
      scale_alpha_discrete(range = c(.5, 1)) +
      ggthemes::theme_map() 
    
    ggsave(p_clusters_single, filename = glue::glue("figures/plot_clusters_{opt$K}_single_sigma{opt$sigma}_tau{opt$tau_theta}.png"),
           width = 7, 
           height = 6,
           dpi = 400,
           type = "cairo")
    
    
    
    # Recompute posteriors of the group means
    mu_samples <- map(1:4, function(x) {
      mus <- rstan::extract(stan_fit, pars = "mus", permuted = F)[, x, ] 
      mus_mat <- reticulate::array_reshape(mus, c(nrow(mus), ncol(mus)/k, k), order = c("C"))
      
      relabel_dict <- dplyr::filter(relabel_df, chain == x) %>% 
        arrange(cluster) %>% 
        pull(relabeled_cluster)
      
      mus_mat <- mus_mat[, , relabel_dict]
      
      return(mus_mat)
    })
    
    mu_stats <- map_df(1:12, function(x) {
      samples <- map(1:4, function(y) {
        mu_samples[[y]][,x,]
      }) %>% 
        do.call(what = rbind)
      
      map_df(1:k, function(i) {
        tibble(cluster = i,
               month = x,
               mean = mean(samples[, i]),
               q025 = quantile(samples[, i], .025),
               q975 = quantile(samples[, i], .975))
      })
    })
    
    res <- list(k = k, 
                loo = loo::loo(loo::extract_log_lik(stan_fit)),
                sf_object = sf_objects3,
                mu_stats = mu_stats)
    
    saveRDS(res, file = out_file)
  }
}
