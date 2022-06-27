# Occurrence probability simulation -----------------------------------------

#' Compute posterior retrodictive checks
#'
#' @param country 
#' @param run_level 
#' @param identifier 
#' @param model 
#' @param redo 
#' @param time_left 
#' @param time_right 
#' @param redo_single 
#' @param nsim 
#' @param nsimbinom 
#' @param do_par 
#' @param n_cores 
#' @param what 
#' @param ... 
#'
#' @return
#'
computePPC <- function(country,
                       run_level,
                       identifier,
                       model = "all",
                       redo = NULL,
                       time_left = NULL,
                       time_right = NULL,
                       redo_single = F,
                       nsim = 200,
                       nsimbinom = 5,
                       do_par = T,
                       n_cores = parallel::detectCores()/2,
                       what = "coverage",
                       ...) {
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo,
                                 verbose = F)  
  }
  
  if (model == "best") {
    model <- getBestModel(country = country,
                          run_level = run_level,
                          identifier = identifier,
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = setGadmLev(country, run_level),
                          redo = redo)
  }  
  
  if (!(what %in% c("coverage", "monthly"))) {
    stop("Coverage not known")
  }
  
  runChecks(country = country,
            run_level = run_level,
            identifier = identifier,
            model = model)
  
  if (is.null(redo)) {
    redo <- !checkAllDone(
      country = country,
      run_level = run_level,
      identifier = identifier,
      time_left = time_left,
      time_right = time_right,
      suffix = "sample_stats",
      file_type = "rds")
  }
  
  ppc_file <- makeStdResName(country = country, 
                             run_level = run_level,
                             identifier = identifier, 
                             model = model,
                             suffix = "ppc_stats",
                             file_type = "rds")
  
  
  if (what == "coverage") {
    ppc_file <- str_replace(ppc_file, "\\.rds", "coverage.rds")
  }
  
  if (file.exists(ppc_file) & !redo) {
    
    res <- readRDS(ppc_file) %>% 
      mutate(country = country,
             run_level = run_level,
             model = model)
    
    return(res)  
  } 
  
  # Fraction of positive observations in each month
  predict_probs <- predictProbs(country = country,
                                run_level = run_level,
                                identifier = identifier,
                                model = model,
                                redo = redo,
                                redo_single = redo_single,
                                nsim = nsim,
                                nsimbinom = nsimbinom)
  
  if (is.null(predict_probs)) {
    return(NULL)
  }
  
  if (what == "monthly") {
    ppc_stats <-  predict_probs %>% 
      dplyr::filter(set == "observations") %>%
      dplyr::rename(obs_len = n_months) %>% 
      dplyr::group_by(month, lev, obs_len, identifier, variant, country) %>%
      dplyr::summarise(n_neg = sum(n_obs-n_pos),
                       n_pos = sum(n_pos),
                       samples_pos = list(map_dbl(1:(nsim*nsimbinom), function(x) {purrr::map_dbl(samples, ~ .[x]) %>% sum()})),
                       samples_neg = list(map_dbl(1:(nsim*nsimbinom), function(x) {purrr::map2_dbl(samples, n_obs, ~ .y - .x[x]) %>% sum()}))) %>%
      dplyr::mutate(mean_n_pos = map_dbl(samples_pos, ~mean(.)),
                    q025_n_pos = purrr::map_dbl(samples_pos, ~ quantile(., .025)),
                    q975_n_pos = purrr::map_dbl(samples_pos, ~ quantile(., .975)),
                    mean_n_neg = purrr::map_dbl(samples_neg, ~mean(.)),
                    q025_n_neg = purrr::map_dbl(samples_neg, ~ quantile(., .025)),
                    q975_n_neg = purrr::map_dbl(samples_neg, ~ quantile(., .975))) %>% 
      tidyr::pivot_longer(cols = contains("n_"),
                          names_to = "set",
                          values_to = "n") %>%
      dplyr::mutate(type = ifelse(str_detect(set, "pos"), "pos", "neg"),
                    what = map_chr(set, ~ str_split(., "(^|_)n_")[[1]][1] %>%
                                     str_remove("_") %>% 
                                     ifelse(. == "", "obs", .))) %>%
      dplyr::select(-set) %>% 
      tidyr::pivot_wider(values_from = "n",
                         names_from = "what") 
    
    
  } else {
    
    coverage <- map_df(c(.05, seq(.1, .9, by = .1), .95),
                       function(x) {
                         bounds <- widthToInt(x)
                         predict_probs %>% 
                           dplyr::filter(set == "observations") %>%
                           dplyr::mutate(obs_low = map_dbl(samples, ~quantile(., bounds[1])),
                                         obs_high = map_dbl(samples, ~quantile(., bounds[2])),
                                         width = x) %>% 
                           dplyr::select(obs, lev, n_pos, obs_low, obs_high, width)
                       })
    
    ppc_stats <- coverage %>% 
      mutate(lev = as.character(lev)) %>% 
      rbind(coverage %>% mutate(lev = "all")) %>% 
      dplyr::group_by(width, lev) %>%
      summarise(coverage = sum(n_pos >= obs_low & n_pos <= obs_high)/n()) %>% 
      mutate(country = country,
             run_level = run_level,
             model = model)
  }
  
  
  saveRDS(ppc_stats, file = ppc_file)
  
  
  return(ppc_stats)
}

#' Compute probabilities
#'
#' @param data 
#' @param init_par 
#' @param do_par 
#' @param n_cores 
#' @param variant 
#'
#' @return
#' 
computeProbs <- function(data, 
                         init_par, 
                         do_par, 
                         n_cores,
                         variant) {
  
  if (length(init_par) == 1){
    probs <- computeSingleProb(data = data, 
                               init_par = init_par,
                               variant = variant)
    return(probs)
  } 
  
  acomb <- function(...) abind::abind(..., along=1)
  
  
  if (do_par) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("computeSingleProb", "log_diff_exp", "log1m_exp", "modelProb", "modelProbNull"))
    n_chunks <- ceiling(length(init_par)/n_cores)
    
    out <- foreach(parset = itertools::ichunk(init_par, n_chunks),
                   .combine = acomb,
                   .inorder = F,
                   .packages = c("qgam", "foreach", "abind", "iterators", "stringr", "itertools")
    ) %dopar% {
      foreach(parsubset = parset,
              .combine = acomb) %do% {
                res <- computeSingleProb(data = data, 
                                         init_par = parsubset, 
                                         variant = variant) 
                out <- array(dim = c(1, dim(res)))
                out[1,,,] <- res
                out
              }
    }
    parallel::stopCluster(cl)
  } else {
    out <- foreach(parsubset = init_par,
                   i = icount(length(init_par)),
                   .combine = acomb) %do% {
                     res <- computeSingleProb(data = data, 
                                              init_par = parsubset, 
                                              variant = variant) 
                     out <- array(dim = c(1, dim(res)))
                     out[1,,,] <- res
                     cat("Done", i, "\n")
                     out
                   }
  }
  return(out)
}


#' Compute simulated probabilities
#' 
#' @description 
#' 
#' @param variant 
#' @param data 
#' @param init_par 
#'
#' @return
#'
computeSingleProb <- function(variant, 
                              data, 
                              init_par) {
  if (length(init_par) == 1)
    init_par <- init_par[[1]]
  
  # Determine whether the model is the speedup version or not
  model_version <- getModelVersion(data) 
  
  if (model_version == "old") {
    if (variant == "null") {
      modelProbNull(data, init_par)
    } else {
      modelProb(data, init_par)
    }
  } else if (model_version == "speedup") {
    if (variant == "null") {
      modelProbNullSpeedup(data, init_par)
    } else {
      modelProbSpeedup(data, init_par)
    }
  }
}


#' Compute mean seasonality
#'
#' @description computes mean seasonlity coefficients for each admin unit based
#'
#' @param betas_g1 mean of betas of group 1 (or of country if not mixture) (month, mean)
#' @param betas_g2 mean of betas of group 2 (if mixture)
#' @param lambdas mixture coefficients (admin_unit, mean)
#'
#' @return a dataframe with the mean betas for each admin unit
computeMeanSeasonality <- function(betas_g1, 
                                   betas_g2 = NULL, 
                                   lambdas = NULL) {
  if (is.null(betas_g2)) {
    return(betas_g1 %>%
             rename(betas = mean) %>% 
             mutate(admin_unit = "country"))
  }
  
  all_betas <- inner_join(betas_g1, betas_g2, by = "month", suffix = c(".g1", ".g2"))
  
  map_df(1:nrow(lambdas), function(x) {
    all_betas %>% 
      mutate(betas = mean.g1 * lambdas$mean[x] + mean.g2 * (1-lambdas$mean[x]),
             admin_unit = lambdas$admin_unit[x]) %>% 
      select(admin_unit, month, betas)
  })
}

#' Model probabilties null model
#'
#' @param data 
#' @param init_par 
#'
#' @return
#'
modelProbNull <- function(data,
                          init_par) {
  
  for (par in names(data)) {
    s <- str_c(par, "<- data$", par) 
    eval(expr = parse(text = s))
  }
  # Transformed
  
  log_unif <- -log(12)
  offset_years <- matrix(nrow = 12, ncol = K); 
  
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  for (i in 0:11) {
    for (k in 1:K) {
      if (do_offsets == 1) {
        offset_years[i+1, k] = ifelse(months[k] + i > 12, years[k] + 1, years[k]);
      } else {
        offset_years[i+1, k] = years[k];
      }
    }
  }
  
  for (par in names(init_par)) {
    s <- str_c(par, "<- init_par$", par) 
    eval(expr = parse(text = s))
  }
  
  # Transformed parameters
  sigma = 1/sqrt(tau); # overall marginal standard deviation
  prob_non_zero <- array(dim = c(N_offset_elems, 1, N_obs))  # probabilities of cholera at the upper spatial level (admin 1 or country)
  
  if (!exists("lambdas")) {
    lambdas <- rep(1, N_units)
  }
  
  # Compute log-liks
  sum_log1p <- rep(NA, N_obs_upper)
  sum_l <- rep(NA, N_obs_upper)
  for (z1 in 1:N_offset_elems) {
    if (N_obs_upper > 0) {
      for (i in 1:N_obs_upper) {
        sum_log1p[i] = 0;
        sum_l[i] = 0;
      }
      for (k in 1:N_map_upper) {
        if (map_2_units[ind_upper[k]] == 0) {
          # For country-level data loop over all units
          for (j in 1:N_units) {
            x = 0;
            x = x + xis[map_2_lev[ind_upper[k]]-1] + etas[offset_years[z1, ind_upper[k]], 1] + convolved_re[j] * sigma;
            
            sum_l[map_2_upper_obs[ind_upper[k]]] = sum_l[map_2_upper_obs[ind_upper[k]]] + x;
            sum_log1p[map_2_upper_obs[ind_upper[k]]] = sum_log1p[map_2_upper_obs[ind_upper[k]]] + qgam::log1pexp(-1 * x);
          }
        } else {
          x = 0;
          # For the rest of the data loop over subunits
          if (map_2_lev[ind_upper[k]] > 1) {
            # for upper-level units apply a reporting probability
            x = x + xis[map_2_lev[ind_upper[k]]-1];
          }
          x = x + etas[offset_years[z1, ind_upper[k]], 1] + convolved_re[map_2_units[ind_upper[k]]] * sigma;
          
          sum_l[map_2_upper_obs[ind_upper[k]]] = sum_l[map_2_upper_obs[ind_upper[k]]]  + x;
          sum_log1p[map_2_upper_obs[ind_upper[k]]] = sum_log1p[map_2_upper_obs[ind_upper[k]]] + qgam::log1pexp(-1 * x);
        }
      }
      
      for (i in 1:N_obs_upper) {
        if (is.nan(log_diff_exp(sum_log1p[i], -sum_l[i]))) {
          cat(i, sum_log1p[i], sum_l[i], "\n")
        }
        prob_non_zero[z1, 1, obs_ind_upper[i]] = log_diff_exp(sum_log1p[i], -sum_l[i]) + sum_l[i];
      }
    }
    # Set lower units
    prob_non_zero[z1, 1, obs_ind_lower] <- etas[offset_years[z1, ind_lower], 1] + convolved_re[map_2_units[ind_lower]] * sigma
  }
  return(prob_non_zero)
}

modelProb <- function(data,
                      init_par) {
  
  for (par in names(data)) {
    s <- str_c(par, "<- data$", par) 
    eval(expr = parse(text = s))
  }
  # Transformed
  
  log_unif <- -log(12)
  offset_years <- matrix(nrow = 12, ncol = K); 
  
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  if (N_groups > 1) {
    do_grouping = 1;
    N_betas_diff = 2;
    N_betas_g2 = 12;
    N_betas_other_g2 = 10;
    N_lambdas = N_units;
    
    if (N_eta_groups > 1) {
      ind_eta_grp2 = 2;
    } else {
      ind_eta_grp2 = 1;
    }
    
  } else  {
    do_grouping = 0;
    N_betas_diff = 0;
    N_betas_g2 = 0;
    N_betas_other_g2 = 0;
    N_lambdas = 0;
    ind_eta_grp2 = 1;
  }
  
  if(do_offsets == 1 && do_grouping == 1){
    N_offset_elems_g2 = 12;
    N_offset_elems = 1;
  } else {
    N_offset_elems_g2 = 1;
  }
  
  # cat("Running model do_offsets:", do_offsets, ", do_grouping:", do_grouping, "\n");
  # cat("N_ocffset_elems:", N_offset_elems, ", N_offset_elems_g2:", N_offset_elems_g2, "\n");
  
  
  for (i in 0:11) {
    for (k in 1:K) {
      if (do_offsets == 1) {
        offset_years[i+1, k] = ifelse(months[k] + i > 12, years[k] + 1, years[k]);
      } else {
        offset_years[i+1, k] = years[k];
      }
    }
  }
  
  for (par in names(init_par)) {
    s <- str_c(par, "<- init_par$", par) 
    eval(expr = parse(text = s))
  }
  
  # Transformed parameters
  sigma = 1/sqrt(tau); # overall marginal standard deviation
  prob_non_zero <- array(dim = c(N_offset_elems, N_offset_elems_g2, N_obs))  # probabilities of cholera at the upper spatial level (admin 1 or country)
  
  if (!exists("lambdas")) {
    lambdas <- rep(1, N_units)
  }
  
  # Compute log-liks
  sum_log1p <- rep(NA, N_obs_upper)
  sum_l <- rep(NA, N_obs_upper)
  for (z1 in 1:N_offset_elems) {
    for (z2 in 1:N_offset_elems_g2) {
      
      if (N_obs_upper > 0) {
        for (i in 1:N_obs_upper) {
          sum_log1p[i] = 0;
          sum_l[i] = 0;
        }
        for (k in 1:N_map_upper) {
          if (map_2_units[ind_upper[k]] == 0) {
            # For country-level data loop over all units
            for (j in 1:N_units) {
              x = 0;
              x = x + xis[map_2_lev[ind_upper[k]]-1] + lambdas[j] * (betas_g1[months[ind_upper[k]]] + etas[offset_years[z1, ind_upper[k]], 1]) + convolved_re[j] * sigma;
              
              if (do_grouping == 1) {
                x = x + (1-lambdas[j]) * (betas_g2[months[ind_upper[k]]] + etas[offset_years[z2, ind_upper[k]], ind_eta_grp2]);
              }
              sum_l[map_2_upper_obs[ind_upper[k]]] = sum_l[map_2_upper_obs[ind_upper[k]]] + x;
              sum_log1p[map_2_upper_obs[ind_upper[k]]] = sum_log1p[map_2_upper_obs[ind_upper[k]]] + qgam::log1pexp(-1 * x);
            }
          } else {
            x = 0;
            # For the rest of the data loop over subunits
            if (map_2_lev[ind_upper[k]] > 1) {
              # for upper-level units apply a reporting probability
              x = x + xis[map_2_lev[ind_upper[k]]-1];
            }
            x = x +  lambdas[map_2_units[ind_upper[k]]] * (betas_g1[months[ind_upper[k]]] + etas[offset_years[z1, ind_upper[k]], 1]) +
              convolved_re[map_2_units[ind_upper[k]]] * sigma;
            
            if (do_grouping == 1) {
              x = x + (1-lambdas[map_2_units[ind_upper[k]]]) * (betas_g2[months[ind_upper[k]]] + etas[offset_years[z2, ind_upper[k]], ind_eta_grp2]);
            }
            sum_l[map_2_upper_obs[ind_upper[k]]] = sum_l[map_2_upper_obs[ind_upper[k]]]  + x;
            sum_log1p[map_2_upper_obs[ind_upper[k]]] = sum_log1p[map_2_upper_obs[ind_upper[k]]] + qgam::log1pexp(-1 * x);
          }
        }
        
        for (i in 1:N_obs_upper) {
          if (is.nan(log_diff_exp(sum_log1p[i], -sum_l[i]))) {
            cat(i, sum_log1p[i], sum_l[i], "\n")
          }
          prob_non_zero[z1, z2, obs_ind_upper[i]] = log_diff_exp(sum_log1p[i], -sum_l[i]) + sum_l[i];
        }
      }
      # Set lower units
      prob_non_zero[z1, z2, obs_ind_lower] <- lambdas[map_2_units[ind_lower]] * (betas_g1[months[ind_lower]] + etas[offset_years[z1, ind_lower], 1]) +
        convolved_re[map_2_units[ind_lower]] * sigma
      
      if (do_grouping == 1) {
        prob_non_zero[z1, z2, obs_ind_lower] =  prob_non_zero[z1, z2, obs_ind_lower] + (1-lambdas[map_2_units[ind_lower]]) * (betas_g2[months[ind_lower]] + etas[offset_years[z2, ind_lower], ind_eta_grp2]);
      }
    }
  }
  return(prob_non_zero)
}

#' Model probabilities speeupd
#'
#' @description Computes modeled probabilities of occurrence in the speeupd version
#' of the code
#'
#' @param data
#' @param init_par
#'
#' @return a tibble with modeled probabilities on the logit scale
#' 
modelProbSpeedup <- function(data, 
                             init_par) {
  
  for (par in names(data)) {
    s <- str_c(par, "<- data$", par) 
    eval(expr = parse(text = s))
  }
  # Transformed
  
  log_unif <- -log(12)
  offset_years <- matrix(nrow = 12, ncol = N_u_lower); 
  
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  if (N_groups > 1) {
    do_grouping = 1;
    N_betas_diff = 2;
    N_betas_g2 = 12;
    N_betas_other_g2 = 10;
    N_lambdas = N_units;
    
    if (N_eta_groups > 1) {
      ind_eta_grp2 = 2;
    } else {
      ind_eta_grp2 = 1;
    }
    
  } else  {
    do_grouping = 0;
    N_betas_diff = 0;
    N_betas_g2 = 0;
    N_betas_other_g2 = 0;
    N_lambdas = 0;
    ind_eta_grp2 = 1;
  }
  
  if(do_offsets == 1 && do_grouping == 1){
    N_offset_elems_g2 = 12;
    N_offset_elems = 1;
  } else {
    N_offset_elems_g2 = 1;
  }
  
  # pre-determing vectors of offsetted yearly random effects for each offset
  for (i in 0:11) {
    for (k in 1:N_u_lower) {
      if (do_offsets == 1) {
        offset_years[i+1, k] = (umonths[k] + i) > 12 ? uyears_full[k] + 1 : uyears_full[k];
      } else {
        offset_years[i+1, k] = uyears_full[k];
      }
    }
  }
  
  for (par in names(init_par)) {
    s <- str_c(par, "<- init_par$", par) 
    eval(expr = parse(text = s))
  }
  
  # Transformed parameters
  sigma = 1/sqrt(tau); # overall marginal standard deviation
  ulogitp <- vector(mode = "numeric", length = N_u_lower)
  prob_non_zero <- array(dim = c(N_offset_elems, N_offset_elems_g2, N_obs))  # probabilities of cholera at the upper spatial level (admin 1 or country)
  
  if (!exists("lambdas")) {
    lambdas <- rep(1, N_units)
  }
  
  # Compute probs
  sum_log1p <- 0;
  sum_l <- 0;
  
  for (z1 in 1:N_offset_elems) {
    for (z2 in 1:N_offset_elems_g2) {
      off_uyears = offset_years[z1, ]; # unique years
      off_uyears_g2 = offset_years[z2, ]; # unique years
      
      # compute unique lower level logitps
      for (i in 1:N_u_lower) {
        ulogitp[i] = lambdas[uunits[i]] * (betas_g1[umonths[i]] + etas[off_uyears[i], 1]) + 
          do_grouping * (1-lambdas[uunits[i]]) * (betas_g2[umonths[i]] + etas[off_uyears_g2[i], N_eta_groups]) +
          convolved_re[uunits[i]] * sigma;
      }
      
      # loop over observations and compute the probability of occurrence
      for (start_idx in 1:N_obs) {
        xi = xis[map_2_lev[start_idx]];
        sum_log1p = 0;
        sum_l = 0;
        
        for (k in starts[start_idx]:ends[start_idx]) {
          x = xi + ulogitp[map_2_u[k]];
          sum_l = sum_l + x;
          sum_log1p = sum_log1p + qgam::log1pexp(-1 * x);   
        }
        
        prob_non_zero[z1, z2, start_idx] = log_diff_exp(sum_log1p, -sum_l) + sum_l;
      }
    }
  }  
  
  return(prob_non_zero)
}



#' Perdict probabilities
#'
#' @description Predicts the probabilities of observing cholera for given model choices
#'
#' @param param
#'
#' @return return
#' 
predictProbs <- function(country,
                         run_level,
                         identifier,
                         time_left = NULL,
                         time_right = NULL,
                         model = "all",
                         nsim,
                         nsimbinom = NULL,
                         redo = F,
                         redo_single = redo,
                         redo_data = F,
                         do_par = T,
                         n_cores = parallel::detectCores()/2,
                         only_fake = F,
                         fake_full = F,
                         verbose = F,
                         rm_gadmlev = T) {
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo_data,
                                 verbose = verbose)  
  }
  
  if (model == "best") {
    model <- getBestModel(country = country,
                          run_level = run_level,
                          identifier = identifier,
                          time_left = time_left,
                          time_right = time_right,
                          redo = T)
  }
  
  runChecks(country = country,
            run_level = run_level,
            identifier = identifier,
            model = model)
  
  sample_file <- makeStdResName(country = country, 
                                run_level = run_level,
                                identifier = identifier, 
                                model = model,
                                time_left = time_left,
                                time_right = time_right,
                                suffix = ifelse(only_fake, ifelse(fake_full, "sample_stats_full_fake", "sample_stats_fake"), "sample_stats"),
                                file_type = "rds")
  
  if (rm_gadmlev) {
    sample_file <- str_remove(sample_file, str_c("_", setGadmLev(country, run_level)))
  }
  
  if (checkFile(sample_file) & !redo & !redo_single) {
    res <- readRDS(sample_file) %>% 
      mutate(run_level = as.double(run_level))
    
    return(res)
  } else {
    
    # Export probs one by one and combined
    if (model == "all" | length(model) > 1) {
      
      if (model == "all") {
        models <- getAllModels()
      } else {
        models <- model
      }
      
      prob_sample_stats <- foreach(mit = models,
                                   .combine = dplyr::bind_rows,
                                   .errorhandling = "remove") %do%
        {
          
          done <- makeStdResName(country = country, 
                                 run_level = run_level,
                                 identifier = identifier, 
                                 model = mit,
                                 time_left = time_left,
                                 time_right = time_right,
                                 suffix = "sample_stats",
                                 file_type = "rds") %>% 
            checkFile()
          
          predictProbs(country = country,
                       run_level = run_level,
                       identifier = identifier,
                       model = mit,
                       time_left = time_left,
                       time_right = time_right,
                       nsim = nsim,
                       nsimbinom = nsimbinom,
                       redo = !done,
                       redo_single = !done,
                       do_par = do_par,
                       n_cores = n_cores,
                       only_fake = only_fake) %>% 
            dplyr::mutate(run_level = as.numeric(run_level))
        }
      
    } else {
      
      cat("Computing probs for", country, identifier, run_level, model, "\n")
      
      res_files <- getResFiles(country = country,
                               run_level = run_level,
                               identifier = identifier,
                               model = model,
                               time_left = time_left,
                               time_right = time_right,
                               what = "fit")
      
      prob_sample_stats <- foreach(res_file = res_files,
                                   .combine = bind_rows) %do% 
        {
          
          # Load result
          res <- readRDS(res_file)
          
          run_info <- parseResFilename(res_file)
          
          # Model variant
          variant <- run_info$variant
          model_attributes <- try(attr(res$chol_stanfit, "stan_args")[[1]])
          
          if (inherits(model_attributes, "try-error")) {
            return(NULL)
          }
          
          n_samples <- (model_attributes$iter - model_attributes$warmup) * 4
          
          # Load data mapping 
          data_mapping_file <- str_remove(res_file, "model_fits/results_marcc_2|model_fits/results_marcc") %>% 
            str_replace("fit", "data_mapping")
          
          if (file.exists(data_mapping_file)) {
            data_mapping <- readRDS(data_mapping_file)
            cat("Reading ", data_mapping_file, "\n")
          } else {
            data_mapping <- list(u_unit = seq(res$data$N_units))
          }
          
          # Make fake data to compute all seasonal variables
          fake_data <- makeFakeData(data = res$data, 
                                    u_unit = data_mapping$u_unit)
          
          # Get names of initial parameters
          init_par_names <- names(model_attributes$init) %>% 
            str_subset("ind_", negate = T)
          
          if (length(init_par_names) == 0) {
            init_par_names <- getParNames(res$chol_stanfit)
          }
          
          init_par_names <- setdiff(init_par_names,
                                    c("etas_tilde", "phi", "theta", "tau_beta", "betas_other_g2", "betas_diff", "lambdas_raw", "theta_lambdas_raw"))
          if (str_detect(variant, "mixture")) {
            init_par_names <- c(init_par_names, "betas_g2", "lambdas")
          }
          init_par_names <- c(init_par_names, "etas", "convolved_re")
          
          init_par_names <- init_par_names %>% 
            unique() %>% 
            str_subset("ll|log_lik|lp__|prob|ulogitp|lp", negate = T)
          
          par_samples <- rstan::extract(res$chol_stanfit, pars = init_par_names) 
          
          pars <- map(1:n_samples, 
                      function(i) {
                        map2(par_samples, names(par_samples), function(x, y) {
                          if (length(dim(x)) == 1) {
                            x[i]
                          } else if (length(dim(x)) == 2) {
                            x[i,,drop=T]
                          } else {
                            abind::adrop(x = x[i,,, drop=F], 
                                         drop = c(T, F, F))
                          }
                        })
                      })
          
          if (is.null(res$data$N_eta_groups)) {
            res$data$N_eta_groups <- 1
            res$data$ind_betas_diff <- c(1,1)
            res$data$beta_months <- rep(1, 10)
            
            fake_data$N_eta_groups <- 1
            fake_data$ind_betas_diff <- c(1,1)
            fake_data$beta_months <- rep(1, 10)
          }
          
          mean_pars <- map(par_samples, function(x) {
            if (length(dim(x)) == 1) {
              mean(x) 
            } else {
              apply(x, 2:max(length(dim(x))), mean)
            }
          })
          
          if (!fake_full) {
            # Compute probabilities for fake data (all unit/time combinations at the lowest admin level)
            prob_samples_fake <- computeSingleProb(data = fake_data,
                                                   init_par = mean_pars,
                                                   variant = variant)
            
            # Get offset estimates
            model_offset <- extractOffset(variant, res$chol_stanfit)
            
            if (is.null(model_offset)) {
              prob_samples_fake <- prob_samples_fake[,1,]
            } else {
              if (variant == "mixture_offset") {
                # Compute weighted mean across offset values
                prob_samples_fake <- matrix(model_offset, nrow = 1) %*% prob_samples_fake[1,,]
              } else {
                # Compute weighted mean across offset values
                prob_samples_fake <- matrix(model_offset, nrow = 1) %*% prob_samples_fake[,1,]
              }
              prob_samples_fake <- prob_samples_fake[1,]
            }   
            
            
            prob_samples_fake <- expit(prob_samples_fake)
            
            # Compute statistics
            prob_sample_fake_stats <- tibble(mean = prob_samples_fake,
                                             median = NA,
                                             q025 = NA,
                                             q975 = NA) %>% 
              mutate(obs = row_number())
          } else {
            # Compute probabilities for data 
            prob_samples_fake <- computeProbs(data = fake_data,
                                              init_par = pars[sample(1:n_samples, nsim)],
                                              do_par = do_par,
                                              n_cores = n_cores,
                                              variant = variant)
            # Get offset estimates
            model_offset <- extractOffset(variant, res$chol_stanfit)
            
            if (is.null(model_offset)) {
              prob_samples_fake <- prob_samples_fake[,1,1,]
            } else {
              if (variant == "mixture_offset") {
                # Compute weighted mean across offset values
                prob_samples_fake <- t(apply(prob_samples_fake, 1, function(x) {
                  matrix(model_offset, nrow = 1) %*% x[1,,]
                }))
              } else {
                # Compute weighted mean across offset values
                prob_samples_fake <- t(apply(prob_samples_fake, 1, function(x) {
                  matrix(model_offset, nrow = 1) %*% x[,1,]
                }))
              }
            }   
            
            # Compute probabilities
            prob_samples_fake <- expit(prob_samples_fake)
            
            
            # Compute statistics
            prob_sample_fake_stats <- apply(prob_samples_fake, 2, 
                                            function(x) tibble(mean = mean(x),
                                                               median = median(x),
                                                               q025 = quantile(x, 0.025),
                                                               q975 = quantile(x, 0.975))) %>% 
              bind_rows() %>% 
              mutate(obs = row_number())
            
            # Expand to samples of observations
            fake_obs_samples <- map(1:ncol(prob_samples_fake), ~ prob_samples_fake[, .])
            
            prob_sample_fake_stats <- prob_sample_fake_stats %>% 
              mutate(samples = fake_obs_samples)
          }
          
          
          # Unpack observations
          fake_data_long <- makeData(fake_data, data_mapping)
          
          
          if (!only_fake) {
            
            # Unpack observations
            data <- makeData(data = res$data, 
                             data_mapping = data_mapping)
            
            # Compute probabilities for data 
            prob_samples <- computeProbs(data = res$data,
                                         init_par = pars[sample(1:n_samples, nsim)],
                                         do_par = do_par,
                                         n_cores = n_cores,
                                         variant = variant)
            
            if (is.null(model_offset)) {
              prob_samples <- prob_samples[,1,1,]
            } else {
              if (variant == "mixture_offset") {
                # Compute weighted mean across offset values
                prob_samples <- t(apply(prob_samples, 1, function(x) {
                  matrix(model_offset, nrow = 1) %*% x[1,,]
                }))
              } else {
                # Compute weighted mean across offset values
                prob_samples <- t(apply(prob_samples, 1, function(x) {
                  matrix(model_offset, nrow = 1) %*% x[,1,]
                }))
              }
            }   
            # Compute probabilities
            prob_samples <- expit(prob_samples)
            
            
            # Compute statistics
            prob_sample_stats <- apply(prob_samples, 2, 
                                       function(x) tibble(mean = mean(x),
                                                          median = median(x),
                                                          q025 = quantile(x, 0.025),
                                                          q975 = quantile(x, 0.975))) %>% 
              bind_rows() %>% 
              mutate(obs = row_number())
            
            # Expand to samples of observations
            obs_samples <- map(
              1:ncol(prob_samples),
              function(y) {
                map(1:nsimbinom, ~rbinom(nrow(prob_samples), 
                                         size = res$data$n_replicas[y], 
                                         prob = prob_samples[, y])) %>% 
                  unlist()
              })
            
            
            # Compute statistics
            prob_sample_stats <- prob_sample_stats %>% 
              mutate(samples = obs_samples,
                     obs_mean = map_dbl(obs_samples, ~ mean(.)),
                     obs_q025 = map_dbl(obs_samples, ~ quantile(., .025)),
                     obs_q975 = map_dbl(obs_samples, ~ quantile(., .975)))
            
            out <- data %>% 
              group_by(obs) %>% 
              slice(1) %>% 
              inner_join(dplyr::filter(prob_sample_stats, obs %in% .$obs), by = "obs") %>% 
              mutate(set = "observations") %>% 
              bind_rows(fake_data_long %>% 
                          inner_join(dplyr::filter(prob_sample_fake_stats, obs %in% .$obs), by = "obs") %>% 
                          mutate(set = "fake")) %>% 
              ungroup() %>% 
              bind_cols(run_info)
          } else {
            out <- fake_data_long %>% 
              inner_join(dplyr::filter(prob_sample_fake_stats, obs %in% .$obs), by = "obs") %>% 
              mutate(set = "fake")%>% 
              ungroup() %>% 
              bind_cols(run_info)
          }
          
          cat("Done prob prediction for ", country, identifier, run_level, model, "\n")
          return(out)
        }
    }
    
    # Get admin unit information
    prob_sample_stats <- dplyr::inner_join(prob_sample_stats %>% 
                                             dplyr::mutate(unit = as.numeric(as.character(unit))), 
                                           getAdminUnitMapping(country, run_level),
                                           by = c("unit" = "serial_id"))
    
    saveRDS(prob_sample_stats, file = sample_file)
    return(prob_sample_stats)
  }
}

# Helper functions --------------------------------------------------------

#' expit
#'
#' @description compute the inverse logit transform
#' 
#' @param x value to transform
#'
#' @return a double
#'
expit <- function(x) {
  1/(1+exp(-1*x))
}

#' Log1 minus exp
#'
#' @param x 
#'
#' @return
#'
log1m_exp <- function(x) {
  # https://matt-graham.github.io/mici/docs/utils.html#mici.utils.log1m_exp
  if (x >= 0.0) {
    return(NaN)
  } else if (x > log(2)){
    return(log(-expm1(x)))
  } else {
    return(log1p(-exp(x)))
  }
}

#' log diff exp
#'
#' @param val1 
#' @param val2 
#'
#' @return
#' @export
#'
#' @examples
log_diff_exp <- function(val1, val2) {
  # https://matt-graham.github.io/mici/docs/utils.html#mici.utils.log1m_exp
  if (val1 == -Inf & val2 == -Inf) {
    return(-Inf)
  } else if (val1 < val2){
    return(NaN)
  } else if (val1 == val2) {
    return(-Inf)
  } else {
    return(val1 + log1m_exp(val2 - val1))
  }
}

#' title Get data source
#' 
#' @description Gets the uniue OC ids that compose the observations of a given admin unit country
#' 
#' @param param
#' 
#' @return a dataframe with the sources and observed case
#' 
getDataSource <- function(country, 
                          gadm_long_id = NULL,
                          gadm_id = NULL, 
                          gadm_lev = NULL,
                          username,
                          month_left = NULL,
                          month_right = NULL) {
  
  conn <- DBI::dbConnect(RPostgres::Postgres(), 
                         user = username,
                         dbname = "taxdat")
  
  if (is.null(gadm_lev)) {
    if (is.null(gadm_id_long))
      stop("please provide either gadm_long_id, or {gadm_id, gadm_lev}")
    
    gadm_lev <- str_count(gadm_id_long, "\\.")
  }
  
  if (is.null(gadm_long_id)) {
    if (is.null(gadm_id))
      stop("please provide either gadm_long_id, or {gadm_id, gadm_lev}")
    
    # Get long id
    gadm_long_id <- DBI::dbGetQuery(
      conn, 
      glue::glue_sql("SELECT gid_{gadm_lev} as id
      FROM admin.gadm_lev{gadm_lev}
                     WHERE id = {gadm_id}", .con = conn))$id[1]
  }
  
  query <- "SELECT * FROM cholera.monthly_data WHERE gadm_id = {gadm_id} AND gadm_lev = {gadm_lev}"
  
  if (!is.null(month_left)) {
    query <- str_c(query, " AND month_left >= {month_left}")
  }
  
  if (!is.null(month_left)) {
    query <- str_c(query, " AND month_right <= {month_right}")
  }
  
  monthly_data <- DBI::dbGetQuery(conn, glue::glue_sql(query, .con = conn))
  
  query2 <- "
  SELECT DISTINCT d.location_period_id, b.in_parts, country, time_left as tl, time_right as tr,
observation_collection_id as oc_id, suspected_cases,
FROM admin.gadm_lev{gadm_lev} a
JOIN
admin.gadm_lev{gadm_lev}_dict b
ON a.id = b.gadm_id
JOIN
admin.location_periods_dict c
ON b.lp_id = c.location_period_id
JOIN
cholera.data d
ON c.observation_id = d.id
WHERE a.gid_{gadm_lev} = {gadm_long_id}
  "
  
  if (!is.null(month_left)) {
    query2 <- str_c(query2, " AND d.time_left >= {month_left}")
  }
  
  time_right <- month_right
  lubridate::month(time_right) <- lubridate::month(time_right) + 1
  time_right <- time_right - 1
  
  if (!is.null(month_left)) {
    query2 <- str_c(query2, " AND d.time_right <= {time_right}")
  }
  
  query2 <- str_c(query2, " ORDER BY oc_id, location_period_id, in_parts, tl, tr")
  full_data <- DBI::dbGetQuery(conn, glue::glue_sql(query2, .con = conn))
  
  
  return(sources)
}


#' Get model version (deprecated)
#'
#' @param data 
#' @param verbose 
#'
#' @return
#'
getModelVersion <- function(data,
                            verbose = F) {
  
  model_version <- ifelse(!is.null(data$N_u_lower), "speedup", "old")
  
  if (verbose)
    cat("Model version:", model_version, "\n")
  
  return(model_version)
}

#' Make data
#'
#' @description make a stan input dataset to run probability computatoins
#' 
#' @param data 
#' @param data_mapping 
#'
#' @return
#'
makeData <- function(data, 
                     data_mapping) {
  
  model_version <- getModelVersion(data)
  
  if (model_version == "old") {
    # Unpack observations
    data <- tibble(
      obs = data$map_2_obs,
      month = data$months,
      year = data$years,
      unit = data$map_2_units,
      lev = data$map_2_lev
    ) %>% 
      inner_join(
        tibble(
          obs = 1:data$N_obs,
          n_pos = data$y,
          n_obs = data$n_replicas
        ), by = "obs"
      ) %>% 
      mutate(unit = factor(unit),
             lev = factor(lev)) %>% 
      group_by(obs) %>% 
      mutate(n_units = length(unique(unit)),
             n_months = length(unique(month)))
  } else {
    # Unpack observations
    data <- tibble(
      obs = map_dbl(1:data$N_obs, function(x) {
        rep(x, data$ends[x] - data$starts[x] + 1)
      }),
      month = map_dbl(1:data$N_obs, function(x) {
        data$umonths[data$starts[x]:data$ends[x]]
      }),
      year = map_dbl(1:data$N_obs, function(x) {
        data$uyears[data$starts[x]:data$ends[x]]
      }),
      unit = map_dbl(1:data$N_obs, function(x) {
        data$uunits[data$starts[x]:data$ends[x]]
      }),
      lev = data$map_2_lev
    ) %>% 
      inner_join(
        tibble(
          obs = 1:data$N_obs,
          n_pos = data$y,
          n_obs = data$n_replicas
        ), by = "obs"
      ) %>% 
      mutate(unit = factor(unit),
             lev = factor(lev)) %>% 
      group_by(obs) %>% 
      mutate(n_units = length(unique(unit)),
             n_months = length(unique(month)))
    
  }
  
  if (!is.null(data_mapping$u_years)) {
    data$datenum <- data_mapping$u_years[data$year] + (data$month-1)/12
    data$year <- factor(data_mapping$u_years[data$year])
  } else {
    data$datenum <- data$year + (data$month-1)/12
    data$year <- factor(data$year)
  }
  
  return(data)
}

#' Make fake data
#'
#' @description Make fake stan input dataset for probability computations
#' 
#' @param data 
#' @param u_unit 
#' @param verbose 
#'
#' @return
#'
makeFakeData <- function(data,
                         u_unit, 
                         verbose = F) {
  
  # Get model version
  model_version <- getModelVersion(data, verbose = verbose)
  
  # Fake data for time series
  fake_data <- data
  
  if (model_version == "old") {
    
    fake_data_combs <- expand.grid(
      unit = seq_along(u_unit),
      years = unique(data$year),
      months = 1:12
    ) %>% 
      as_tibble() %>% 
      mutate(obs = row_number())
    
    n_fake_obs <- nrow(fake_data_combs)
    fake_data$years <- fake_data_combs$years
    fake_data$months <- fake_data_combs$months
    fake_data$map_2_lev <- rep(1, n_fake_obs)
    fake_data$map_2_obs <- 1:n_fake_obs
    fake_data$map_2_units <- fake_data_combs$unit
    fake_data$K <- n_fake_obs
    fake_data$N_obs <- n_fake_obs
    fake_data$N_obs_lower <- n_fake_obs
    fake_data$N_obs_upper <- 0
    fake_data$N_components <- rep(1, n_fake_obs)
    fake_data$n_replicas <- rep(1, n_fake_obs)
    fake_data$N_map_lower <- n_fake_obs
    fake_data$N_map_upper <- 0
    fake_data$ind_lower <- 1:n_fake_obs
    fake_data$obs_ind_lower <- 1:n_fake_obs
    fake_data$ind_upper <- array(dim = c(0))
    fake_data$obs_ind_upper <- array(dim = c(0))
    fake_data$map_2_upper_obs <- rep(0, n_fake_obs)
    fake_data$y <- rep(0, n_fake_obs)
  } else {
    
    fake_data_combs <- tibble(
      unit = data$uunits,
      years = data$uyears,
      months = data$umonths
    ) %>% 
      mutate(obs = row_number())
    
    n_fake_obs <- data$N_u_lower
    fake_data$map_2_lev <- rep(1, n_fake_obs)
    fake_data$map_2_u <- 1:n_fake_obs
    fake_data$starts <- 1:n_fake_obs
    fake_data$ends <- 1:n_fake_obs
    fake_data$K <- n_fake_obs
    fake_data$N_obs <- n_fake_obs
    fake_data$n_replicas <- rep(1, n_fake_obs)
    fake_data$y <- rep(0, n_fake_obs)
  }
  return(fake_data)
}

#' Width to interval
#'
#' @param wd 
#'
#' @return
#'
widthToInt <- function(wd) {
  halfwd <- (1 - wd)/2
  c("low" = halfwd, "high" = 1 - halfwd)
}

