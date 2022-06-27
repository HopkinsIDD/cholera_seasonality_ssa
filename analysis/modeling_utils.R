# Filenames ---------------------------------------------------------------

makeRunDataFile <- function(opt, config) {
  rundata_dir <- str_glue("{opt$cholera_dir}/generated_data/run_data")
  
  if (!dir.exists(rundata_dir)) {
    dir.create(rundata_dir)
  }
  
  glue::glue("{rundata_dir}/run_data_{config$country_iso3}_{makeCholeraThresh(config)}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}.rds")
}

makeAdjacencyFile <- function(opt, config, run_level) {
  adj_dir <- str_glue("{opt$cholera_dir}/generated_data/adjacency")
  
  if (!dir.exists(adj_dir)) {
    dir.create(adj_dir)
  }
  
  glue::glue("{adj_dir}/adjacency_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}.rds")
}

makePriorRhoFile <- function(opt, config, run_level) {
  rho_dir <- str_glue("{opt$cholera_dir}/generated_data/rho_priors")
  
  if (!dir.exists(rho_dir)) {
    dir.create(rho_dir)
  }
  
  glue::glue("{rho_dir}/prior_rho_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}.rds")
}

makeMappingFile <- function(opt, config, run_level) {
  mapping_dir <- str_glue("{opt$cholera_dir}/generated_data/data_mapping")
  
  if (!dir.exists(mapping_dir)) {
    dir.create(mapping_dir)
  }
  
  glue::glue("{mapping_dir}/data_mapping_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}_dmy{config$drop_multi_yr}.rds")
}

makeResultsFile <- function(opt, cntry_iso3, run_level, drop_multi_yr)  {
  res_dir <- str_glue("{opt$cholera_dir}/generated_data/model_fits")
  
  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }
  
  glue::glue("{res_dir}/fit_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}_dmy{config$drop_multi_yr}.rds")
}

makePriorFile <- function(opt, cntry_iso3, run_level, drop_multi_yr)  {
  res_dir <- str_glue("{opt$cholera_dir}/generated_data/priors")
  
  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }
  
  glue::glue("{res_dir}/prior_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}_dmy{config$drop_multi_yr}.rds")
}

makeOutputFile <- function(res_file) {
  
  if (!dir.exists("generated_data/outputs")) {
    dir.create("generated_data/outputs")
  }
  
  str_replace(res_file, "model_fits", "outputs") %>% 
    str_replace("fit", "output")
}

makePriorOutputFile <- function(res_file) {
  
  if (!dir.exists("generated_data/outputs_prior")) {
    dir.create("generated_data/outputs_prior")
  }
  
  str_replace(res_file, "priors", "outputs_prior") %>% 
    str_replace("prior_", "output_prior_")
}

makeOutputSimFile <- function(res_file) {
  
  if (!dir.exists("generated_data/outputs_sim")) {
    dir.create("generated_data/outputs_sim")
  }
  
  str_replace(res_file, "model_fits", "outputs_sim") %>% 
    str_replace("fit", "output_sim")
}

makeCmdstanrFile <- function(res_file, what = "model_fits") {
  
  if (!dir.exists("generated_data/cmdstanr")) {
    dir.create("generated_data/cmdstanr")
  }
  str_replace(res_file, what, "cmdstanr")
}


makeCmdstanrSimFile <- function(res_file, what = "model_fits") {
  
  if (!dir.exists("generated_data/cmdstanr_sim")) {
    dir.create("generated_data/cmdstanr_sim")
  }
  
  str_replace(res_file, what, "cmdstanr_sim") %>% 
    str_replace("fit", "sim")
}

makeGAMParmFile <- function(opt, cntry_iso3, run_level)  {
  gam_dir <- str_glue("{opt$cholera_dir}/generated_data/gam_par")
  
  if (!dir.exists(gam_dir)) {
    dir.create(gam_dir)
  }
  
  glue::glue("{gam_dir}/{cntry_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}_runlev{run_level}_init_pars.png")
}

makeJsonFile <- function(opt, config) {
  json_dir <- str_glue("{opt$cholera_dir}/generated_data/json_data")
  
  if (!dir.exists(json_dir)) {
    dir.create(json_dir)
  }
  
  glue::glue("{json_dir}/stan_data_{config$country_iso3}_{makeCholeraThresh(config)}_{config$variant}_{makePeriod(config)}_k{config$keep}_l{config$gadm_levels}.json")
}


makeIntermediateFileBase <- function(res_file, split = "fits/") {
  int_dir <- "interm/cmdstanr_int/"
  
  if (!dir.exists(int_dir)) {
    dir.create(int_dir, recursive = T)
  }
  
  res_file %>% 
    str_remove("\\.rds") %>% 
    str_split(split) %>% 
    .[[1]] %>% 
    .[2] %>% 
    str_c(int_dir, .)
}

# Data processing ---------------------------------------------------------

#' @title Compute adjacency
#' @description computes the adjacency matrix for the conutry at a given GADM level
#'
#' @param cntry.sf An sf object of the country GADM2 unists
#' @param run_level GADM level at which to run the model
#'
#' @return a list with the modified sf object and the adjacency information
computeAdjacency <- function(country,
                             run_level) {
  
  cntry.sf <- getGadmSf(country = country, gadm_lev = run_level)
  
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


#' @title Compute ICAR priors
#' @description Computes ICAR complexity-adjusted prior for the spatial variance contribution
#'
#' @param N number of nodes in adjacency matrix
#' @param node1 origin nodes of sparse adjacency matrix
#' @param node2 destination nodes of sparse adjacency matrix
#' @param U_rho parameter controlling the prior of the spatial variance contribution
#' @param alpha_rho parameter controlling the prior of the spatial variance contribution
#' 
#' @details Complexity-adjusted hyper priors follow the specifications in Riebler et al. 2016
#'  
#' @return a list with hyper-priors
computeICARPriors <- function(N,
                              node1, 
                              node2,
                              U_rho = .5,
                              alpha_rho = 2/3
){
  
  # Use this for parallel computing
  # suppressMessages(INLA:::inla.dynload.workaround())
  
  # Build the adjacency matrix using INLA library functions
  adj.matrix <- sparseMatrix(i = node1, j = node2, x = 1, symmetric = TRUE)
  # The ICAR precision matrix (note! This is singular)
  Q <- Diagonal(N, rowSums(adj.matrix)) - adj.matrix
  # Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert <- Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  Q_scaled <- suppressMessages(INLA::inla.scale.model(Q_pert, constr = list(A = matrix(1, 1, N), e = 0)))
  # Compute the diagonal elements of the covariance matrix subject to the
  # constraint that the entries of the ICAR sum to zero.
  # See the inla.qinv function help for further details.
  Q_inv <- suppressMessages(INLA::inla.qinv(Q_pert, constr = list(A = matrix(1, 1, N), e = 0)))
  Q_inv_full <- suppressMessages(INLA:::inla.ginv(Q_scaled))
  
  # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <- exp(mean(log(diag(Q_inv))))
  
  # Prior for the fraction of variance that comes from connectedness
  prior_rho_tab <- suppressMessages(INLA:::inla.pc.bym.phi(Q = Q_scaled, rankdef = 1, u = U_rho, alpha = alpha_rho, return.as.table = T))
  prior_rho_num <- str_split(prior_rho_tab, " ")[[1]][-1] %>% as.numeric()
  n_pc_rho <- length(prior_rho_num)/2
  log_prior <- prior_rho_num[(n_pc_rho + 1):(2 * n_pc_rho)] 
  phi_intern <- prior_rho_num[1:n_pc_rho]
  prior_rho_df <- tibble(rho = phi_intern, ll = log_prior + phi_intern + 2*log(1+exp(-phi_intern)))
  
  res <- list(
    scaling_factor = scaling_factor,
    prior_rho_df = prior_rho_df,
    Q_inv_full = Q_inv_full
  )
  return(res)
}

#' Get run data
#' 
#' @description Gets the data to run the seasonality model
#'
#' @param path_to_monthly_data path to the monthly data file. This is pre-computed in XX.R
#' @param gadm_dicts list of dictionaries between GADM ids and shapefile ids
#' @param gadm_levels which gadm levels in data to keep
#' @param frac_days_thresh fraction of the days in a month for which data is removed if cholera == 0
#' @param keep_national_yearly whether to keep yearly national level data
#' @param start_date start date to dplyr::filter data
#' @param end_date
#'
#' @return a dataframe with monthly cholera occurences
#' 
getRunData <- function(cholera_directory = getCholeraDirectory(),
                       path_to_monthly_data = getMonthlyDataPath(),
                       country_iso3,
                       gadm_levels = c(0, 1, 2),
                       frac_days_thresh = .8,
                       keep_national_yearly = T,
                       start_date = NULL,
                       end_date = NULL,
                       scraps = F,
                       incid_tresh = 1e-3,
                       case_thresh = 10,
                       thresh = "none",
                       verbose = T) { 
  
  # ---
  # Input parsing
  # ---
  if (is.character(gadm_levels)) {
    gadm_levels <- map_dbl(1:nchar(gadm_levels), 
                           ~str_sub(gadm_levels, . ,.) %>% 
                             as.integer())
  }
  
  if (!is.null(start_date)) {
    if (start_date == "all") {
      start_date <- NULL
    }
  }
  
  if (!is.null(end_date)) {
    if (end_date == "all") {
      end_date <- NULL
    }
  }
  
  
  # ---
  # A. Get GADM data 
  # ---
  
  gadm_data <- getGadmData(cholera_directory = cholera_directory) %>% 
    dplyr::filter(country %in% country_iso3) %>% 
    dplyr::select(-popsum_2010_2016) %>% 
    # Pivot to have year as variable
    tidyr::pivot_longer(cols = contains("popsum"),
                        names_to = "year",
                        values_to = "popsum") %>%
    dplyr::mutate(year = as.integer(stringr::str_extract(year, "[0-9]{4}")))
  
  # List of GADM units in country
  gadm_list <- getGadmList(country = country_iso3)
  
  # Load monthly data of country
  monthly_data <- vroom::vroom(path_to_monthly_data, 
                               col_types = cols()) %>%
    # Filter for known GADM units
    dplyr::inner_join(gadm_list) %>%
    # Tidy for processing
    dplyr::mutate(
      month = lubridate::month(month_left),
      year = year_left,
      n_days_month = lubridate::days_in_month(month_left),
      frac_days = n_days / n_days_month
    ) %>% 
    # Remove data that is 0 and does not cover the whole month
    dplyr::filter(!(frac_days < frac_days_thresh & pmax(suspected_cases, confirmed_cases) == 0)) %>% 
    dplyr::left_join(gadm_data) %>% 
    # Fill missing population data for years before 2000
    dplyr::group_by(gid) %>% 
    dplyr::group_modify(function(x, y) {
      replacement <- gadm_data %>% 
        dplyr::filter(year == 2000, gid == y$gid[1]) %>% 
        dplyr::select(cholera_cases_2010_2016, cholera_rates_2010_2016, popsum) %>% 
        dplyr::slice(1) %>% 
        as.list()
      tidyr::replace_na(x, replacement)
    }) %>% 
    dplyr::rowwise() %>% 
    # Months in observations 
    dplyr::mutate(
      months = list(seq.Date(min(month_left), max(month_right), by = "1 months")),
      n_months = length(seq.Date(min(month_left), max(month_right), by = "1 months"))) %>% 
    dplyr::ungroup()
  
  # Filter by modelling start and end dates
  if (!is.null(start_date)) {
    monthly_data <- monthly_data %>%
      dplyr::filter(monthly_data >= as.Date(start_date))
  }
  
  if (!is.null(end_date)) {
    monthly_data <- monthly_data %>%
      dplyr::filter(month_right <= as.Date(end_date))
  }
  
  
  # ---
  # B. Data processing 
  # ---
  
  # This step involves filtering the data and processing as follows:
  
  # Setup observations vector: sum of months with or without cholera per admin unit and per year
  # each observation collection x period id combination constitutes a separate observation
  run_data <- monthly_data %>%
    # Remove multi-year observations
    dplyr::filter(!str_detect(period_type, "multi_year")) 
  
  # Compute outcome variables for each unique month x location across OCs
  run_data <- run_data %>%
    dplyr::group_by(month_left, month_right, year, gid, gadm_lev, in_lower, in_upper) %>%
    dplyr::summarise(
      n_obs = n(),
      maxcases = max(pmax(suspected_cases, confirmed_cases, na.rm = T)),
      maxincid = max(maxcases/popsum),
      cholera = computeOutcomeVar(suspected_cases,
                                  confirmed_cases,
                                  popsum,
                                  n_months,
                                  cholera_rates_2010_2016,
                                  thresh = thresh, 
                                  incid_tresh = incid_tresh,
                                  case_thresh = case_thresh)
    ) %>% 
    dplyr::arrange(month_left, desc(gadm_lev), gid, year) %>% 
    dplyr::ungroup() %>% 
    # Set a unique observation id
    dplyr::mutate(obs_id = row_number())
  
  
  # ---
  # C. Data filtering and cleaning
  # ---
  
  # Filter out years for which we only have yearly data
  multimonth_years <- run_data %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      months = list(seq.Date(min(month_left), max(month_right), by = "1 months")),
      n_months = length(seq.Date(min(month_left), max(month_right), by = "1 months"))) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(year) %>%
    dplyr::summarise(n_month = min(purrr::map_dbl(months, ~ length(.)))) %>%
    dplyr::filter(n_month < 12) %>%
    .[["year"]]
  
  if (scraps) {
    res <- list(multimonth = dplyr::filter(run_data, !(year %in% multimonth_years)))
  }
  
  run_data <- dplyr::filter(run_data, year %in% multimonth_years) %>% 
    dplyr::filter(gadm_lev %in% gadm_levels)
  
  # If flag require national-level yearly data
  if (!keep_national_yearly) {
    if (scraps) {
      res <- append(res, list(yearly = dplyr::filter(run_data, !(gadm_lev == 0 & purrr::map_dbl(months, ~ length(.)) >= 12))))
    }
    run_data <- dplyr::filter(run_data, !(gadm_lev == 0 & purrr::map_dbl(months, ~ length(.)) >= 12))
  }
  
  if (scraps) {
    return(append(res, list(data = run_data)))
  }
  
  # Remove data that is not reliable defined as:
  # 0 counts for location_periods that are smaller than the gadm unit
  # or >0 counts for location_periods that are larger than the gadm units
  dropped_df <- run_data %>%
    dplyr::filter(!((in_upper & in_lower) | (!in_lower & cholera > 0) | (!in_upper & cholera == 0)))
  
  run_data <- run_data %>% 
    dplyr::filter(!(obs_id %in% dropped_df$obs_id))
  
  if (verbose) {
    cat("-- Dropped", nrow(dropped_df), "observations:", "\n")
    
    dropped_df %>% 
      dplyr::group_by(gadm_lev) %>% 
      dplyr::summarise(n_units = length(unique(gid)),
                       n_obs = n(),
                       n_pos = sum(cholera>0)) %>% 
      print()
    
    cat("-- Running with:\n")
    run_data %>% 
      dplyr::group_by(gadm_lev) %>% 
      dplyr::summarise(n_units = length(unique(gid)),
                       n_obs = n(),
                       n_pos = sum(cholera>0)) %>% 
      print()
  }
  
  
  # ---
  # D. Consolidate and return
  # ---
  
  # Regroup information from different in_upper/in_lower combinations
  run_data <- run_data %>%
    dplyr::group_by(month_left, month_right, year, gid, gadm_lev) %>%
    dplyr::summarise(
      n_obs = sum(n_obs),
      cholera = sum(cholera),
      months = list(seq.Date(min(month_left), max(month_right), by = "1 months"))
    ) %>% 
    dplyr::ungroup()
  
  return(run_data)
}



#' @title Instantiate data for stan
#' @description Instanticates the data for stan
#'
#' @param adjacency
#' @param data_mapping
#' @return a list with the data
instantiateStanData <- function(adjacency, 
                                data_mapping,
                                prior_rho,
                                U_tau,
                                alpha_tau) {
  
  N_units <- nrow(adjacency$sf_object)
  N_edges <- adjacency$adj_list$N_edges
  node1 <- adjacency$adj_list$node1
  node2 <- adjacency$adj_list$node2
  
  # Instantiate data
  data <- list(
    K = data_mapping$K,
    N_units = N_units,
    N_edges = N_edges,
    N_u_lower = data_mapping$N_u_lower,
    node1 = node1,
    node2 = node2,
    N_lev = data_mapping$N_lev,
    country_level = as.integer(0 %in% data_mapping$levels),
    N_obs = data_mapping$N_obs,
    N_years = length(data_mapping$u_years),
    N_years_obs = length(unique(data_mapping$uyears)),
    map_2_lev = data_mapping$map_2_lev,
    map_2_year = data_mapping$map_2_year,
    map_2_u = data_mapping$map_2_u,
    uunits = data_mapping$uunits,
    umonths = data_mapping$umonths,
    uyears = data_mapping$uyears,
    uyears_full = data_mapping$uyears_full,
    u_years_ind = unique(data_mapping$uyears_full),
    starts = data_mapping$starts,
    ends = data_mapping$ends,
    y = data_mapping$y,
    n_replicas = data_mapping$n_replicas,
    scaling_factor = prior_rho$scaling_factor,
    tabulated_rho = prior_rho$prior_rho_df$rho, 
    tabulated_pc_prior_rho = prior_rho$prior_rho_df$ll,
    alpha = alpha_tau,
    U = U_tau,
    N_pc_rho = nrow(prior_rho$prior_rho_df))
  
  return(data)
}

#' @title Observations to administrative unit id mapping
#' @description Computes the mapping between observations and administrative unit ids for Stan
#'
#' @param run_data data on which to run the model
#' @param gadm_dicts dictionaries of GADM units
#' @param run_level GADM level at which to run the model
#'
#' @return return
obsToMappingFlat <- function(country,
                             run_data,
                             run_level,
                             variant,
                             version = "base",
                             drop_multi_yr = F) {
  
  gadm_list <- getGadmList(country)
  gadm_mapping <- getGadmMapping(country)
  
  u_gadm1 <- gadm_list$gid[gadm_list$gadm_lev == 1] %>% sort() # unique gadm level 1 ids
  n_gadm1 <- length(u_gadm1) # number of gadm units
  n_subunits_all <- gadm_mapping %>% count(gid_1) %>% arrange(gid_1) %>% pull(n)
  u_gadm2 <- gadm_list$gid[gadm_list$gadm_lev == 2] %>% sort()
  n_gadm2 <- length(u_gadm2)
  u_levels <- unique(run_data$gadm_lev) %>% sort(decreasing = T)
  
  # Unique years
  u_years <- map(run_data$months, ~year(.)) %>%
    unlist() %>%
    unique() %>% 
    sort()
  
  if ("offset" %in% variant) {
    years_padded <- padYears(u_years)  
    ref_years <- years_padded
    n_years <- length(years_padded)
    cat("\n\nRUNNING WITH OFFSET \n")
  } else {
    n_years <- length(u_years)
    ref_years <- u_years
  }
  
  # Initialize maps
  map_2_units <- c()
  map_2_lev <- c()
  map_2_obs <- c()
  months <- c()
  years <- c()
  
  if (drop_multi_yr) {
    dropped_obs <- run_data %>% 
      dplyr::filter(lubridate::year(month_left) != lubridate::year(month_right))
    
    cat("-- Dropping ", nrow(dropped_obs), 
        str_c("(", formatC(nrow(dropped_obs)/nrow(run_data)*100, digits = 2), "%)"),
        "observations spanning multiple years", "\n")
    
    run_data <- run_data %>% 
      dplyr::filter(lubridate::year(month_left) == lubridate::year(month_right))
  }
  
  if (run_level == 2) {
    # RUN LEVEL : 2
    # Create mappings from the admin 2 to admin 1 level
    # Reorder for efficiency
    run_data <- run_data %>% 
      dplyr::group_by(gadm_lev, gid, year) %>% 
      dplyr::arrange(gadm_lev, gid, year, month_left, month_right, .by_group = T)
    
    for (i in 1:nrow(run_data)) {
      n_months <- length(run_data$months[[i]])
      lev <- run_data$gadm_lev[i]
      
      # Get subunits if level 1
      if (lev == 1) {
        lev2_units <- which(u_gadm2 %in% gadm_mapping$gid_2[gadm_mapping$gid_1 == run_data$gid[i]])
      } else if (lev == 2) {
        lev2_units <- which(u_gadm2 == run_data$gid[i])
      } else if (lev == 0){
        lev2_units <- 0  # case for country-level data which is handled in the stan code
      }
      n_units <- length(lev2_units)
      n_elems <- n_months * n_units
      
      # Concatenate maps
      map_2_units <- c(map_2_units, rep(lev2_units, each = n_months))
      map_2_lev <- c(map_2_lev, rep(which(u_levels == lev), n_months * n_units))
      map_2_obs <- c(map_2_obs, rep(i, n_months * n_units))
      months <- c(months, rep(lubridate::month(run_data$months[[i]]), times = n_units))
      years <- c(years, rep(lubridate::year(run_data$months[[i]]), times = n_units))
    }
  }  else {
    # RUN LEVEL : 1
    # Modify data to account for admin2 level observations (if any) at the 
    # admin 1 level if cholera was observed
    if (sum(run_data$gadm_lev == 2) > 0) {
      run_data_mod <- run_data %>%
        # Get cholera occurrences from admin 2 data
        dplyr::filter(gadm_lev == 2, cholera == 1) %>% 
        dplyr::mutate(gid = purrr::map_chr(gid, ~gadm_mapping$gid_1[gadm_mapping$gid_2 == .]),
                      gadm_lev = 1) %>% 
        distinct() %>% 
        # Get absence of cholera from admin 2 data
        rbind(
          run_data %>% 
            dplyr::filter(gadm_lev == 2, cholera == 0) %>% 
            dplyr::mutate(gid_1 = purrr::map_chr(gid, ~gadm_mapping$gid_1[gadm_mapping$gid_2 == .])) %>% 
            dplyr::group_by(month_left, month_right, gid_1) %>% 
            dplyr::summarise(gid_2 = list(gid)) %>% 
            dplyr::mutate(all_gadm2 = map2_lgl(
              gid_1, 
              gid_2, 
              function(x, y){
                id2 <- gadm_mapping$gid_2[gadm_mapping$gid_1 == x]
                # Check if all gadm2 level ids are have data for given gadm1 level id
                !any(!purrr::map_lgl(id2, ~. %in% y))
              } 
            )) %>% 
            ungroup() %>% 
            dplyr::filter(all_gadm2) %>% 
            select(month_left, month_right, gid_1) %>% 
            dplyr::inner_join(run_data %>% 
                                dplyr::filter(gadm_lev == 2, cholera == 0) %>% 
                                dplyr::mutate(gid = purrr::map_chr(gid, ~gadm_mapping$gid_1[gadm_mapping$gid_2 == .]),
                                              gadm_lev = 1) %>% 
                                distinct(),
                              by = c("month_left", "month_right", "gid_1" = "gid")) %>% 
            rename(gid = gid_1)
          
        ) %>% 
        rbind(
          run_data %>% 
            dplyr::filter(gadm_lev != 2)) %>% 
        dplyr::group_by(month_left, month_right, year, gadm_lev, gid) %>% 
        dplyr::summarise(n_obs = sum(n_obs),
                         cholera = sum(cholera),
                         months = list(as.Date(unique(unlist(months)), origin = "1970-01-01"))) %>% 
        ungroup() 
    } else {
      run_data_mod <- run_data
    }
    
    # set level to 2 for the stan codes to run
    run_data_mod <- run_data_mod %>% 
      dplyr::mutate(gadm_lev = case_when(gadm_lev == 1 ~ 2, T ~ 0))
    
    run_data <- run_data_mod
    
    # Reorder for efficiency
    run_data <- run_data %>% 
      dplyr::group_by(gadm_lev, gid, year) %>% 
      dplyr::arrange(gadm_lev, gid, year, month_left, month_right, .by_group = T)
    
    u_levels <- unique(run_data$gadm_lev) %>% sort(decreasing = T)
    
    for (i in 1:nrow(run_data)) {
      n_months <- length(run_data$months[[i]])
      lev <- run_data$gadm_lev[i]
      
      # Get subunits if level 1
      if (lev == 2) {
        lev2_units <- which(u_gadm1 == run_data$gid[i])
      } else if (lev == 0){
        lev2_units <- 0  # case for country-level data which is handled in the stan code
      }
      n_units <- length(lev2_units)
      n_elems <- n_months * n_units
      
      # Concatenate maps
      map_2_units <- c(map_2_units, rep(lev2_units, each = n_months))
      map_2_lev <- c(map_2_lev, rep(which(u_levels == lev), n_months * n_units))
      map_2_obs <- c(map_2_obs, rep(i, n_months * n_units))
      months <- c(months, rep(lubridate::month(run_data$months[[i]]), times = n_units))
      years <- c(years, rep(lubridate::year(run_data$months[[i]]), times = n_units))
    }
  }
  
  # Ungroup
  run_data <- run_data %>% dplyr::ungroup()
  
  # set years to indices
  u_years <- unique(years) %>% sort()
  years_ind <- purrr::map_dbl(years, ~which(ref_years == .))
  u_years_ind <- unique(years_ind) %>% sort()
  
  # Compute number of components that contribute to each observation (months x admin units)
  N_components <- purrr::map_dbl(1:nrow(run_data), 
                                 ~ ifelse(map_2_units[map_2_obs == .][1] == 0, 
                                          n_units, 
                                          length(map_2_units[map_2_obs == .])) * length(months[map_2_obs == .]))
  
  # First pass at output
  res <- list(
    run_data = run_data,
    u_years = ref_years,
    u_gadm2= u_gadm2, 
    u_gadm1 = u_gadm1,
    levels = unique(run_data$gadm_lev) %>% sort(decreasing = T),
    N_lev = length(unique(run_data$gadm_lev)),
    N_obs = nrow(run_data),
    N_components = N_components,
    K = length(map_2_units),
    y = run_data$cholera,
    n_replicas = run_data$n_obs,
    map_2_units = map_2_units,
    map_2_lev = map_2_lev,
    map_2_lev_mapping = map_2_lev,
    map_2_obs = map_2_obs,
    months = months, 
    years = years_ind
  )
  
  if (run_level == 2) {
    res$u_unit <- u_gadm2
  }  else {
    res$u_unit <- u_gadm1
  }
  
  if (version == "speedup") {
    
    # New unique unit x year x month combinations
    uunit <- 1:length(res$u_unit)
    uyears <- 1:length(unique(u_years_ind))
    umonths <- 1:12
    
    # Expand to unique combinations
    ucomb <- expand.grid(
      month = umonths,
      unit = uunit,
      year = uyears
    ) %>% 
      dplyr::mutate(uind = dplyr::row_number(),
                    uyears_full = u_years_ind[year]) %>% 
      tibble::as_tibble()
    
    print(nrow(ucomb))
    
    # Reorder observations to be in the order: country-admin1-admin2 and 
    # yearly-sub-yearly-monthly
    orig_mapping <- 1:length(map_2_units)
    N_obs <- res$N_obs
    
    # Determine the starts and ends of each observation
    starts <- vector(mode = "numeric", length = N_obs)
    ends <- vector(mode = "numeric", length = N_obs)
    map_2_u <- c()
    
    for (i in 1:N_obs) {
      # Get unique combinations
      mapping_ind <- orig_mapping[map_2_obs == i]
      subdat <- tibble::tibble(
        year = years_ind[mapping_ind],
        month =  months[mapping_ind],
        unit = map_2_units[mapping_ind]
      ) %>% 
        dplyr::distinct()
      
      if (subdat$unit[1] == 0) {
        subdat <- purrr::map_df(uunit, function(x) subdat %>% dplyr::mutate(unit = x))
      }
      
      subucomb <- suppressMessages(ucomb %>% 
                                     dplyr::inner_join(subdat,
                                                       by = c("unit" = "unit",
                                                              "month" = "month",
                                                              "uyears_full" = "year")))
      
      starts[i] <- length(map_2_u) + 1
      ends[i] <- length(map_2_u) + length(subucomb$uind) 
      map_2_u <- c(map_2_u, subucomb$uind)
    }
    
    # Observation-level map 2 lev
    map_2_lev_obs <- purrr::map_dbl(run_data$gadm_lev, ~ which(u_levels == .))
    map_2_year_obs <- purrr::map_dbl(run_data$year, ~ which(u_years == .))
    
    res$K = length(map_2_u)
    res$map_2_lev_mapping = res$map_2_lev 
    res$map_2_lev = map_2_lev_obs
    
    # Concatenate results
    res <- c(res, 
             list(N_u_lower = nrow(ucomb),
                  uunits = ucomb$unit,
                  uyears = ucomb$year,
                  uyears_full = ucomb$uyears_full,
                  umonths = ucomb$month,
                  starts = starts,
                  ends = ends,
                  map_2_u = map_2_u,
                  map_2_year = map_2_year_obs))
    
  }
  return(res)
}

#' @title Run level
#' @description Determines the GADM level at which to run the model (either level 1
#' or 2) depending on the availability of spatial data
#'
#' @param run_data data on which to run the model
#' @param gadm_list list of GADM units
#' @param min_month_count minimum number of distinct available months at the GDAM 2 level to be considered in the computation of spatial coverage
#' @param mean_cov_thresh threshold on the mean spatial coverage
#' @param max_cov_thresh threshold on the maximal spatial coverage
#'
#' @return the run level
runLevel <- function(country_iso3,
                     run_data,
                     cholera_directory = "./",
                     min_month_count = 4,
                     mean_cov_thresh = .1,
                     max_cov_thresh = .4,
                     verbose = T) {
  
  gadm_data <- getGadmList(country_iso3) 
  
  # Get the fraction of admin units covered at each level for each month
  stats_spatial_coverage <- run_data %>%
    dplyr::filter(purrr::map_dbl(run_data$months, ~ length(.)) < min_month_count) %>%
    dplyr::group_by(gadm_lev, year, month_left) %>%
    dplyr::summarise(frac = sum(unique(gid) %in% gadm_data$gid[gadm_data$gadm_lev == gadm_lev[1]]) /
                       sum(gadm_data$gadm_lev == gadm_lev[1])) %>%
    dplyr::group_by(gadm_lev) %>%
    dplyr::summarise(
      mean = mean(frac),
      max = max(frac)
    )
  
  l2_ind <- which(stats_spatial_coverage$gadm_lev == 2)
  
  if (length(l2_ind) == 0) {
    run_level <- 1
  } else {
    # Choose the level at which to run
    if (stats_spatial_coverage$mean[l2_ind] < mean_cov_thresh & stats_spatial_coverage$max[l2_ind] < max_cov_thresh) {
      run_level <- 1
    } else {
      run_level <- 2
    }
  }
  
  if (run_level == 2 & country_iso3 != "NGA") {
    # If lev 1 admin unit sizes is small enough use admin 1 lev
    gadm_areas <- getGadmSf(country = country_iso3,
                            gadm_lev = 1) %>% 
      # Area in sqkm
      mutate(area = sf::st_area(geom) %>% as.numeric()/1e6)
    
    mean_areas <- getMeanGadmAreas()
    
    admin1_area_thresh <-  mean_areas$mean_area[mean_areas$gadm_lev == 1] * 1.1
    
    if (mean(gadm_areas$area) < admin1_area_thresh) {
      if (verbose) {
        cat("Selecting runlevel 1 despite having sufficient data at level 2 based on area threshold \n")
      }
      run_level <- 1
    }
  }
  
  return(run_level)
}

# Initatial values --------------------------------------------------------

#' @title Draw initial values
#' @description Function to draw initial values
#'
#' @param n_chain number of chains to us
#' @param data the data passed to Stan
#' @param Q_inv_full the covarianvce matrix
#' @param variant the variant of the model
#' @param random 
#' @param res_file
#'
#' @return return
drawInitValues <- function(n_chain, 
                           data, 
                           Q_inv_full, 
                           variant = NULL, 
                           random = T, 
                           res_file = NULL,
                           skip = F) {
  
  if (!random & is.null(res_file)) 
    stop("Need to specify a results file")
  
  mixture <- str_detect(variant, "mixture")
  offset <- str_detect(variant, "offset")
  
  inv_logit <- function(x) 1/(1+exp(-x))
  mixture_pars <- c("betas_g2", "betas_diff", "betas_other_g2", "lambdas_raw", "tau_lambdas")
  model_pars <- c("tau", "tau_betas", "sigma_etas", "xis", "etas_tilde", "theta", "phi", "rho", "betas_g1") %>% 
    {if (mixture) {c(., mixture_pars)} else {.}}
  
  
  Q_beta_full_inv <- makeQbetaFullInv()
  
  # Set vector of probabilities
  prob_non_zero <- rep(NA, data$N_obs)
  
  init_list <- list()
  max_attemps <- 10000
  iter <- 0
  
  while (length(init_list) < n_chain) {
    
    if(iter > max_attemps) stop("Reached maximum number of attempts")
    
    if (iter %% 10 == 0) {
      cat("i: ", iter, "\n")
      phi <- MASS::mvrnorm(1, rep(0, data$N_units), as.matrix(Q_inv_full))
    }
    
    init_pars <- list(
      tau = rgamma(1, 1, .1),
      sigma_etas =  1/sqrt(rgamma(1, 1, .1)),
      etas_tilde = matrix(rnorm(data$N_years * data$N_eta_groups), ncol = data$N_eta_groups),
      sub_xis = rnorm(data$N_lev-1, -min(7,data$N_units/15), 1),
      rho = rbeta(1, 1, 1),
      theta = rnorm(data$N_units),
      phi = phi
    )
    
    if (variant != "null") {
      init_pars <- append(
        init_pars, 
        list(tau_betas = rgamma(1, 1, .1))
      )
      if(mixture) {
        betas.list <- list(betas_g1 = generateBetas(Q_beta_full_inv, 1), 
                           betas_g2 = generateBetas(Q_beta_full_inv, 2))
      } else {
        betas.list <- list(betas_g1 = MASS::mvrnorm(1, rep(0, 12), as.matrix(Q_beta_full_inv)))
      }
    } else {
      betas.list <- list(betas_g1 = rep(0, 12))
    }
    
    sigma <- 1/sqrt(init_pars$tau)
    convolved_re <- sqrt(1 - init_pars$rho) * init_pars$theta + sqrt(init_pars$rho / data$scaling_factor) * init_pars$phi;
    
    test <- T
    
    for (betas in betas.list) {
      for (g in 1:data$N_groups) {
        for (j in 1:data$N_eta_groups) {
          etas <- init_pars$etas_tilde[, j] * init_pars$sigma_etas
          ind_upper <- unique(data$map_2_obs[data$map_2_lev>1])
          if (length(ind_upper) > 0) {
            for (i in ind_upper) {
              inds <- which(data$map_2_obs == i)
              lev <- data$map_2_lev[inds][1]
              xi <- ifelse(lev > 1, init_pars$xis[lev-1], 0)
              tmp <- 1
              for (ind in inds) {
                if (data$map_2_units[ind] == 0) {
                  for (k in 1:data$N_units) {
                    tmp = tmp * (1-inv_logit(xi + betas[data$months[ind]] + etas[data$years[ind]] + convolved_re[k] * sigma));
                  }
                } else {
                  tmp = tmp * (1-inv_logit(xi + betas[data$months[ind]] + etas[data$years[ind]] + convolved_re[data$map_2_units[ind]] * sigma));
                }
              }
              prob_non_zero[i] <- 1 - tmp
            }
            test <- test & sum(prob_non_zero[!is.na(prob_non_zero)] > (1-1e-10) | prob_non_zero[!is.na(prob_non_zero)] < 1e-10) == 0
          } 
        }
      }
    }  
    
    if (test | skip) {
      if(mixture) {
        init_list[[length(init_list) + 1]] <- append(
          init_pars, 
          list(betas_g1 = betas.list$betas_g1,
               betas_other_g2 = betas.list$betas_g2[-c(1,7)],
               betas_diff = log(abs(c(betas.list$betas_g2[c(1,7)] - betas.list$betas_g1[c(1,7)]))),
               ind_betas_diff = c(1, 7),
               lambdas_raw = MASS::mvrnorm(1, rep(0, data$N_units), as.matrix(Q_inv_full)),
               tau_lambdas =  as.array(rgamma(1, 1, .1)))) 
      } else if (variant != "null") {
        init_list[[length(init_list) + 1]] <- append(init_pars, 
                                                     list(betas_g1 = betas.list$betas_g1,
                                                          ind_betas_diff = c(1, 1)))
      } else {
        init_pars$etas_tilde <- init_pars$etas_tilde[, 1]
        init_pars$ind_betas_diff <- c(1, 1)
        init_list[[length(init_list) + 1]] <- init_pars
      }
      cat("Found initial values\n")
    }
    iter <- iter + 1
  }
  return(init_list)
}  

initPars <- function(n_chain, 
                     data_mapping,
                     data,
                     cntry.sf,
                     scaling_factor, 
                     variant = NULL,
                     cntry_iso3,
                     run_level,
                     Q_inv_full,
                     do_plot = F,
                     plot_file = NULL, 
                     ind_betas_diff = NULL) {
  
  variant <- str_c(variant, collapse = "_")
  
  if (is.null(variant))
    stop("Specify model variant")
  
  # Initialize with random start for null model
  if (variant  == "null") {
    init_part_rnd <- drawInitValues(n_chain, 
                                    data, 
                                    Q_inv_full, 
                                    variant = variant)
    
    cat("-- Initialized with Random draw fit\n")
    return(init_part_rnd)
  } 
  
  if (!str_detect(variant, "mixture")) {
    # Try to initialize with base fit
    try(init_from_base <- initValuesFromBase(opt = opt,
                                             cntry_iso3 = cntry_iso3,
                                             run_lev = run_lev,
                                             n_chain = n_chain,
                                             data = data, 
                                             variant = variant))
    
    # Return if successful
    if (!is.null(init_from_base)) {
      cat("-- Initialized with Base model fit\n")
      return(init_from_base)
    } else {
      cat("-- Failed base model fit, trying GAM \n")
    }
  } else {
    cat("-- Mixture-type model, trying GAM directly \n")
  }
  
  # Initialize with GAM
  init_pars_glm <- try(initialParamGLM(n_chain = n_chain, 
                                       data_mapping = data_mapping,
                                       data = data,
                                       cntry.sf = cntry.sf,
                                       scaling_factor = scaling_factor, 
                                       variant = variant,
                                       cntry_iso3 = cntry_iso3,
                                       run_level = run_level,
                                       do_plot = do_plot,
                                       plot_file = plot_file,
                                       ind_betas_diff = ind_betas_diff),
                       silent = F)
  
  if (!inherits(init_pars_glm, "try-error")) {
    cat("-- Initialized with GAM fit\n")
    return(init_pars_glm)
    
    # } else if (str_detect(variant, "mixture")) {
    #   warning("Failed GAM initialization for mixture-type model")
  } else {
    # Initialize with random start
    init_part_rnd <- drawInitValues(n_chain, 
                                    data, 
                                    Q_inv_full, 
                                    variant = variant)
    
    cat("-- Initialized with Random draw fit\n")
    return(init_part_rnd)
  }
}

#' @title Try to load initial values from base fit
#'
#' @param n_chain number of chains to us
#' @return return
initValuesFromBase <- function(opt,
                               cntry_iso3,
                               run_lev,
                               n_chain,
                               data, 
                               variant = NULL) {
  
  if (is.null(variant))
    stop("Specify model variant")
  
  mixture <- str_detect(variant, "mixture")
  offset <- str_detect(variant, "offset")
  
  # Make base result file
  opt2 <- opt
  opt2$variant <- "base"
  base_res_file <- makeResultsFile(opt2, cntry_iso3, run_level, drop_multi_yr = FALSE)
  
  if (!file.exists(base_res_file)) {
    cat("-- Could not find base results file \n")
    return(NULL)
  } else {
    cat("-- Found base results file \n")
  }
  
  # Load Base result file
  chol_stanfit <- readRDS(base_res_file)
  chol_stanfit <- chol_stanfit$chol_stanfit
  
  if (length(chol_stanfit@stan_args) == 0) {
    cat("-- Base file does not contain samples \n")
    return(NULL)
  }
  
  # Total number of sample
  n_samp <- chol_stanfit@stan_args[[1]]$iter - chol_stanfit@stan_args[[1]]$warmup
  # random number of samples
  inds <- sample(1:n_samp, n_chain)
  
  # Names of base model parameters
  model_pars_base <-  c("tau", "tau_betas", "sigma_etas", "xis", "etas_tilde", "theta", "phi", "rho", "betas_g1")
  
  # Initialize list
  init_list <- list()
  
  # Get samples from the base model
  pars <- rstan::extract(chol_stanfit, pars = model_pars_base)
  
  for (ind in inds) {
    ipars <- list()
    for (p in model_pars_base) {
      if (length(dim(pars[[p]])) == 3) {
        ipars <- append(ipars, set_names(list(pars[[p]][ind, ,]), p))
      } else if(length(dim(pars[[p]])) == 2) {
        ipars <- append(ipars, set_names(list(pars[[p]][ind, ]), p))
      } else {
        ipars <- append(ipars, set_names(list(pars[[p]][ind]), p))
      }
    }
    
    if (mixture) {
      # Generate betas for second group
      betas_g2 <- generateBetas(makeQbetaFullInv(), 2)
      
      betas_2 <- list(betas_other_g2 = betas_g2[-c(1,7)],
                      betas_diff = log(abs(c(betas_g2[c(1,7)] - ipars$betas_g1[c(1,7)]))))
      
      ipars <- append(ipars, betas_2)
    }
    
    if (offset) {
      # Regenerate yearly random effects
      etas_tilde <- matrix(rnorm(data$N_years * data$N_groups), ncol = data$N_groups)
      ipars$etas_tilde <- etas_tilde
    }    
    
    if (mixture) {
      ipars$etas_tilde <- matrix(ipars$etas_tilde, ncol = 1)
      if (offset) {
        ipars$etas_tilde <- cbind(ipars$etas_tilde,  ipars$etas_tilde)
      }
    } else {
      ipars$etas_tilde <- matrix(ipars$etas_tilde, ncol = 1)
    }
    
    if (!is.null(ipars$xis)) {
      ipars$xis <- as.array(ipars$xis)
    }
    init_list[[length(init_list) + 1]] <- ipars
  }
  rm("chol_stanfit")
  cat("-- Starting with ")
  return(init_list)
}


#' @title Initialize parameters
#'
#' @description Initial guess for parameters using a GLM
#'
#' @param n_chain number of chains to us
#' @param data the data passed to Stan
#' @param Q_inv_full the covarianvce matrix
#' @param variant the variant of the model
#' @param random 
#' @param res_file
#'
#' @return the list of parameters
initialParamGLM <- function(n_chain, 
                            data_mapping,
                            data,
                            cntry.sf,
                            scaling_factor, 
                            variant = NULL,
                            cntry_iso3,
                            run_level,
                            do_plot = F,
                            plot_file = NULL, 
                            ind_betas_diff = NULL) {
  
  cat("- Running GLM for initial parameter values \n")
  
  mixture <- str_detect(variant, "mixture")
  offset <- str_detect(variant, "offset")
  
  # Setup data for GAM
  ind_runlev <- which(data_mapping$map_2_lev_mapping == 1)
  
  dat <- tibble (
    y = data_mapping$y[data_mapping$map_2_obs[ind_runlev]],
    n = data_mapping$n_replicas[data_mapping$map_2_obs[ind_runlev]],
    month = data_mapping$months[ind_runlev],
    year = data_mapping$years[ind_runlev],
    unit = data_mapping$map_2_units[ind_runlev],
  ) %>% 
    dplyr::mutate(
      month = factor(month, levels = c(1:12)),
      year = factor(year),
      gadm_id = data_mapping$u_unit[unit],
      unit = factor(unit),
      failures = n - y,
      full_year = data_mapping$u_years[year]
    )
  
  # Compute coordinates of admin unit centroids
  centroids <- sf::st_centroid(cntry.sf)
  xy_coords <- sf::st_coordinates(centroids) %>% 
    as.data.frame() %>%
    dplyr::mutate(id = centroids$gid)
  
  # Aggregate data by gadm id and months
  dat_bymonth <- dat %>% 
    dplyr::group_by(unit, gadm_id, month) %>% 
    dplyr::summarise(pos = sum(y)/sum(n)) %>% 
    dplyr::group_by(gadm_id, unit) %>% 
    tidyr::complete(month = factor(1:12), fill = list(pos = 0)) %>%
    dplyr::mutate(pos_scaled = (pos - mean(pos, na.rm = T))/sd(pos, na.rm = T),
                  pos_scaled = ifelse(is.nan(pos_scaled), 0, pos_scaled)) %>% 
    dplyr::inner_join(xy_coords, by = c("gadm_id" = "id")) %>% 
    pivot_wider(values_from = "pos_scaled", names_from = "month", id_cols = c("unit", "gadm_id", "X", "Y")) %>% 
    ungroup() %>% 
    na.omit() %>%
    dplyr::mutate(X = (X-min(X))/(max(X)-min(X)),
                  Y = (Y-min(Y))/(max(Y)-min(Y)))
  
  dat_summary <- dat %>% 
    dplyr::group_by(unit, gadm_id) %>% 
    dplyr::summarise(pos = sum(y)/sum(n),
                     n_tot =sum(n)) 
  
  # SF object
  cntry.geom <- select(cntry.sf, gid, geom) %>% 
    .[order(.$gid), ] 
  
  cntry.sf <- cntry.sf %>% 
    as.data.frame() %>% 
    select(-geom) %>% 
    left_join(dat_bymonth, by = c("gid" = "gadm_id")) %>%  
    left_join(dat_summary, by = c("gid" = "gadm_id", "unit" = "unit"))
  
  cntry.sf$unit <- purrr::map_dbl(cntry.sf$gid, ~ which(data_mapping$u_unit == .))
  cntry.geom$unit <- purrr::map_dbl(cntry.geom$gid, ~ which(data_mapping$u_unit == .))
  
  if (mixture) {
    # Compute grouping with kmeans on monthly fraction of occurrence
    group_res <- kmeans(dat_bymonth %>% select(-unit, -gadm_id), centers = 2, nstart = 1e2)
    dat_bymonth$grp <- group_res$cluster
    dat <- dat %>% left_join(dat_bymonth %>% select(unit, grp)) %>% dplyr::mutate(grp = factor(grp))
    
    # Fit GLM model using kmeans grouping 
    eq <- ~ month:grp + year + unit - 1
    res <- fitGLMNET(dat %>% dplyr::filter(!is.na(grp)), eq)
    
    # Get monthly coefficients
    betas <- cbind(res[str_detect(names(res), "month")][1:12], res[str_detect(names(res), "month")][13:24])
    mean_betas <- apply(betas, 2, mean)
    betas <- betas - rep(mean_betas, each = 12)
    
    beta_mean <-  mean_betas[which.max(abs(mean_betas))]
    beta_offset <- mean_betas[which.min(abs(mean_betas))] - mean_betas[which.max(abs(mean_betas))]
    grp_high_betas <- which.max(abs(mean_betas))
    
    # Get yearly random effects
    etas <- cbind(res[str_detect(names(res), "year")] + beta_mean, res[str_detect(names(res), "year")] + beta_mean)
    
    # Get spatial random effects
    ind_grp1 <- str_replace_all(names(res)[str_detect(names(res), "unit")], "unit", "") %>% 
      as.numeric() %>%
      purrr::map_lgl(~ . %in% dat_bymonth$unit[dat_bymonth$grp == 1])
    sp_re1 <- res[str_detect(names(res), "unit")][ind_grp1]
    sp_re2 <- res[str_detect(names(res), "unit")][!ind_grp1]
    
    # Combine spatial random effects correcting for offset between groups
    if (grp_high_betas == 2) {
      sp_re1 <- sp_re1 + beta_offset
    } else {
      sp_re2 <- sp_re2 + beta_offset
    }
    sp_re <- c(sp_re1, sp_re2)
    
    # Define which group corresponds to group 1 in the Stan code 
    # By convention group 1 has larger values closer to July and lower in January
    # Groups are determined based on the month of maximal discrepancy between estimates
    
    inds <- c(which.max(betas[, grp_high_betas]), 
              which.min(betas[, 1] - betas[, 2]))
    
    ind_higherg1 <- inds[which.min(abs(7 - inds))] 
    grp1_ind <- which.max(betas[ind_higherg1,])
    
    all_betas_diff <- betas[, grp1_ind] - betas[, -grp1_ind]
    # First index corresponds to month where group 2 is larger than group 1
    if (is.null(ind_betas_diff)) {
      if (grp1_ind == grp_high_betas) {
        ind_betas_diff <- c(which.min(all_betas_diff),
                            which.max(betas[, grp1_ind]))
      } else {
        ind_betas_diff <- c(which.max(betas[, grp_high_betas]),
                            which.max(all_betas_diff))
      }
    }
    
    if (ind_betas_diff[1] == ind_betas_diff[2]) {
      # inds[2] <- inds[1] + 1
      stop("GLM found same initial ind value")
    }
    
    betas_diff <- all_betas_diff[ind_betas_diff]
    
    # Reorder betas and etas matrices
    betas <- cbind(betas[, grp1_ind], betas[, -grp1_ind])
    etas <- cbind(etas[, grp1_ind], etas[, -grp1_ind])
    
    if (!offset) {
      etas <- etas[,1]
    }
    
    # Spatial RE model for mixture coefficients using parameterization in Riebler et al. 2016 to fit with INLA
    rep_grp <- 2
    grp_data <- dat_bymonth %>% 
      dplyr::mutate(grp = as.numeric(grp == grp1_ind) * rep_grp, 
                    unit = as.numeric(as.character(unit)),
                    n = rep_grp) %>% 
      tidyr::complete(unit = cntry.sf$unit) %>% 
      replace_na(list(grp = 0, n = 1))
    
    spatial_grp <- getSpatialRE(cntry.sf = cntry.geom, 
                                data = grp_data %>% rename(y = grp), 
                                grouping = T,
                                scaling_factor = scaling_factor,
                                cntry_iso3 = cntry_iso3,
                                variant = variant,
                                run_level = run_level)
    
  } else {
    if (variant %in% c("base", "offset")) {
      
      # Fit GLM model
      eq <- ~ month + year + unit - 1
      coefs <-  fitGLMNET(dat, eq)
      
      # Get monthly coefficients
      betas <- coefs[str_detect(names(coefs), "month")]
      mean_betas <- mean(betas)
      betas <- betas - mean_betas
      ind_betas_diff <- NULL
      
    } else if (variant == "null")  {
      eq <- ~ year + unit - 1
      coefs <- fitGLMNET(dat, eq)
    }
    
    # Get yearly random effects
    etas <- coefs[str_detect(names(coefs), "year")]
    # Get spatial random effects
    sp_re <- coefs[str_detect(names(coefs), "unit")]
    mean_spre <- mean(sp_re)
    sp_re <- sp_re - mean_spre
    
    # Get yearly random effects
    etas <- etas + mean_betas + mean_spre
  }
  
  if (is.null(dim(etas))) {
    eta_names <- names(etas)
  } else {
    eta_names <- rownames(etas)
  }
  
  # Expand etas vector to all available years
  etas_padded <- etas %>% as.data.frame() %>% 
    dplyr::mutate(year = str_replace_all(eta_names, "year", ""),
                  year = as.integer(year)) %>%
    tidyr::complete(year = seq_along(data_mapping$u_years)) %>% 
    as.matrix() %>%  .[, -1] %>% 
    as.matrix(ncol = 1)
  
  # Get unique years for which no information is available
  which_missing_etas <- which(is.na(etas_padded[, 1]))
  sd_etas <- ifelse(length(etas) > 1, sd(etas), abs(etas)/3)
  rnd_etas <- rnorm(length(which_missing_etas), 0, sd_etas)
  for (i in 1:ncol(etas_padded)) {
    etas_padded[which_missing_etas, i] <- rnd_etas
  }
  
  # Spatial RE using parametrization in Riebler et al. 2016 to fit with INLA
  spat_dat <- tibble(y = sp_re, 
                     idx = as.numeric(str_replace_all(names(sp_re), "unit", ""))) %>% 
    tidyr::complete(idx = unique(cntry.sf$unit)) %>% 
    tidyr::replace_na(list(y = 0))
  
  spatial_re <- getSpatialRE(cntry.sf = cntry.geom, 
                             data = spat_dat %>% rename(unit = idx), 
                             grouping = F,
                             scaling_factor = scaling_factor,
                             cntry_iso3 = cntry_iso3,
                             variant = variant,
                             run_level = run_level)
  
  # Set initial parameter draws for Stan
  ref_params <- list(
    tau = spatial_re$tau,
    tau_beta = 1/sd(betas)^2,
    sigma_etas = sd_etas,
    etas_tilde = (etas_padded/sd_etas)[, 1:data$N_eta_groups, drop = F], # standardize etas
    rho = spatial_re$rho,
    theta = spatial_re$theta,
    phi = spatial_re$phi
  )
  
  if (mixture) {
    ref_params <- c(
      ref_params, 
      list(
        betas_g1 = betas[, grp1_ind],
        betas_diff = log(abs(betas_diff)), # get differences in January and July
        betas_other_g2 = betas[-ind_betas_diff, -grp1_ind],
        ind_betas_diff = ind_betas_diff,
        lambdas_raw = spatial_grp$phi,
        theta_lambdas_raw = spatial_grp$theta,
        rho_lambdas = as.array(spatial_grp$rho),
        tau_lambdas = as.array(spatial_grp$tau)
      ))
  } else if (!(variant == "null")) {
    ref_params <- c(ref_params, list(betas_g1 = betas))
  }
  
  # Parameters with random initial values (xis)
  ref_params <- c(ref_params, 
                  list(sub_xis = as.array(rnorm(data$N_lev-1, -min(7,max(data$N_units/15, 5)), 1))))
  
  # Plot results of GLM initialization
  if(do_plot) {
    plotPars(pars = ref_params, 
             cntry.sf = cntry.geom, 
             variant =  variant, 
             plot_file = plot_file, 
             scaling_factor = scaling_factor, 
             data_mapping = data_mapping)
  }
  
  if (n_chain == 1) {
    return(list(ref_params))
  } else {
    # Initialize list of initial parameters
    init_list <- list()
    
    if (mixture) {
      ind_betas_diff <- as.integer(ref_params$ind_betas_diff)
      ref_params <- ref_params[-which(names(ref_params) == "ind_betas_diff")]
    } else {
      ind_betas_diff <- c(1, 1)
    }
    
    # Perturb parameters
    for (i in 1:n_chain) {
      ipars <- lapply(
        ref_params, 
        function(x) x + rnorm(n = length(x), mean = 0, sd = mean(abs(x))/5)
      )
      init_list[[length(init_list) + 1]] <- c(applyBounds(ipars),
                                              list(ind_betas_diff = ind_betas_diff))
    }
    return(init_list)
  }
}


# Stan model --------------------------------------------------------------

getNumThreads <- function(config) {
  if (config$stan$run_threads) {
    # Check multi-threading
    n_threads <- round(as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))/config$stan$n_chains) 
    
    cat("-- Running with", n_threads, "threads per chain \n")
  } else {
    n_threads <- NULL
    cat("-- Running with single threading \n")
  }
  n_threads
}

getStanFile <- function(config,
                        run_level) {
  
  # Compile the model
  if (any(str_detect(config$variant, "null"))) {
    stan_file <- "binomial_bym2_fast_null_speedup.stan"
  } else if (any(str_detect(config$variant, "offset"))) {
    stan_file <- "binomial_bym2_fast_convolved_lambda_flexmixture_speedup_map_rect.stan"
  } else {
    stan_file <- "binomial_bym2_fast_convolved_lambda_flexmixture_speedup.stan"
  }
  
  if (config$country_iso3 %in% c("NGA", "COD")) {
    if (run_level == 2) {
      if (config$variant == "null") {
        stan_file <- "binomial_bym2_fast_null_speedup_map_rect_v2.stan"
      } else {
        stan_file <- "binomial_bym2_fast_convolved_lambda_flexmixture_speedup_map_rect_v2.stan"
      }
    } 
  }
  
  str_c("analysis/stan/", stan_file)
}


getStanModel <- function(config,
                         run_level) {
  set_cmdstan_path()
  
  stan_file <- getStanFile(config, run_level)
  
  if (config$stan$run_threads) {
    chol_stan <- cmdstan_model(stan_file, 
                               cpp_options = list(stan_threads = TRUE),
                               quiet = FALSE,
                               force_recompile = F)
  } else {
    chol_stan <- cmdstan_model(stan_file,
                               quiet = FALSE,
                               force_recompile = F)
  }
  chol_stan
}

# Logs --------------------------------------------------------------------

#' @title Get error type
#'
#' @description gets the type of error from the log
#'
#' @param err_msg the error message
#'
#' @return a string with the error type
#' 
getErrorType <- function(err_msg) {
  if (str_detect(err_msg, "Killed")) {
    err <- "Memory limit"
  } else if (str_detect(err_msg, "foreign")) {
    err <- "All obs are 0"
  } else if (str_detect(err_msg, "tidyverse")) {
    err <- "No data"
  } else if (str_detect(err_msg, "CANCELLED")) {
    err <- "Time limit"
  }   else if (str_detect(err_msg, "RUNNING")) {
    err <- "Still running"
  } else if (str_detect(err_msg, "already")) {
    err <- "Already ran"
  } else if (str_detect(err_msg, "diver")) {
    err <- "Ran with divergences"
  } else if (err_msg == "") {
    err <- "No error"
  } else if (err_msg == "Failed initialization") {
    err <- "Failed initialization"
  } else if (err_msg == "Failed sampling") {
    err <- "Failed sampling"
  } else {
    err <- "no error recognised"
  }
  return(err)
}

#' @title Model run status
#'
#' @description Gets the status of model runs for a given country, run level and
#' identifier
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param paths what paths to look into
#' @param redo_logs whether to redo log parsing
#'
#' @return a tibble with the run status for each model
#' 
getModelRunStatus <- function(country,
                              run_level,
                              identifier,
                              time_left = NULL,
                              time_right = NULL,
                              gadm_lev = setGadmLev(country, run_level),
                              paths = getResPaths(),
                              log_dir = getLogDir(),
                              redo_logs = F,
                              redo_data = F,
                              verbose = F,
                              ...) {
  
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo_data,
                                 verbose = verbose) 
  }
  
  logs <- runLogParsing(log_dir = log_dir,
                        redo = redo_logs) %>% 
    dplyr::group_by(identifier, country, model, run_level) %>%
    dplyr::filter(start_time == max(start_time))  %>% 
    dplyr::summarise(ran = ifelse(ran, "ran",
                                  ifelse(any(status == "running"), "running", 
                                         ifelse(err_type == "Already ran", "ran", err_type))),
                     percent_complete = max(stan_chain_iter_min)/2000) %>% 
    ungroup() %>% 
    dplyr::mutate(run_level = as.numeric(run_level))
  
  # Result files to parse
  available_results <- try(getResFiles(country = country,
                                       run_level = run_level, 
                                       identifier = identifier, 
                                       paths = paths,
                                       time_left = time_left,
                                       time_right = time_right,
                                       gadm_lev = gadm_lev) %>% 
                             purrr::map_df(function(x) {
                               parseResFilename(x) %>% 
                                 mutate(file = x,
                                        fit_avail = purrr::map_lgl(file, ~checkFile(.)))
                             })) 
  
  if(!inherits(available_results, "try-error")) {
    
    if (verbose) {
      if (all(purrr::map_lgl(getAllModels(), ~ . %in% available_results$variant[available_results$ran == "ran"]))) {
        cat("All model runs available for", country, run_level, identifier)
      } 
    }
    
    logs <- logs %>% 
      dplyr::right_join(available_results %>% 
                          dplyr::mutate(run_level = as.numeric(run_level)) %>% 
                          dplyr::right_join(
                            tidyr::expand(available_results, 
                                          variant = getAllModels(),
                                          country = country,
                                          run_level = run_level,
                                          identifier = identifier)),
                        by = c("country" = "country", 
                               "run_level" = "run_level",
                               "model" = "variant",
                               "identifier" = "identifier"))
    
  } else {
    logs <- dplyr::filter(logs, 
                          identifier == .env[["identifier"]],
                          country == .env[["country"]],
                          run_level == .env[["run_level"]])
  }
  return(logs)
}


#' @title Parse mars log
#'
#' @description Parse a model run log
#'
#' @param log_file the log file
#'
#' @return a tibble with the parsed log
parseLog <- function(log_file) {
  log_txt <- readLines(log_file)
  
  # Timing 
  times <- log_txt[str_detect(log_txt, "[0-9]{2}:[0-9]{2}:[0-9]{2}")] %>% 
    str_subset("INLA", negate = T) 
  
  tz <- str_extract(times, "[A-Z]{3}")[1]
  
  times <- times %>% 
    str_replace_all(tz, "") %>% 
    str_replace_all(" ", "-") %>% 
    str_replace_all("--", "-") %>% 
    str_split("-") %>% 
    map_chr(~paste(11, paste(.[3:5], collapse = " ")))
  
  times <- strptime(times, format = "%m %d %H:%M:%S %Y", tz = tz)
  
  if (length(times) > 1) {
    whole_code_time <- difftime(times[2], times[1], units = "hours") %>% 
      formatC(digits = 2) %>% 
      str_c(" hours") %>% 
      str_trim()
    
  } else {
    whole_code_time <- difftime(Sys.time(), times[1], units = "hours") %>% 
      formatC(digits = 2) %>% 
      str_c(" hours") %>% 
      str_trim()
  }
  
  if (!any(str_detect(log_txt, "Running:\t"))) {
    return(NULL)  
  }
  
  # Get country name
  cntry <- str_extract(log_txt[str_detect(log_txt, "Running:\t")], "[A-Z]{3}") 
  
  # Get gadm run level
  run_level <- str_extract(log_txt[str_detect(log_txt, "GADM level")], "[1-2]") 
  
  # Get model
  model_name <- str_split(log_txt[str_detect(log_txt, "with model")] %>% str_trim(), " ")[[1]] %>% last()
  
  # Get cholera definition
  identifier <- str_split(log_txt[str_detect(log_txt, "with cholera")] %>% str_trim(), " ")[[1]] %>% last()
  
  # Get job ID
  job_id <- str_extract(log_file, "(?<=_)[0-9]+(?=_)")
  
  # Get array ID
  array_id <-  str_extract(log_file, "(?<=_)[0-9]+(?=\\.)") 
  
  
  stan_chain_iter <- map_df(1:4, function(x) {log_txt %>%
      str_subset(str_c("Chain ", x)) %>% 
      str_extract_all("(?<= )[0-9]+(?= /)") %>% 
      unlist() %>% 
      as.numeric() %>% 
      unique() -> iters
    
    tibble(chain = x, max_iter = max(iters))}) 
  
  # Did the code run
  ran <- any(str_detect(log_txt, "End of script")) & 
    !any(str_detect(log_txt, "Execution halted")) &
    !any(str_detect(tail(log_txt, 10), "Error") & !str_detect(tail(log_txt, 10), "parallel"))
  
  
  if (!ran) {
    if (str_detect(last(log_txt), "Chain") | last(log_txt) == "" ) {
      err_msg <- "RUNNING"
    } else {
      err_ind <- which(str_detect(tail(log_txt, 20), "Error|Killed") %>% unlist()) + length(log_txt) - 20
      if (length(err_ind) == 0) 
        err_ind <- length(log_txt)
      err_msg <- log_txt[err_ind:length(log_txt)] %>% 
        paste(collapse = "\n")
    }
    stan_chain_times <- ""
    stan_chain_times_num <- list(NA)
  } else {
    
    # Get stan code runtimes in seconds
    stan_chain_times <- log_txt[str_detect(log_txt, "(Total)")] %>% 
      str_extract(., "(?<=  ).*?(?=s)") %>% 
      str_trim() %>% 
      as.numeric()
    
    stan_chain_times_num <- list(stan_chain_times/3600)
    stan_chain_times <- range(stan_chain_times/3600) %>% 
      formatC(digits = 2) %>% 
      paste(collapse = "-") %>% 
      str_trim() %>% 
      str_c(" hours")
    
    
    if (any(str_detect(log_txt, "Initialization failed"))) {
      ran <- F
      err_msg <- "Failed initialization"
    } else if (any(str_detect(log_txt, "divergent"))){
      err_msg <- str_c("ran with ", 
                       str_subset(log_txt, "divergent") %>% 
                         str_extract("(?<=were )[0-9]+(?= divergent)") %>% 
                         .[1], " divergent transitions")
      
    } else if (any(str_detect(log_txt, "' does not contain samples"))){
      err_msg <- "Failed sampling"
      
    } else {
      err_msg <- ""
    }
  }
  
  if (any(str_detect(log_txt, "Gradient evaluation took 0 seconds"))) {
    grad0 <- T
  } else {
    grad0 <- F
  }
  
  return(tibble(country = cntry, 
                start_time = times[1],
                ran = ran, 
                err_msg = err_msg,
                run_level = run_level,
                model = model_name,
                job_id = job_id,
                array_id = array_id,
                whole_code_time = whole_code_time,
                stan_chain_times = stan_chain_times,
                stan_chain_times_num = stan_chain_times_num,
                identifier = identifier,
                grad0 = grad0,
                chain_iter = list(stan_chain_iter$max_iter),
                log_file = log_file)
  )
}

#' @title Run log parsing
#'
#' @description Run log parsing
#'
#' @param log_dir log directory
#'
#' @return a tibble with the parsed logs
#' 
runLogParsing <- function(log_dir = getLogDir(),
                          redo = T) {
  
  err_file <- paste0("generated_data/model_outputs/parsed_logs.rds")
  
  if (file.exists(err_file) & ! redo) {
    errors <- readRDS(err_file)
  } else {
    
    log_files <- dir(log_dir,
                     full.names = T) %>% 
      sort()
    
    errors <- log_files %>% 
      map_df(~parseLog(.))
    
    errors$stan_chain_times_min <- map_dbl(errors$stan_chain_times_num, ~ifelse(any(is.na(.)) | length(.) == 0, as.numeric(NA), min(.)))
    errors$stan_chain_times_max <-  map_dbl(errors$stan_chain_times_num, ~ifelse(any(is.na(.)) | length(.) == 0, as.numeric(NA), max(.)))
    
    errors$stan_chain_iter_min <- map_dbl(errors$chain_iter, ~ifelse(any(is.infinite(.)) | length(.) == 0, as.numeric(NA), min(.)))
    errors$stan_chain_iter_max <-  map_dbl(errors$chain_iter, ~ifelse(any(is.infinite(.)) | length(.) == 0, as.numeric(NA), max(.)))
    
    errors <- errors  %>% 
      mutate(err_type = map_chr(err_msg, getErrorType),
             code_time = as.numeric(str_replace_all(whole_code_time, " hours", "")),
             divergent = str_detect(err_msg, "diver"),
             n_diverge = as.numeric(str_extract(err_msg, "(?<=with )[0-9]+(?= divergent)")),
             status = case_when(ran & (!divergent | n_diverge < 300) & (err_msg != "Failed sampling") ~ "ran",
                                !ran & err_msg == "RUNNING" ~ "running",
                                ran & (divergent & n_diverge > 300) ~ "diverged",
                                ran & err_msg != "failed sampling" ~ "failed sampling",
                                T ~ "failed"))
    
    saveRDS(errors, file = err_file)
  }
  return(errors)
}

# Misc --------------------------------------------------------------------

# Function to apply constraints on parameter values
applyBounds <- function(pars) {
  tol <- 1e-3
  pars[str_detect(names(pars), "rho")] <- pars[str_detect(names(pars), "rho")] %>% 
    pmin(1-tol) %>% 
    pmax(0+tol)
  
  pars[str_detect(names(pars), "tau|sigma")] <- pars[str_detect(names(pars), "tau|sigma")] %>% 
    pmax(0+tol)
  
  if (any(str_detect(names(pars), "lambda"))) {
    pars$rho_lambdas <- as.array(pars$rho_lambdas)
    pars$tau_lambdas <- as.array(pars$tau_lambdas)
  }
  return(pars)
}


classify_date = function(date){
  ifelse(date <= 2,
         "day",
         ifelse(
           (date > 5) & (date <= 8),
           'week',
           ifelse(
             (date > 12) & (date <= 16),
             'biweek',
             ifelse(
               (date > 26) & (date <= 32),
               'month',
               ifelse(
                 (date > 56) & (date <= 64),
                 'bimonth',
                 ifelse(
                   (date > 360) & (date <= 370),
                   'year',
                   ifelse(
                     (date > 720) & (date <= 740),
                     'biyear',
                     classify_date_sub(date)
                   )
                 )
               )
             )
           )
         )
  )
}

classify_date_sub <- function(date){
  ifelse(date <= 2,
         "day",
         ifelse(
           (date > 2) & (date <= 8),
           'Oweek',
           ifelse(
             (date > 8) & (date <= 32),
             'Omonth',
             ifelse(
               (date > 32) & (date <= 185),
               'Osemester',
               ifelse(
                 (date > 185) & (date <= 370),
                 'Oyear',
                 # as.character((floor(as.numeric(date)/365)+1))
                 'multi_year'
               )
             )
           )
         )
  )
}

#' Compute outcome variable
#'
#' @param dfa 1 row dataframe
#' @param thresh 
#' @param incid_tresh 
#'
#' @details 
#'
#' @return an integer with the value of number of occurrences
#'
#' @examples
#' 
computeOutcomeVar <- function(suspected_cases,
                              confirmed_cases,
                              popsum,
                              n_months,
                              cholera_rates_2010_2016,
                              thresh, 
                              incid_tresh = 1e-5,
                              case_thresh = 1) {
  
  if (thresh == "occurrence") {
    sum(pmax(suspected_cases, confirmed_cases, na.rm = T) > 0)
  } else if (thresh == "cases") {
    sum(pmax(suspected_cases, confirmed_cases, na.rm = T)/popsum/n_months > incid_tresh)
  } else if (thresh == "mean_annual_incidence") {
    sum(pmax(suspected_cases, confirmed_cases, na.rm = T)/popsum*12/n_months > cholera_rates_2010_2016)
  } else if (thresh == "annual_incidence") {
    sum(pmax(suspected_cases, confirmed_cases, na.rm = T) > case_thresh)
  } else {
    stop("Uknown threshold to determine cholera presence/absence")
  }
}

fitGLMNET <- function(dat, eq) {
  X <- model.matrix(eq, data = dat) %>% as.matrix()
  res <- glmnet::glmnet(x = X, y = factor(dat$y >= 1), family = "binomial", alpha = 0, intercept = F, nlambda = 20)
  return(coef(res)[, 20])
}

#' @description Function to generate multivariate normals
#' @param Q_beta_full_inv the covariance matrix
#' @return return
generateBetas <- function(Q_beta_full_inv, grp = 1) {
  betas <- MASS::mvrnorm(1, rep(0, 12), as.matrix(Q_beta_full_inv))
  if (grp == 2) {
    while (betas[6] > 0 | betas[7] > 0 | betas[8] > 0 | betas[12] < 0 | betas[1] < 0) {
      betas <- MASS::mvrnorm(1, rep(0, 12), as.matrix(Q_beta_full_inv))
    }
  } else {
    while (betas[6] < 0 | betas[7] < 0 | betas[8] < 0 | betas[12] > 0 | betas[1] > 0) {
      betas <- MASS::mvrnorm(1, rep(0, 12), as.matrix(Q_beta_full_inv))
    }
  }
  return(betas)
}


get_ind_betas_diff <- function(str) {
  
  if(is.null(str))
    return(NULL)
  
  ind_betas_diff <- str_split(str, "_")[[1]] %>% as.numeric()
  
  if (length(ind_betas_diff) != 2)
    stop("ind_betas_diff option not in proper format")
  
  return(ind_betas_diff)
}

getSpatialRE <- function(cntry.sf,
                         data,
                         grouping,
                         scaling_factor,
                         cntry_iso3,
                         variant,
                         run_level) {
  
  
  cntrynb <- poly2nb(cntry.sf, row.names = cntry.sf$unit)
  
  dir.create("generated_data/re_graph", showWarnings = F)
  graph_file <- glue::glue("generated_data/re_graph/re_graph_{cntry_iso3}_{variant}.gra")
  nb2INLA(nb=cntrynb, file = graph_file)
  
  formula <- y ~ -1 + f(unit,
                        model = "bym2",
                        graph=graph_file,
                        scale.model = TRUE, 
                        constr = TRUE, 
                        hyper=list(
                          phi = list(
                            prior = "pc",
                            param = c(.5, 2/3),
                            initial = -1), 
                          prec = list(
                            prior = "pc.prec",
                            param = c(1, .9),
                            initial = 5)))
  if (grouping) {
    r <- try(inla(formula = formula,
                  Ntrials = n,
                  data = data,
                  family = "binomial", 
                  control.predictor = list(compute=TRUE)))
  } else {
    r <- try(inla(formula = formula, 
                  data = data,
                  family = "gaussian", 
                  control.predictor = list(compute=TRUE)))
    
  }
  
  if (!inherits(r, "try-error")) {
    # Unpack results
    conv_re <- as_tibble(r$summary.random$unit) %>% slice(1:nrow(cntry.sf)) %>% .[["mean"]]
    phi <- as_tibble(r$summary.random$unit) %>% slice((nrow(cntry.sf)+1):nrow(r$summary.random$unit)) %>% .[["mean"]]
    tau <- r$summary.hyperpar[ifelse(grouping, 1, 2),1]
    rho <- r$summary.hyperpar[ifelse(grouping, 2, 3),1]
    theta <- (sqrt(tau) * conv_re - sqrt(rho/scaling_factor) * phi)/(sqrt(1-rho))
  } else {
    # Unpack results
    conv_re <- rnorm(nrow(cntry.sf))
    phi <- rnorm(nrow(cntry.sf))
    tau <- .5
    rho <- .3
    theta <- (sqrt(tau) * conv_re - sqrt(rho/scaling_factor) * phi)/(sqrt(1-rho))
  }
  
  res <- list(conv_re = conv_re,
              phi = phi,
              theta = theta,
              rho = rho,
              tau = tau)
  return(res)
}

makeCholeraThresh <- function(opt) {
  if (opt$thresh == "cases") {
    str_c(opt$case_thresh, "cases")
  } else {
    opt$thresh
  }
}

makeQbetaFullInv <- function() {
  #Build the adjacency matrix using INLA library functions
  beta.graph <- rbind(diff(diag(12), diff = 1), c(1, rep(0, 10), -1))#sparseMatrix(i= seq(12), c(12, seq(11)), x=1,symmetric=F)
  Q_beta <- t(beta.graph) %*% beta.graph
  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_beta_pert <- Q_beta + Diagonal(12) * sqrt(.Machine$double.eps)
  Q_beta_full_inv <- INLA:::inla.ginv(inla.scale.model(Q_beta_pert, constr=list(A = matrix(1,1,12),e=0)))
  return(Q_beta_full_inv)
}

#' @title Pad years
#' @description Pads years to use with offset
#'
#' @param obs_years vector f observed years
#'
#' @return a vector with padded year ids
padYears <- function(obs_years) {
  u_years <- sort(unique(obs_years))
  has_next <- purrr::map_lgl(u_years, ~ (.+1) %in% u_years)
  padded_years <- sort(c(u_years, u_years[!has_next] + 1))
  return(padded_years)
}

plotPars <- function(pars,
                     cntry.sf, 
                     variant, 
                     plot_file, 
                     scaling_factor, 
                     data_mapping) {
  
  library(gridExtra)
  
  # Function to plot initial parameter values
  betas <- tibble(g1 = pars$betas_g1)
  if (str_detect(variant, "mixture")) {
    betas$g2 <- c(pars$betas_g1[1] + exp(pars$betas_diff[1]),
                  pars$betas_other_g2[1:5], 
                  pars$betas_g1[7] - exp(pars$betas_diff[2]),
                  pars$betas_other_g2[6:10])
  }
  betas <- betas %>% 
    dplyr::mutate(month = 1:12) %>% 
    gather(grp, value, - month)
  
  p_betas <- ggplot(betas, aes(x = month, y = value, color = factor(grp))) + 
    geom_line() + 
    theme_bw() + 
    ylab("betas")
  
  if (!is.null(pars$ind_betas_diff)) {
    p_betas <- p_betas + 
      geom_vline(data = tibble(month = pars$ind_betas_diff,
                               grp = c("g2", "g1")),
                 aes(xintercept = month, color = factor(grp)),
                 lty = 2)
  }
  
  etas <- tibble(value = pars$etas_tilde[, 1] * pars$sigma_etas, year = data_mapping$u_years)
  p_etas <- ggplot(etas, aes(x = year, y = value)) + geom_bar(stat = "identity") + theme_bw() + ylab ("etas")
  
  cntry.sf <- cntry.sf[order(cntry.sf$id), ]
  
  if (str_detect(variant, "mixture")) {
    
    cntry.sf$theta_lambdas <- pars$theta_lambdas_raw
    cntry.sf$phi_lambdas <- pars$lambdas_raw
    cntry.sf$lambda <- 1/sqrt(pars$tau_lambdas[1]) * (sqrt(1-pars$rho_lambdas[1]) * pars$theta_lambdas_raw + sqrt(pars$rho_lambdas[1]/scaling_factor) * pars$lambdas_raw)
    cntry.sf$lambda  <- 1/(1 + exp(-cntry.sf$lambda))
  }
  
  cntry.sf$theta <- pars$theta
  cntry.sf$phi <- pars$phi
  cntry.sf$conv_re <- 1/sqrt(pars$tau) * (sqrt(1-pars$rho) * pars$theta + sqrt(pars$rho/scaling_factor) * pars$phi)
  
  p_spat <- cntry.sf %>% select(theta, phi, conv_re) %>% gather(var, value, -geometry) %>%
    ggplot(aes(fill = value)) +
    geom_sf() + 
    facet_wrap(~var) + 
    scale_fill_gradient2() + 
    theme_bw()
  
  if (str_detect(variant, "mixture")) {
    
    p_grp <- cntry.sf %>% select(theta_lambdas, phi_lambdas) %>% gather(var, value, -geometry) %>%
      ggplot(aes(fill = value)) +
      geom_sf() + 
      facet_wrap(~var) + 
      scale_fill_gradient2(low = "#1B5BA8", high = "#F5B540") +
      theme_bw()
    
    p_lambda <- cntry.sf %>% select(lambda) %>% gather(var, value, -geometry) %>%
      ggplot(aes(fill = value)) +
      geom_sf() + 
      facet_wrap(~var) +
      scale_fill_gradient2(midpoint = .5, low = "#1B5BA8", high = "#F5B540") +
      theme_bw()
    
    p <- arrangeGrob(grobs = list(p_betas, p_etas, p_spat, p_grp, p_lambda), layout_matrix = rbind(c(1, 1, 2, 2), c(3, 3, 3, 3), c(4, 4, 5, 5)))
    
  } else {
    p <- arrangeGrob(grobs = list(p_betas, p_etas, p_spat), layout_matrix = rbind(c(1, 1, 2, 2), c(3, 3, 3, 3)))
  }
  
  ggsave(p, filename = plot_file, width = 10, height = 15, type = "cairo")
  
}

# Workaround for MARCC since lwgeom is cannot be installed
st_dist <- function(x, y) {
  cx <- st_coordinates(x)
  cy <- st_coordinates(y)
  distmat <- apply(cx, 1, function(j) apply(cy, 1, function(i)
    sqrt((i[2] - j[1])^2 + (i[1] - j[1])^2)))
  distmat
}

#' Set do map version 2
#' @description Set flag to run version 2 of threaded stan code
#' 
#' @param config 
#' 
#' @details This dropping multi-year observations
#' @return logical
#' 
setDoMapV2 <- function(config) {
  do_map_v2 <- FALSE 
  
  if (config$country_iso3 %in% c("NGA", "COD")) {
    if (!is.null(config$run_level) & !is.na(config$run_level)) {
      if (config$run_level == 2) {
        config$gadm_levels <- "12"
        do_map_v2 <- T
      }
    }
  }
  
  do_map_v2
}
