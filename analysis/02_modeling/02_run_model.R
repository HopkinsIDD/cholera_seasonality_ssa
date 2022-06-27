# This script runs the seasonality model
# 
# Preamble ----------------------------------------------------------------

library(rstan)   
library(cmdstanr)
library(INLA)
library(tidyverse)
library(magrittr)
library(sf)
library(spdep)
library(foreach)
library(itertools)
library(lubridate)
library(optparse)

# User-supplied options
option_list <- list(
  make_option(c("-c", "--config"), default = "./analysis/02_modeling/configs/default_config.yml", 
              action ="store", type = "character", help = "Configuration file"),
  make_option(c("-d", "--cholera_dir"), default = "./",
              action ="store", type = "character", help = "Directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load model configuration file
config <- yaml::read_yaml(opt$config)

# Arguments
variant <- str_split(config$variant, "_")[[1]]
config$run_level <- as.character(config$run_level)

# Source helper functions
source(str_c(opt$cholera_dir, "analysis/modeling_utils.R"))
source(str_c(opt$cholera_dir, "analysis/utils.R"))
source(str_c(opt$cholera_dir, "analysis/nbgraph_utils.R"))

set.seed(1234)

# Setup model -------------------------------------------------------------

do_map_v2 <- setDoMapV2(config)

# Get GADM shapefile
cntry.sf <- getGadmSf(country = config$country, 
                      gadm_lev = 2)

# Choice of admin level --------------------------------------------------------
run_data_file <- makeRunDataFile(opt = opt,
                                 config = config)

if (file.exists(run_data_file) & !config$redo_data) {
  run_data <- readRDS(run_data_file)
  cat("-- Loading pre-computed run data \n")
} else {
  # Get the data to run
  run_data <- getRunData(
    path_to_monthly_data = str_c(opt$cholera_dir, config$monthly_data), 
    country_iso3 = config$country_iso3,
    gadm_levels = config$gadm_levels,
    keep_national_yearly = config$keep,
    start_date = config$start_date,
    end_date = config$end_date,
    thresh = config$thresh,
    case_thresh = config$case_thresh
  )
  saveRDS(run_data, run_data_file)
}

# Get the GADM level on which to run the model, this depends on spatial data availability
if (is.null(config$run_level) | is.na(config$run_level) | config$run_level == "best") {
  run_level <- runLevel(country = config$country_iso3,
                        run_data = run_data)
} else{
  run_level <- as.numeric(config$run_level)
}

cat("############################## \n\n\n", 
    "Running:\t",  config$country_iso3,
    "\n period:\t", str_c(ifelse(is.null(config$start_date), "-Inf", config$start_date), ifelse(is.null(config$end_date), "Inf", config$end_date), sep = " - "),
    "\n at GADM level:\t", run_level ,
    "\n with cholera:\t", config$thresh,
    "\n with model:\t", config$variant,
    "\n data file:\t", config$monthly_data,
    "\n\n\n##############################\n")

# Check if results have already been computed
res_file <- makeResultsFile(opt, config$country_iso3, run_level, do_map_v2)

if (file.exists(res_file)) {
  # Load result file to check if it contains samples
  chol_stanfit <- readRDS(res_file)
  chol_stanfit <- chol_stanfit$chol_stanfit
  
  if (length(chol_stanfit@stan_args) == 0) {
    cat("-- Previous results file does not contain samples, re-running model \n")
    rm(chol_stanfit)
  } else if (config$redo_comp) {
    cat("-- Redoing data modeling \n")
  } else {
    stop("---- Data already modeled")
  }
}

# Data preparation --------------------------------------------------------

# Compute adjacency
adjacency_file <- makeAdjacencyFile(opt, config, run_level = run_level)

if (file.exists(adjacency_file) & !config$redo_dat) {
  adjacency <- readRDS(adjacency_file)
  cat("-- Loading previously computed adjacency \n")
} else {
  adjacency <- computeAdjacency(country = config$country_iso3, 
                                run_level = run_level)
  
  saveRDS(adjacency, file = adjacency_file)
}

# Comput prior of spatial autocorrelatoin
prior_rho_file <- makePriorRhoFile(opt, config, run_level = run_level)

if (file.exists(prior_rho_file) & !config$redo_data) {
  prior_rho <- readRDS(prior_rho_file)
  cat("-- Loading previously computed prior rho \n")
} else {
  # Compute prior of spatial variance contribution
  prior_rho <- computeICARPriors(N = adjacency$adj_list$N,
                                 node1 = adjacency$adj_list$node1,
                                 node2 = adjacency$adj_list$node2,
                                 U_rho = .5,
                                 alpha_rho = 2/3)
  
  saveRDS(prior_rho, file = prior_rho_file)
}

# Admin unit and time mappings ---------------------------------------
# Save mapping file
mapping_file <- makeMappingFile(opt, config, run_level)

if (file.exists(mapping_file) & !config$redo_data) {
  data_mapping <- readRDS(mapping_file)
  cat("-- Loading previously computed data mapping \n")
} else {
  # Observations to admin unit id mapping
  data_mapping <- obsToMappingFlat(country = config$country_iso3,
                                   run_data = run_data, 
                                   run_level = run_level, 
                                   variant = variant,
                                   drop_multi_yr = do_map_v2,
                                   version = "speedup")
  saveRDS(data_mapping, file = mapping_file)
}

rm("run_data")

# Data instanciation -------------------------------------------------------
data <- instantiateStanData(adjacency = adjacency,
                            data_mapping = data_mapping,
                            prior_rho = prior_rho,
                            # Parameters for prior of spatial marginal variance
                            U_tau  = 1,
                            alpha_tau = .75)

# Set variant-specific parameter
data <- append(
  data, 
  list(N_offsets = ifelse(str_detect(config$variant, "offset"), 1, 0),
       N_groups =  ifelse(str_detect(config$variant, "mixture"), 2, 1),
       N_eta_groups =  ifelse(str_detect(config$variant, "mixture"), 
                              ifelse(str_detect(config$variant, "mixture") &
                                       str_detect(config$variant, "offset"), 1, 1), 1))
)


# Parameter initialization ------------------------------------------------

# If available get initial values from base model, if not generate values
# that give reasonable probabilities
init_pars <- initPars(n_chain = config$stan$n_chains, 
                      data_mapping = data_mapping,
                      data = data,
                      cntry.sf = adjacency$sf_object,
                      variant = config$variant,
                      scaling_factor = prior_rho$scaling_factor,
                      cntry_iso3 = config$country_iso3,
                      run_level = run_level,
                      Q_inv_full = prior_rho$Q_inv_full,
                      do_plot = F,
                      plot_file = NULL,
                      ind_betas_diff = config$ind_betas_diff)

# For mixture-type models extract the indices for seasonality coefficient ordering
# between groups
# if (any(str_detect(config$variant, "mixture")) | is.null(init_pars[[1]]$ind_betas_diff)) {
data$ind_betas_diff <- config$ind_betas_diff
data$beta_months <- setdiff(1:12, data$ind_betas_diff) %>% sort()
# } else {
#   data$beta_months <- rep(1, 10)
# }

init_pars <- map(init_pars, ~ .[-which(names(.) == "ind_betas_diff")])

# Run stan ----------------------------------------------------------------

# Select model based on config and run level
chol_stan <- getStanModel(config, run_level)

interm_fbase <- makeIntermediateFileBase(res_file)
# Write data to json
json_data_file <-  makeJsonFile(opt, config)
cmdstanr::write_stan_json(data, json_data_file)

if (file.exists(str_glue("{interm_fbase}-1.csv")) & !config$redo_comp) {
  
  # Transform back to stanfit object
  cmdstan_fit <- cmdstanr::read_cmdstan_csv(dir(str_c(opt$cholera_dir, "/interm/cmdstanr_int/"),
                                                full.names = T,
                                                pattern = str_remove(interm_fbase, "interm/cmdstanr_int/"), 
                                                recursive = T) %>% 
                                              str_subset("csv"))
  
  cmdstan_fit$save_object(file = makeCmdstanrFile(res_file))
  
  # Transform back to stanfit object
  chol_stanfit <- rstan::read_stan_csv(dir(str_c(opt$cholera_dir, "/interm/cmdstanr_int/"),
                                           full.names = T,
                                           pattern = str_remove(interm_fbase, "interm/cmdstanr_int/"), 
                                           recursive = T) %>% 
                                         str_subset("csv"))
  
} else {
  cmdstan_fit <- chol_stan$sample(data = json_data_file,
                                  seed = 1234,
                                  init = init_pars,
                                  sig_figs = 5,
                                  chain_ids = seq_len(config$stan$n_chains),
                                  chains = config$stan$n_chains,
                                  parallel_chains = config$stan$n_chains,
                                  threads_per_chain = getNumThreads(config),
                                  iter_warmup = config$stan$n_warmup,
                                  iter_sampling = config$stan$n_iter,
                                  max_treedepth = config$stan$control$max_treedepth,
                                  metric = config$stan$control$metric,
                                  adapt_delta = config$stan$control$adapt_delta,
                                  save_warmup = F,
                                  refresh = 50,
                                  save_latent_dynamics = F,
                                  output_dir = opt$cholera_dir,
                                  output_basename = interm_fbase)
  
  
  cmdstan_fit$save_object(file = makeCmdstanrFile(res_file))
  
  # Transform back to stanfit object
  chol_stanfit <- rstan::read_stan_csv(cmdstan_fit$output_files())
}

# Save output
saveRDS(list(config = config, 
             run_level = run_level,
             sf_object = adjacency$cntry.sf,
             data = data,
             chol_stanfit = chol_stanfit,
             init_pars = init_pars),
        file = res_file)


