# This script extracts model fit results for plotting and analysis
# 
# Preamble ---------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(sf)
library(spdep)
library(lubridate)
library(rstan)
library(optparse)

# User-supplied options
option_list <- list(
  make_option(c("-c", "--config"), default = "./analysis/02_modeling/configs/default_config.yml", 
              action ="store", type = "character", help = "Configuration file"),
  make_option(c("-d", "--cholera_dir"), default = "./",
              action ="store", type = "character", help = "Directory")
)

# Parse options
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

source(str_c(opt$cholera_dir, "/analysis/utils.R"))
source(str_c(opt$cholera_dir, "/analysis/modeling_utils.R"))
source(str_c(opt$cholera_dir, "/analysis/postprocessing_utils.R"))

# Load model configuration file
config <- yaml::read_yaml(opt$config)

# Extract -----------------------------------------------------------------

if (is.na(as.numeric(config$run_level))) {
  stop("Please specify run level")
}

fit_file <- makeResultsFile(opt = opt, 
                            cntry_iso3 = config$country_iso3, 
                            run_level = as.numeric(config$run_level), 
                            drop_multi_yr = setDoMapV2(config))

output_file <- makeOutputFile(fit_file)

if (!file.exists(output_file) | file.info(output_file)$ctime < file.info(fit_file)$ctime) {
  
  print(fit_file)
  
  # Extract run info
  run_info <- parseResFilename(fit_file)
  
  cntry <- run_info$country
  run_level <- run_info$run_level
  model <- run_info$variant
  
  # Load fit
  res <- readRDS(fit_file)
  
  if (length(res$chol_stanfit@stan_args) == 0) {
    print("-- Does not contain samples")
    return(NULL)
  }
  # Extract parameter traces to assess convergence
  par_conv_df <- extractPars(model_fit = res$chol_stanfit, 
                             model_type = model, 
                             parset = "convergence")
  print("Done conv par")
  # Extract parameter posterior statistics
  par_mod_df <- extractPars(model_fit = res$chol_stanfit, 
                            model_type = model, 
                            parset = "modeling",
                            trace = F) 
  print("Done mod par")
  
  # Compute offset probability
  offset_probs <- extractOffset(model_type = model,
                                chol_stanfit = res$chol_stanfit)
  print("Done offsets")
  
  # Model performance
  # If available load post-fitting simulations of loglik
  cmdstanfit_file <- makeCmdstanrFile(fit_file)
  # Check if output file for simulations already exists
  llsim_file <- makeCmdstanrSimFile(cmdstanfit_file) %>% 
    str_replace("cmdstanr", "cmdstanr_sim") %>% 
    str_replace("sim_", "sim_loglik_")
  
  if (file.exists(llsim_file)) {
    llsim <- readRDS(llsim_file)
    ll_mat <- llsim$draws(variables = "log_lik")
    ll_mat <- matrix(ll_mat, prod(dim(ll_mat)[1:2]), dim(ll_mat)[3])
    rm(llsim)
    gc()
  } else {
    ll_mat <- loo::extract_log_lik(res$chol_stanfit)
  }
  
  # Compute LOO and WAIC
  nrow_old <- nrow(ll_mat)
  # Remove samples if there are NAs
  ll_mat[is.infinite(ll_mat)] <- NA
  ll_mat <- na.omit(ll_mat)
  
  # Remove all samples which haven't completed running
  ll_mat <- ll_mat[ll_mat[,2]!=0, ]
  
  cat("!! Removing", nrow_old - nrow(ll_mat), "samples for which the ll is NA \n")
  
  model_loo <- loo::loo(ll_mat)
  print("Done loo")
  
  model_waic <- loo::waic(ll_mat)
  print("Done WAIC")
  
  output <- list(
    fit_file = fit_file,
    cntry = cntry,
    run_level = run_level,
    model = model,
    time_left = run_info$time_left,
    time_right = run_info$time_right,
    par_conv_df = par_conv_df,
    par_mod_df = par_mod_df,
    offset_probs = offset_probs,
    model_IC = list(waic = model_waic,
                    loo = model_loo),
    data = res$data,
    sf_object = res$sf_object,
    timings = rstan::get_elapsed_time(res$chol_stanfit),
    stan_stats = list(n_max_tree_depth = rstan::get_num_max_treedepth(res$chol_stanfit),
                      n_diverge = rstan::get_num_divergent(res$chol_stanfit))
  )
  saveRDS(output, file = output_file)
}
