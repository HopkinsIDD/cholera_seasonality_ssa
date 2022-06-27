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

fit_file <- makePriorFile(opt = opt, 
                            cntry_iso3 = config$country_iso3, 
                            run_level = as.numeric(config$run_level), 
                            drop_multi_yr = setDoMapV2(config))

output_file <- makePriorOutputFile(fit_file)

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
  
  output <- list(
    fit_file = fit_file,
    cntry = cntry,
    run_level = run_level,
    model = model,
    time_left = run_info$time_left,
    time_right = run_info$time_right,
    par_conv_df = par_conv_df,
    par_mod_df = par_mod_df,
    timings = rstan::get_elapsed_time(res$chol_stanfit),
    stan_stats = list(n_max_tree_depth = rstan::get_num_max_treedepth(res$chol_stanfit),
                      n_diverge = rstan::get_num_divergent(res$chol_stanfit))
  )
  saveRDS(output, file = output_file)
}
