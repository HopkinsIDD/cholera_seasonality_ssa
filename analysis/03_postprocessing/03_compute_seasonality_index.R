# This script computes the seasonality index for a given config
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

if (config$variant == "null") {
  stop("No seasonality index for null model")
}

# Extract -----------------------------------------------------------------

if (is.na(as.numeric(config$run_level))) {
  stop("Please specify run level")
}

fit_file <- makeResultsFile(opt = opt, 
                            cntry_iso3 = config$country_iso3, 
                            run_level = as.numeric(config$run_level), 
                            drop_multi_yr = setDoMapV2(config))

output_file <- makeOutputSimFile(fit_file) %>% 
  str_replace("\\.rds", "_seas_index.rds")

# Data mapping with additional information on unique units
data_mapping <- readRDS(makeMappingFile(opt = opt, 
                                        config = config, 
                                        run_level = as.numeric(config$run_level)))

# Get admin-level seasonality coefficients
betas <- getBetas(country = config$country_iso3,
                  run_level = config$run_level,
                  model = config$variant,
                  identifier = config$thresh,
                  time_left = config$start_date,
                  time_right = config$end_date,
                  gadm_lev = config$gadm_levels,
                  redo = getRedoDefaults(),
                  verbose = F)

# Get peak month
peak_month <- betas %>%
  dplyr::group_by(country, gid, variant, run_level) %>% 
  dplyr::summarise(peak_month = which.max(mean),
                   peak_val = max(mean)) 

# Get admin-level seasonality coefficients
probs <- getProbSamples(country = config$country_iso3,
                        run_level = config$run_level,
                        model = config$variant,
                        identifier = config$thresh,
                        time_left = config$start_date,
                        time_right = config$end_date,
                        gadm_lev = config$gadm_levels,
                        redo = getRedoDefaults(),
                        verbose = F)

seas_index <- map_df(1:length(unique(betas$gid)), function(x) {
  inds <- which(data_mapping$uunits == x)
  ps <- probs[,,inds]
  ps <- matrix(ps, prod(dim(ps)[1:2]), dim(ps)[3])
  pm <- peak_month$peak_month[peak_month$gid == data_mapping$u_unit[x]]
  
  map_df(1:nrow(ps), function(y) {
    tibble(
      seas_index = calcSeasIndex2(
      probs = ps[y, ], 
      months = data_mapping$umonths[inds], 
      peak_month = pm),
      draw = y)
  }) %>% 
    mutate(gid = data_mapping$u_unit[x])
})

saveRDS(seas_index, file = output_file)