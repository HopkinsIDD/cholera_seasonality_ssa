# This script postprocesses the outputs of a given set of model runs
# 
# Preamble ---------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(sf)
library(spdep)
library(lubridate)
library(rstan)
library(optparse)
library(foreach)

# User-supplied options
# Parse options
opt_list <- list(
  make_option(c("-d", "--cholera_dir"), default = "./",
              action ="store", type = "character", help = "Directory"),
  make_option(opt_str = c("-r", "--redo"), type = "logical",
              default = T, help = "redo computations"),
  make_option(opt_str = c("-s", "--redo_single"), type = "logical",
              default = T, help = "redo computations"),
  make_option(opt_str = c("-i", "--identifier"), type = "character",
              default = "mean_annual_incidence", help = "identifier")
)

opt <- parse_args(OptionParser(option_list = opt_list))

source(str_c(opt$cholera_dir, "/analysis/utils.R"))
source(str_c(opt$cholera_dir, "/analysis/modeling_utils.R"))
source(str_c(opt$cholera_dir, "/analysis/postprocessing_utils.R"))
source(str_c(opt$cholera_dir, "/analysis/plotting_utils.R"))

# Model selection ---------------------------------------------------------

Sys.setenv("GEOPKG" = "TRUE")

# Get best models
best_models <- getAllBestModels(countries = "all",
                                run_levels = "best",
                                identifier = opt$identifier,
                                redo = opt$redo,
                                redo_single = opt$redo_single)

# Save best models
saveRDS(best_models, file = str_glue("generated_data/model_outputs/best_models_{opt$identifier}.rds"))

# Seasonality estimates ---------------------------------------------------

# Best simulated seasonality coefficients
best_seas <- getBestSeas(countries = "all",
                         run_level = "best",
                         identifier = opt$identifier,
                         redo = list(main = opt$redo,
                                     best_seas = opt$redo,
                                     best_seas_single = opt$redo,
                                     data = F,
                                     best_model = F))

saveRDS(best_seas, file = str_glue("generated_data/model_outputs/best_seas_{opt$identifier}.rds"))

best_seas_sf <- getSfObjects(best_seas %>% 
                               replace_na(list(mean = 0)))

saveRDS(best_seas_sf, file ="generated_data/clustering_sf_objects.rds")

# Get seasonality estimates by group
best_seas_grp <- getBestSeasByGroup(countries = "all",
                                    identifier = opt$identifier,
                                    run_level = "best",
                                    redo  = list(main = opt$redo,
                                                 best_seas = opt$redo,
                                                 best_seas_single = opt$redo)) 

saveRDS(best_seas_grp, file = str_glue("generated_data/model_outputs/best_seas_by_grp_{opt$identifier}.rds"))


# Seasonality index -------------------------------------------------------

# Get seasonality estimates by group
seas_index <- computeAllSeasIndex(countries = "all",
                                  identifier = opt$identifier,
                                  run_level = "best",
                                  redo = list(main = opt$redo,
                                              probs = T)) 

saveRDS(seas_index, file = str_glue("generated_data/model_outputs/seas_index_{opt$identifier}.rds"))

