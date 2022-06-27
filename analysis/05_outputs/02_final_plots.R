# This script produces the final plots
# 
# Preamble ----------------------------------------------------------------
library(tidyverse)
library(sf)
library(lubridate)
library(foreach)
library(itertools)
library(tidync)
library(optparse)
library(kableExtra)

# Source helper functions
source("analysis/utils.R")
source("analysis/postprocessing_utils.R")
source("analysis/plotting_utils.R")
source("analysis/simulation_utils.R")

# Parse options
opt_list <- list(
  optparse::make_option(opt_str = c("-r", "--redo"), type = "logical",
                        default = F, help = "redo computations"),
  optparse::make_option(opt_str = c("-s", "--redo_single"), type = "logical",
                        default = F, help = "redo computations"),
  optparse::make_option(opt_str = c("-i", "--identifier"), type = "character",
                        default = "mean_annual_incidence", help = "identifier"),
  optparse::make_option(opt_str = c("-d", "--data_path"), type = "character",
                        default = "generated_data/monthly_data_20220215.csv", help = "identifier")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

Sys.setenv("GEOPKG" = "TRUE")

# Get estimates -----------------------------------------------------------

# Get seasonality estimates
best_seas <- getBestSeas(countries = "all",
                         identifier = opt$identifier,
                         run_level = "best",
                         redo = list(main = F)) 

sf_objects <- getSfObjects(best_seas)

# Get names for export
sf_objects <- sf_objects %>% 
  group_by(gadm_lev) %>% 
  group_modify(function(x, y) {
    inner_join(x, getGadmName(x$gid, y$gadm_lev)) 
  }) %>% 
  ungroup() %>% 
  st_as_sf()

# Save for sharing
best_seas <- best_seas %>% 
  inner_join(sf_objects %>% 
               ungroup() %>% 
               st_drop_geometry() %>% 
               select(gid, name))

best_seas %>% 
  mutate(month = str_extract(variable, "(?<=\\,)[0-9]+")) %>% 
  select(country, gid, name, run_level, variant, month, mean, q2.5, q97.5) %>% 
  rename(q025 = q2.5,
         q975 = q97.5,
         model = variant) %>% 
  arrange(country, gid) %>% 
  write_csv("generated_data/model_outputs/best_seasonality_estimates.csv")

# Get seasonality estimates
best_seas_grp <- getBestSeasByGroup(countries = "all",
                                    identifier = opt$identifier,
                                    run_level = "best",
                                    redo = list(main = F)) 

# Figure 1 ----------------------------------------------------------------

# Plot of seasonality coefficients as a map
p_maps <- plotAllSeasMaps(identifier= opt$identifier,
                          run_level = "best",
                          best_seas = best_seas,
                          sf_objects = NULL, 
                          redo = opt$redo,
                          save_fig = T,
                          verbose = T)

# Plot the seasonality coefficient
p_seasindex <- plotAllSeasInd2(identifier = opt$identifier,
                               run_level = "best",
                               redo = opt$redo,
                               save_fig = T)

# Plot combined plot
p_seas_combined <- plotAllSeasCombined2(p_seas_index = p_seasindex,
                                        p_maps = p_maps,
                                        save_fig = T)

ggsave(p_seas_combined, filename = "figures/pub_figures/Figure_1.svg", width = 9, height = 6)

# Figure 2: Seas groups ----------------------------------------------------------------

p_seas_groups <- plotAllSeasGroups2(countries = "all",
                                    identifier= opt$identifier,
                                    run_level = "best",
                                    best_seas = best_seas,
                                    sf_objects = sf_objects, 
                                    n_groups = 5,
                                    redo = opt$redo,
                                    save_fig = T,
                                    verbose = F)

ggsave(p_seas_groups, filename = "figures/pub_figures/Figure_2.svg", width = 7.5, height = 6)

# Figure 3: covariate correlation ----------------------------------------
p_corr <- plotAllCovarCorr(identifier = opt$identifier,
                           run_level = "best", 
                           covariates = c("Precipitation", "FloodMean", "Mean_T"),
                           best_seas = best_seas,
                           sf_objects = NULL, 
                           method = "spearman",
                           lags = 0:2, 
                           redo = opt$redo,
                           save_fig = T)

ggsave(p_corr, filename = "figures/pub_figures/Figure_3.svg", width = 7.5, height = 6)

p_corr <- plotAllCovarCorr(identifier = opt$identifier,
                           run_level = "best", 
                           covariates = "all",
                           best_seas = best_seas,
                           sf_objects = NULL, 
                           method = "spearman",
                           lags = 0, 
                           redo = opt$redo,
                           save_fig = T)

p_corr2 <- p_corr + facet_wrap(~covar, labeller = labeller(covar = getCovarDict()))

ggsave(p_corr2, 
       filename = "figures/pub_figures/corr_all_single.png", width = 8, height = 5,
       dpi = 400)

# Figure S1: data ---------------------------------------------------------

# Plot of seasonality coefficients
p_data_maps <- plotAllDataMaps(countries = "all",
                               identifier= opt$identifier,
                               run_level = "best",
                               redo = opt$redo,
                               redo_single = opt$redo_single,
                               redo_data = F,
                               save_fig = T,
                               verbose = F)

# Figure S2: Prior shrinkage ---------------------------------------------
p_shrink <- plotAllPriorShrinkage(countries = "all",
                                  identifier = opt$identifier,
                                  redo = opt$redo,
                                  redo_single = opt$redo_single,
                                  save_fig = T)

# Figure S3: Retoridctive check coverage ---------------------------------
p_coverage <- plotAllCoverage(countries = "all",
                              identifier = opt$identifier,
                              redo = opt$redo,
                              save_fig = T,
                              verbose = T) 

# Figure S4: Seasonality time series  --------------------------------------------------------------
p_seas_ts <- plotAllSeasTS(identifier = opt$identifier,
                           run_level = "best",
                           best_seas_grp = best_seas_grp,
                           sf_objects = sf_objects, 
                           redo = opt$redo,
                           save_fig = T)

# Figure S5: Uncertainty in seasonality ----------------------------------
p_seas_uncertainty <- plotAllSeasMapsUncertainty(identifier= opt$identifier,
                                                 run_level = "best",
                                                 best_seas = best_seas,
                                                 sf_objects = NULL, 
                                                 redo = opt$redo,
                                                 save_fig = T,
                                                 verbose = T)

# Figure S6: Seas index uncertainty --------------------------------------------------------------
p_seas_index_uncertainty <- plotAllSeasIndexUncertainty(identifier = opt$identifier,
                                                        run_level = "best",
                                                        redo = opt$redo,
                                                        save_fig = T)

# Figure S7: Seas index scatter --------------------------------------------------------------
p_seas_index_scatter <- plotAllSeasIndexScatter(identifier = opt$identifier,
                                                run_level = "best",
                                                redo = opt$redo,
                                                save_fig = T)

# Figure S8: Peak month --------------------------------------------------------------
p_peak_month <- plotAllPeakMonth(identifier = opt$identifier,
                                 run_level = "best",
                                 best_seas = best_seas,
                                 sf_objects = sf_objects, 
                                 redo = opt$redo,
                                 save_fig = T)

# Figure S9: Offset ---------------------------------------------------------
p_offset <- plotAllOffsets(identifier = opt$identifier,
                           run_level = "best",
                           redo = opt$redo,
                           save_fig = T)

# Figure S10: Model comparison clustering ---------------------------------
p_clustercomp <- plotClusteringModelComp(identifier = opt$identifier)

# Figure S11: Alternative groupings ---------------------------------------
# Plot of seasonality coefficients
p_seas_other_groups <- map(2:7,
                           ~plotAllSeasGroups(countries = "all",
                                              identifier= opt$identifier,
                                              run_level = "best",
                                              best_seas = best_seas,
                                              sf_objects = sf_objects, 
                                              n_groups = as.character(.),
                                              redo = opt$redo,
                                              save_fig = T,
                                              verbose = F,
                                              label = NULL))

cowplot::plot_grid(plotlist = p_seas_other_groups, 
                   ncol = 2,
                   labels = "auto") %>% 
  ggsave(filename = str_c(getFigDir(), "/seas_group_map_all_k2-3-4-5-6-7_rlbest_mean_annual_incidence_all-all.png"),
         width = 12,
         height = 12)

# Figure S12: Correlation statistics  -----------------------------
p_corrstats <- plotCorrStats(countries = "all",
                             identifier = opt$identifier,
                             run_level = "best", 
                             covariates = c("Precipitation", "FloodMean", "Mean_T"),
                             redo_single = opt$redo,
                             method = "spearman",
                             save_fig = T,
                             verbose = T)

# Figure S13: Correlation with other covariates ---------------------------
p_pcorr <- plotAllCovarCorr(identifier = opt$identifier,
                            run_level = "best", 
                            covariates = "all",
                            best_seas = best_seas,
                            sf_objects = sf_objects, 
                            lags = 0:3, 
                            redo = opt$redo,
                            method = "pearson",
                            save_fig = T)

p_scorr <- plotAllCovarCorr(countries = "all",
                            identifier = opt$identifier,
                            run_level = "best", 
                            covariates = "all",
                            best_seas = best_seas,
                            sf_objects = sf_objects, 
                            lags = 0:3, 
                            redo = opt$redo,
                            method = "spearman",
                            save_fig = T)




# Figure S14:  Other models ----------------------------------------------------
p_comp <- plotDefinitionComp(redo = opt$redo,
                             redo_single = opt$redo_single)


# Table S1: Model comparison ---------------------------------------------
all_loo <- runAll(countries = "all",
                  run_levels = "best",
                  identifiers = opt$identifier,
                  fun = getModelIC,
                  fun_name = "model_IC",
                  redo = opt$redo)

country_dict <- read_csv("data/all_country_geonaming.csv",
                         col_types = cols()) %>% 
  select(name, alpha_3) %>% 
  {
    vec <- .$name
    names(vec) <- .$alpha_3
    vec
  }

country_dict["COD"] <- "Democratic Republic of the Congo"
country_dict["TCD"] <- "Tchad"
country_dict["TZA"] <- "Tanzania"

# Make table
all_loo_for_table <- all_loo %>% 
  select(country, variant, elpd_diff, se_diff, looic) %>% 
  group_by(country) %>% 
  mutate(rank = rank(-elpd_diff),
         deltaloo = looic - min(looic),
         weight = exp(-.5 * deltaloo)/sum(exp(-.5 * deltaloo))) %>% 
  ungroup() %>% 
  mutate(country = country_dict[country]) %>% 
  arrange(country, rank)

all_loo_for_table %>% 
  select(-country, -deltaloo) %>% 
  mutate(weight = formatC(weight, digits = 2, format = "f")) %>% 
  mutate_if(is.numeric, formatC, digits = 1, format = "f") %>% 
  mutate(rank = as.integer(rank)) %>% 
  mutate_if(is.numeric, as.character) %>% 
  mutate(elpd_diff = ifelse(elpd_diff == "0.0", "-", elpd_diff),
         se_diff = ifelse(se_diff == "0.0", "-", se_diff)) %>% 
  kable(col.names = c("model", "ELPD diff", "SE", "LOO-IC", "rank", "weight")) %>%
  kable_classic(bootstrap_options = "striped", 
                full_width = F, 
                position = "left") %>% 
  pack_rows(index = table(all_loo_for_table$country)) %>% 
  row_spec(0, bold = T)



