# This script produces the final statistics for
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

# Sink messages to log file
log_file <- str_glue("generated_data/model_outputs/stats_log_{opt$identifier}.txt")

sink(log_file, type = "output", append = F)

cat("################################################################\n",
    "# Data and outputs statistics for model runs \n",
    "# Occurrence definition:", opt$identifier, "\n",
    "# Produced at", as.character(now()),
    "\n################################################################\n\n\n")

# Statement: data ---------------------------------------------------------

cat("# Data \n")

# Load all data used in models that ran
all_data <- vroom::vroom(opt$data_path, 
                            col_types = cols(),
                         progress = F) %>% 
  mutate(type = case_when(period_type %in% c("day", "week", "biweek", "Oweek", "month") ~ "sub-monthly",
                          period_type == "year" ~ "year", 
                          T ~ "other")) %>% 
  mutate(cholera = pmax(suspected_cases, confirmed_cases)) %>% 
  # Keep valid observations
  dplyr::filter(country %in% getAllCountries() &
                  !(country %in% getExcludedCountries()) &
                  ((in_upper & in_lower) | (!in_lower & cholera > 0) | (!in_upper & cholera == 0)))

# num obs
cat("Number of observations:", formatC(nrow(all_data), big.mark = "'"), "\n\n")

cat("Observations by gadm level: \n")
# num obs per admin level
all_data %>% 
  count(gadm_lev) %>% 
  pwalk(~ cat("GADM lev", ..1, ":" , formatC(..2, big.mark = "'"), "observations \n"))

cat("\n")

# date range
cat("Observations span", min(all_data$year_left), "-", max(all_data$year_left), "\n\n")

cat("Of which", formatC(sum(all_data$year_left >= 2000)/nrow(all_data)*100,
                        digits = 3), "% are after the year 2000 \n\n")


# gadm units
distinct_admin_units <- all_data %>% 
  distinct(gid, gadm_lev)

cat("Units by gadm level: \n")
distinct_admin_units %>% 
  count(gadm_lev) %>% 
  pwalk(~ cat("GADM lev", ..1, ":" , formatC(..2, big.mark = "'"), "admin units \n"))

# % monthly/sub-monthly
obs_by_type <- all_data %>% 
  count(type)

n_sub_monthly <- obs_by_type$n[obs_by_type$type == "sub-monthly"]
n_yearly <- obs_by_type$n[obs_by_type$type == "year"]

cat("There are", formatC(n_sub_monthly, big.mark = "'"), "(",
    formatC(n_sub_monthly/sum(obs_by_type$n)*100, digits = 3), 
    "%) sub-monthly observations \n")

cat("There are", formatC(n_yearly, big.mark = "'"), "(",
    formatC(n_yearly/sum(obs_by_type$n)*100, digits = 3), "% of total,",
    formatC(n_yearly/(sum(obs_by_type$n) - n_sub_monthly)*100, digits = 3),
    "% of > monthly) yearly observations \n\n")

# Statement: Model selection ----------------------------------------------

cat("# Models \n")

# Get best models
best_models <- getAllBestModels(countries = "all",
                                run_levels = "best",
                                identifier = opt$identifier,
                                redo = F,
                                redo_single = F)

cat("No results are available for", length(setdiff(getAllSSACountries(), best_models$country)), " countries:",
    setdiff(getAllSSACountries(), best_models$country),"\n")
# cat("Excluded", nrow(getDropCountries()), "countries because of lack of data:",
#     getDropCountries()$country,"\n")
cat("There were", sum(best_models$best_model != "null"), "/", nrow(best_models), "countries with seasonality \n")
cat("Countries without seasonality:", best_models$country[best_models$best_model == "null"] %>% sort() %>% str_c(collapse = ", "), "\n")
cat("There were", sum(str_detect(best_models$best_model, "mixture")), "/", nrow(best_models), "countries with mixture \n")

# Population count
all_ssa_pop <- getGadmData() %>% 
  dplyr::filter(country %in% getAllSSACountries()) %>% 
  dplyr::mutate(seasonal = country %in% best_models$country[best_models$best_model != "null"]) %>% 
  dplyr::summarise(frac_seasonal = sum(popsum_2020[seasonal])/sum(popsum_2020))

cat("A total of", all_ssa_pop$frac_seasonal %>% {. * 100} %>% format(digits = 3), 
    "% of the population in SSA lives in countries with seasonality.\n")


closeAllConnections()
