# This script generates simulated outputs based on model fits
# 
# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(optparse)

# User-supplied options
option_list <- list(
  make_option(c("-c", "--config"), default = "./analysis/02_modeling/configs/default_config.yml", 
              action ="store", type = "character", help = "Configuration file"),
  make_option(c("-d", "--cholera_dir"), default = "./",
              action ="store", type = "character", help = "Directory"),
  make_option(c("-r", "--redo"), default = FALSE, type = "logical")
)

# Parse options
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

source(str_c(opt$cholera_dir, "/analysis/utils.R"))
source(str_c(opt$cholera_dir, "/analysis/modeling_utils.R"))
source(str_c(opt$cholera_dir, "/analysis/postprocessing_utils.R"))

# Load model configuration file
config <- yaml::read_yaml(opt$config)

cat("Running ", opt$config, "\n")

# Extract -----------------------------------------------------------------

if (is.na(as.numeric(config$run_level))) {
  stop("Please specify run level")
}
# Rstan fit file 
rstan_file <- makeResultsFile(opt = opt, 
                              cntry_iso3 = config$country_iso3, 
                              run_level = as.numeric(config$run_level), 
                              drop_multi_yr = setDoMapV2(config)) 

fit_file <- makeCmdstanrFile(rstan_file)

# Check if output file for simulations already exists
sim_file <- makeCmdstanrSimFile(rstan_file)
output_sim_file <- makeOutputSimFile(rstan_file)

if (file.exists(output_sim_file) & !opt$redo) {
  stop("Generated data already computed")
}

# Data mapping with additional information on unique units
data_mapping <- readRDS(makeMappingFile(opt = opt, 
                                        config = config, 
                                        run_level = as.numeric(config$run_level)))

data <- rjson::fromJSON(file = makeJsonFile(opt, config))

if (!file.exists(fit_file)) {
  stop("Stan fit object ", fit_file, " not found.")
}

if (file.exists(sim_file) & !opt$redo) {
  gen_quant <- readRDS(sim_file)
} else {
  
  # Load cmdstanr fit file
  res <- readRDS(fit_file)
  
  # Extract offset probabilities if any
  if (str_detect(fit_file, "offset")) {
    
    data$offset_probs <- makeOutputFile(rstan_file)%>%
      readRDS(.) %>% 
      .$offset_probs
    
    # Handle case where offset probabilities were misspecified
    if (round(sum(data$offset_probs)) == 12) {
      data$offset_probs <- data$offset_probs/12
    }
    
  } else {
    data$offset_probs <- rep(1, 12)
  }
  
  # Generate observations ---------------------------------------------------
  
  if (str_detect(fit_file, "null")) {
    if ("tau" %in% res$metadata()$stan_variables) {
      stan_model <- "analysis/stan/binomial_bym2_fast_null_speedup_tau_generate.stan"
    } else {
      stan_model <- "analysis/stan/binomial_bym2_fast_null_speedup_sigma_generate.stan"
    }
  } else {
    if (config$country_iso3 %in% c("NGA", "COD") & config$run_level == 2) {
      stan_model <- "analysis/stan/binomial_bym2_fast_convolved_lambda_flexmixture_speedup_map_rect_v2_generate.stan"
    } else {
      stan_model <- "analysis/stan/binomial_bym2_fast_convolved_lambda_flexmixture_speedup_tau_generate.stan"
    }
  }
  
  gen_model <- cmdstan_model(stan_model,
                             quiet = FALSE,
                             force_recompile = F)
  
  if ((config$country_iso3 == "NGA" & str_detect(config$variant, "offset")) |
      (config$country_iso3 == "COD" & str_detect(config$variant, "mixture_o"))) {
    fitted_params <- res$post_warmup_draws
  } else if (config$country_iso3 == "NGA") {
    fitted_params <- res$draws()
    fitted_params <- fitted_params[sample(1:1000, 85), c(1, 2, 4),]
  } else if (config$country_iso3 == "COD") {
    fitted_params <- res$draws()
    fitted_params <- fitted_params[sample(1:1000, 85), ,]
  }  else {
    fitted_params <- res
  }
  
  gen_quant <- gen_model$generate_quantities(fitted_params = fitted_params,
                                             data = data,
                                             parallel_chains = 4)
  
  # Save cmdstanr object
  gen_quant$save_object(file = sim_file)
}
# Outputs -----------------------------------------------------------------

# Simulated quantiles
sim_stats <- gen_quant$summary(variables = "sim", ~ 
                                 c(mean = mean(.), 
                                   posterior::quantile2(., probs = c(.025, seq(.05, .95, by = .05), .975)))) %>% 
  mutate(obs = data$y)

# Time series of predicted probabilities
ts_prob_stats <- gen_quant$summary(variables = "ts_prob") %>% 
  mutate(unit = data_mapping$u_unit[data_mapping$uunit],
         month = data_mapping$umonths,
         year = data_mapping$u_years[data_mapping$uyears])

# Seaonality coefficients by admin unit
if (!str_detect(fit_file, "null")) {
  beta_stats <- gen_quant$summary(variables = "betas", ~ 
                                    c(mean = mean(.), 
                                      median = median(.), 
                                      posterior::quantile2(., probs = c(.025, .25, .75, .975)))) %>% 
    mutate(unit = data_mapping$u_unit[str_extract(variable, "(?<=\\[)[0-9]+(?=\\,)") %>% as.numeric()],
           month = str_extract(variable, "(?<=\\,)[0-9]+(?=\\])") %>% as.numeric())
} else {
  beta_stats <- NULL
}

# Combine and save
output <- list(
  sim_stats = sim_stats,
  ts_prob_stats = ts_prob_stats,
  beta_stats = beta_stats
)

saveRDS(output, file = output_sim_file)
