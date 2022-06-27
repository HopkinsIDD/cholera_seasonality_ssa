# This is a short script to make the repo directory tree for the README
# based on https://stackoverflow.com/questions/36094183/how-to-build-a-dendrogram-from-a-directory-tree

library(data.tree)
library(tidyverse)

paths <- list.dirs() %>% 
  # Keep only non-hidden dirs
  str_subset("\\./\\.", negate = T) %>%
  # Add all files
  c(., dir(full.names = T, recursive = T)) %>% 
  unique()

# Tibble with paths
df <- tibble(path_string = paths)

# Comments describing contents
comments <- tribble(
  ~path_string, ~comment,
  "./analysis/", "",
  "./analysis/01_data_preparation", "-- (not on this repo).",
  "./analysis/02_modeling/configs", "--- Configs used for model running",
  "./analysis/02_modeling/configs/cases", "---- Cases as thershold for cholera presence",
  "./analysis/02_modeling/configs/mean_annual_incidence", "---- Mean annual incidence as threshold",
  "./analysis/02_modeling/configs/occurrence", "---- Occurrence as thershold",
  "./analysis/02_modeling/configs/default_config.yml", "---- Annotated example config",
  "./analysis/02_modeling/01_run_model_prior.R", "--- Runs the seasonality model with priors only",
  "./analysis/02_modeling/02_run_model.R", "--- Runs the seasonality model",
  "./analysis/02_modeling/03_generate_simulations.R", "--- Generates simulated outputs from model fits",
  "./analysis/03_postprocessing/01_extract_prior_results_config.R", "- Extracts outputs from prior model fits",
  "./analysis/03_postprocessing/02_extract_fitting_results_config.R", "- Extracts outputs from model fits",
  "./analysis/03_postprocessing/03_compute_seasonality_index.R", "- Computes the seasonality index",
  "./analysis/03_postprocessing/04_postprocess_outputs.R", "- Postprocesses outputs for reporting",
  "./analysis/04_clustering/01_compute_grouping.R", "- Runs clustering model",
  "./analysis/04_clustering/02_extract_clustering_results.R", "- Extracts clustering outputs",
  "./analysis/05_outputs/01_final_stats.R", "- Makes statistics statements used in manuscript",
  "./analysis/05_outputs/02_final_plots.R", "- Makes all figures in manuscript",
  "./analysis/stan", "- All stan code",
  # "./analysis/stan/binomial_bym2_fast_convolved_lambda_flexmixture_speedup.stan", "--- Main stan model for analysis (see 'Making of' notebook for details)",
  "./analysis/clustring_utils.R", "- Helper functions for clustering",
  "./analysis/make_output_metadata.R", "-- Script to make output metadata",
  "./analysis/make_repo_tree.R", "-- Script to make the repo tree",
  "./analysis/modeling_utils.R", "- Helper functions for seasonality models",
  "./analysis/nbgraph_utils.R", "- Helper functions for spatial adjacency",
  "./analysis/plotting_utils.R", "- Helper functions for plotting",
  "./analysis/postprocessing_utils.R", "- Helper functions for postprocessing",
  "./analysis/simulation_utils.R", "- Helper functions for simulations",
  "./analysis/distribution_utils.R", "- Helper functions for preparing this repo",
  "./analysis/utils.R", "- Misc helper functions",
  "./generated_data/best_seasonality_estimates_schema.json", "-- Final seasonality estimates metadata",
  "./generated_data/best_seasonality_estimates.csv", "-- Final seasonality estimates",
  "./generated_data/run_data", "-- Data used to run the models",
  "./generated_data/run_data/all_run_data_mean_annual_incidence.csv", "--- Compilation of all run data for the main analysis",
  "./generated_data/run_data/all_run_data_schema.json", "--- Metadata of the run data",
  "./generated_data/data_mapping", "-- Formatted data for stan ",
  "./generated_data/geodata", "-- Geographic data used in the analysis",
  "./generated_data/postprocessing", "-- Postprocessed data to reproduce figures and tables",
  "./generated_data/model_outputs", "-- Modeling outputs used to reproduce figures and tables",
  "./figures", "",
  "./figures/pub_figures", "Figures in main text and supplement",
  "./manuscript", "",
  "./manuscript/supplement/", "",
  "./manuscript/supplement/supplement_cholera_seasonality.Rmd", "--- RMD to produce supplementary material pdf",
  "./manuscript/supplement/supplement_cholera_seasonality.pdf", "--- Supplementary material pdf",
  "./manuscript/supplement/refs.bib", "--- References used in the supplement",
  "./notebooks/project_details.Rmd", "--- RMD of details on project implementation",
  "./notebooks/project_details.html", "--- HTML of details on project implementation"
) %>%
  rowwise() %>% 
  # Left align
  mutate(comment = str_remove_all(comment, "-") %>% str_trim(),
         comment = str_c(str_c(rep("-", str_count(path_string, "/")), collapse = ""), " ", comment)
  ) %>% 
  ungroup() %>% 
  mutate(comment = str_remove_all(comment, "- $"),
         comment = str_pad(comment,
                           width = max(map_dbl(comment, nchar)),
                           side = "right") %>%
           str_c("   ", .))

# Cast into data.tree
path_tree <- df %>% 
  inner_join(comments) %>% 
  arrange(path_string) %>% 
  data.tree::as.Node(pathName = "path_string")

# Print to copy-paste into README
print(path_tree, "comment")

