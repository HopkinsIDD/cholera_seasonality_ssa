# This is a short script to make the repo directory tree for the README
# based on https://stackoverflow.com/questions/36094183/how-to-build-a-dendrogram-from-a-directory-tree

library(data.tree)
library(tidyverse)

paths <- list.dirs() %>% 
  # Keep only non-hidden dirs
  str_subset("\\./\\.", negate = T) %>%
  # Add all files
  c(., dir(full.names = T, recursive = T)) %>% 
  unique() %>% 
  str_subset("Rproj|README|LICENSE", negate = T)

# Tibble with paths
df <- tibble(path_string = paths)

# Comments describing contents
comments <- tribble(
  ~path_string, ~comment,
  "./analysis", "# Code for analysis",
  "./analysis/make_output_metadata.R", "## Script to make output metadata",
  "./analysis/make_repo_tree.R", "## Script to make the repo tree",
  "./generated_data", "# All analysis outputs",
  "./generated_data/best_seasonality_estimates_schema.json", "## Final seasonality estimates metadata",
  "./generated_data/best_seasonality_estimates.csv", "## Final seasonality estimates"
) %>%
  # Left align
  mutate(comment = str_pad(comment, 
                           width = max(map_dbl(comment, nchar)), 
                           side = "right") %>% 
           str_c("   ", .))

# Cast into data.tree
path_tree <- df %>% 
  left_join(comments) %>% 
  data.tree::as.Node(pathName = "path_string")

# Print to copy-paste into README
print(path_tree, "comment")
