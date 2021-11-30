# This script is used to produce metadata for all project outputs
# The metadata is based on the Frictionless Standards table template (https://specs.frictionlessdata.io/table-schema/)
# implemented in the R package tableschema.r 


# Preamble ----------------------------------------------------------------

library(tableschema.r)
library(jsonlite)
library(tidyverse)


#' Build Schema
#' Function to build a metadata schema using the Frictionless Standards template
#' 
#' @param path path to the data
#' @param tibble list containing the details 
#' 
#' @details The function writes a json file with the metadata based on the 
#' file path
#' 
#' @return nothing
#' 
buildSchema <- function(path, 
                        details_dict) {
  
  # Load table
  df <- Table.load(path) %>% 
    future::value()
  
  df$infer(limit = 2)
  df$schema$descriptor['missingValues'] = 'NA'
  
  # Add title and comment field
  df$schema$descriptor$fields <- df$schema$descriptor$fields %>% 
    map(~ c(., list(title = "", description = "")))
  
  df$schema$commit()
  
  # Insert values
  df$schema$descriptor$fields <- df$schema$descriptor$fields %>% 
    map(function(x) {
      x$title[1] <- details_dict$title[details_dict$name == x$name] %>% unlist()
      x$description[1] <- details_dict$description[details_dict$name == x$name] %>% unlist()
      x
    })
  
  schema_fname <- str_replace(path, "\\.csv", "_schema.json")
  
  tableschema.r::write_json(df$schema$descriptor, schema_fname)
  
  cat("Wrote schema to", schema_fname, "\n")
}

# Seasonality estimates ---------------------------------------------------

seas_dict <- tribble(
  ~name, ~title, ~description,
  "country", "Country ISOcode", "Country ISO 3 code",
  "gid", "GADM id", "Unique id in the GADM dataset",
  "name", "GADM name", "Name in the GADM dataset",
  "run_level", "Run level", "Administrive level at which the model was run depending on data availability (either 1 or 2)",
  "model", "Best model", "Model variant which was retained by model comparison (one of null, base, offset, mixture or mixture_offset). See Supplementary material for descriptions.",
  "month", "Month", "Month of the year",
  "mean", "Mean estimate", "Posterior mean of the seasonality coefficient (odds ratio of observing excess cholera in a given month. For the null model values were set to NA. Estimates are GADM unit-specific for the mixture and mixture_offset models, and country-specific for the base and offset models.",
  "q025", "2.5% quantile", "Posterior 2.5% quantile of the seasonality coefficient",
  "q975", "97.5% quantile", "Posterior 2.5% quantile of the seasonality coefficient"
)

buildSchema(path = "generated_data/best_seasonality_estimates.csv",
            details_dict = seas_dict)



# Update README -----------------------------------------------------------
file.remove("generated_data/README.md")
sink(file = "generated_data/README.md", append = T)
cat("# Metadata \n")
cat("This file specifies metadata for all generated outputs found in `./generated_data`.",
"Metadata follows the Frictionless Standards template in https://specs.frictionlessdata.io/table-schema/ \n\n")

dir("generated_data", pattern = "json", full.names = T) %>% 
 walk(function(x) {
   cat("## ", x, "\n\n ```json \n", sep = "")
   Schema.load(x) %>% 
       future::value() %>% 
       .$descriptor %>% 
       toJSON(pretty = T) %>%
     print()
   cat("``` \n")
   })

closeAllConnections()

