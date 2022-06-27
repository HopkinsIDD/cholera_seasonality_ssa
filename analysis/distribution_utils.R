# These are functions that were not used in the analysis but only for 
# preparing the repo for distribution


compileData <- function(data_path = str_c(getDataPaths(), "run_data")) {
  
  # Get all data files
  files <- dir(data_path, full.names = T, pattern = "rds")
  
  # Compile
  all_data <- map_df(files, function(x) {
    dat <- readRDS(x)
    dat %>% 
      mutate(country = stringr::str_extract(x, "[A-Z]{3}"),
             identifier = stringr::str_extract(x, 
                                               stringr::str_c(getAllIdentifiers(), collapse = "|")))
  })
  
  walk(getAllIdentifiers(), 
       function(x) {
         all_data %>% 
           dplyr::filter(identifier == x) %>% 
           dplyr::select(country, identifier, gid, gadm_lev, year, month_left, month_right, n_obs, cholera) %>% 
           readr::write_csv(stringr::str_glue("generated_data/all_run_data_{x}.csv"))
       })
}

