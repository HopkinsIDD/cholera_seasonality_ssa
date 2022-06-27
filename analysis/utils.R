# Paths -------------------------------------------------------------------

getCholeraDirectory <- function(default_path = c("./")) {
  if (!exists("cholera_directory")) {
    return(default_path)
  } else {
    return(cholera_directory)
  }
}

#' Get directorry to geodata in shp format
getGeodataDir <- function() {
  str_c(getCholeraDirectory(), "/generated_data/geodata/")
}

#' Get path to Geodata
getGeodataPath <- function() {
  str_c(getCholeraDirectory(), "/generated_data/geodata.gpkg")
}

#' @title Get data paths
#'
#' @description Get the paths of run data
#'
#' @param default_paths default paths to data
#'
#' @return a vector of paths
getDataPaths <- function(default_paths = str_c(getCholeraDirectory(), "/generated_data/")) {
  if (!exists("data_paths")) {
    return(default_paths)
  } else {
    return(data_paths)
  }
}

getLogDir <- function() {
  "./logs"
}

#' Get path to monthly data
getMonthlyDataPath <- function(default_path = str_c(getCholeraDirectory(), 
                                                     "/generated_data/monthly_data_20211215.csv")) {
  if (!exists("path_to_monthly_data")) {
    return(default_path)
  } else {
    return(path_to_monthly_data)
  }
}

#' @description Get result paths
#' 
#' @param default_paths the vector of default paths to results
#'
#' @return return
getResPaths <- function(what = "outputs",
                        default_path = str_c(getCholeraDirectory(), "/generated_data/")) {
  
  if (!exists("res_paths"))  {
    return(str_c(default_path, what))
  } else {
    return(res_paths)
  }
}

# Listings ----------------------------------------------------------------

#' Get all countries
#'
#' @return a vector of country ISO3 codes in SSA
#' 
getAllCountries <- function() {
  
  readRDS(str_c(getCholeraDirectory(), "/data/ssa_countries.rds"))
}

#' Get all countries
#'
#' @return a vector of country ISO3 codes in SSA
#' 
getAllSSACountries <- function() {
  
  if (checkGeopckg()) {
    conn <- connectToGeodata()
    
    dat <- DBI::dbGetQuery(conn, "SELECT country FROM gadm36_lev0") %>% 
      dplyr::pull(country)
    
    DBI::dbDisconnect(conn)
  } else {
    dat <- st_read(glue::glue("{getGeodataDir()}/gadm36_lev0.geojson")) 
    dat <- dat$country
  }
  return(dat)
}

getAllCountryNames <- function() {
  getGadmNames(tibble::tibble(country = getAllCountries(filter_data = F), 
                              mean = NA,
                              variant = NA,
                              run_level = NA)) %>% 
    dplyr::pull(name)
  
}

#' @title Get all Covar
#'
#' @description Get all covariates
#'
#' @return df
getAllCovar <- function(type = c("climatology", "time_covar")) {
  conn <- connectToDB()
  covars <- DBI::dbGetQuery(conn, 
                            glue::glue_sql("
                            SELECT DISTINCT covar, type 
                            FROM covar_metadata
                            WHERE type IN ({type*});", .con = conn))
  DBI::dbDisconnect(conn)
  return(covars)
}

#' @title get all identifiers
#'
#' @description get all identifiers
#'
#' @return vector of model names
getAllIdentifiers <- function() {
  c("occurrence", "cases", "mean_annual_incidence")
}

#' @title get all models
#'
#' @description get all models names
#'
#' @return vector of model names
getAllModels <- function() {
  c("null", "base", "mixture", "offset", "mixture_offset")
}

getAllRunLevels <- function() {
  c(1, 2)
}

getCountryIndex <- function(country, 
                            filter_data = F) {
  which(getAllCountries(filter_data = filter_data) == country)
}



getDataMapping <- function(country,
                           variant,
                           run_lev,
                           identifier,
                           time_left = NULL,
                           time_right = NULL, 
                           gadm_lev = setGadmLev(country, run_level)) {
  # Make fake config
  config <- list(country_iso3 = country,
                 start_date = time_left,
                 end_date = time_right,
                 gadm_levels = gadm_lev,
                 run_level = run_lev,
                 variant = variant,
                 thresh = identifier,
                 keep = T)
  
  config$drop_multi_yr <- setDoMapV2(config)
  
  # Read data_mapping
  makeMappingFile(opt = list(cholera_dir = getCholeraDirectory()), 
                  config = config, 
                  run_level = run_lev) %>% 
    readRDS()
  
}

getLambdas <- function(country,
                       identifier,
                       run_level, 
                       variant,
                       time_left = NULL,
                       time_right = NULL, 
                       gadm_lev = setGadmLev(country, run_level),
                       redo = F) {
  
  run_level <- as.integer(run_level)
  
  data_mapping <- getDataMapping(country = country,
                                 variant = variant,
                                 run_lev = run_level,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 gadm_lev = gadm_lev)
  
  # Get betas to define group labelling
  betas_grp <- getParamPosteriors(country = country,
                                  run_level = run_level,
                                  model = variant,
                                  identifier = identifier,
                                  time_left = time_left,
                                  time_right = time_right,
                                  gadm_lev = gadm_lev,
                                  parameter = "betas_g", 
                                  redo = list(main = redo),
                                  verbose = F)
  
  order_groups <- getOrderGroups(betas_grp) 
  
  lambdas <- getParamPosteriors(country = country,
                                run_level = run_level,
                                model = variant,
                                identifier = identifier,
                                time_left = time_left,
                                time_right = time_right,
                                gadm_lev = gadm_lev,
                                parameter = "lambda", 
                                redo = list(main = redo),
                                verbose = F)   %>% 
    mutate(unit = str_extract(param, "(?<=\\[)[0-9]+") %>% as.numeric(),
           group = ifelse(mean > .5, order_groups[1], order_groups[2]),
           gid = data_mapping$u_unit[unit]) %>%
    select(mean, group, unit, gid) %>% 
    rename(lambda = mean)
  
  return(lambdas)
}
#' @title Get metadata
#'
#' @description gets covariate metadata
#'
#' @param redo
#'
#' @return df
getMetadata <- function(redo = F) {
  
  if (redo == T) {
    makeMetadata()
  }
  
  DBI::dbReadTable(conn, "covar_metadata")
}

# Spatial operations ------------------------------------------------------

checkGeopckg <- function() {
  Sys.getenv("GEOPKG") == "TRUE"
}

#' Connect to database
#'
#' @description Connects to the postgres database for data processing
#'
#' @param username
#' @param database
#'
#' @return a connection obejct
#' 
connectToDB <- function(username = "perez",
                        database = "taxdat") {
  
  conn <- DBI::dbConnect(RPostgres::Postgres(), 
                         user = username, 
                         dbname = database)
  return(conn)
}


#' Connect to Geodata
#'
#' @return an SQLITE connection
#'
connectToGeodata <- function() {
  DBI::dbConnect(RSQLite::SQLite(), getGeodataPath())
}


#' Get GADM data
#'
#' @description Get population and cholera case data aggregated at admin unit level
#'
#' @param cholera_directory
#'
#' @return return
#' 
getGadmData <- function(cholera_directory = getCholeraDirectory()) {
  
  gadm_data <- readRDS(paste0(cholera_directory, "/generated_data/gadm_data.rds")) 
  
  return(gadm_data)
}

#' Get GADM sf object
#'
#' @param country 
#' @param gadm_lev 
#' @param cols 
#'
#' @return an sf object
#'
getGadmSf <- function(country, 
                      gadm_lev, 
                      cols = c("country", "gid"),
                      geodata = getGeodataPath(),
                      geodatadir = getGeodataDir()) {
  
  if (checkGeopckg()) {
    conn <- connectToGeodata()
    
    dat <- st_read(dsn = getGeodataPath(),
                   query = glue::glue_sql(
                     "SELECT *
            FROM gadm36_lev{gadm_lev}
            WHERE country IN ({country*})",
                     .con =  conn))
    
    DBI::dbDisconnect(conn)
  } else {
    cntry <- country
    dat <- st_read(glue::glue("{getGeodataDir()}/gadm36_lev{gadm_lev}.geojson")) 
    dat <- dat[dat$country %in% cntry,]
    names(dat)[names(dat) == attr(dat, "sf_column")] <- "geom"
    st_geometry(dat) <- "geom"
  }
  
  dat %>% 
    {
      if (cols[1] != "all") {
        dplyr::select(., one_of(cols))
      } else {
        .
      }
    }
}


#' Get GADM list
#'
#' @param country 
#' @param gadm_lev
#' @param geodata 
#'
#' @return a dataframe
#' 
getGadmList <- function(country,
                        gadm_lev = 0:2,
                        geodata = getGeodataPath()) {
  
  if (checkGeopckg()) {
    conn <- connectToGeodata()
    
    gadm_list <- purrr::map_df(
      gadm_lev,
      ~ DBI::dbGetQuery(conn,
                        glue::glue_sql(
                          "SELECT country, {.} as gadm_lev, gid
                        FROM gadm36_lev{.}
                        WHERE country IN ({country*})",
                          .con = conn) 
      )) %>% 
      tibble::as_tibble()
    
    DBI::dbDisconnect(conn)
    
  } else {
    cntry <- country
    gadm_list <- purrr::map_df(
      gadm_lev,
      function(x) {
        st_read(glue::glue("{getGeodataDir()}/gadm36_lev{x}.geojson")) %>% 
          st_drop_geometry() %>% 
          as_tibble() %>% 
          filter(country %in% cntry) %>% 
          mutate(gadm_lev = x) %>% 
          select(country, gadm_lev, gid)
      })
  }
  
  gadm_list
}

#' Get mapping
#'
#' @param country 
#' @param geodata 
#'
#' @returna dataframe
#'
getGadmMapping <- function(country, 
                           geodata = getGeodataPath()) {
  
  if (checkGeopckg()) {
    conn <- connectToGeodata()
    
    df <- DBI::dbGetQuery(conn,
                          glue::glue_sql(
                            "SELECT *
                        FROM gadm36_lev1lev2_mapping
                        WHERE country IN ({country*})",
                            .con = conn)) %>% 
      dplyr::as_tibble()
    
    DBI::dbDisconnect(conn)
  } else {
    cntry <- country
    df <- readRDS(str_glue("{getGeodataDir()}/gadm36_lev1lev2_mapping.rds")) %>% 
      filter(country %in% cntry)
  }
  
  return(df)
}

getGadmNames <- function(df) {
  conn <- connectToDB()
  res <- df %>% 
    dplyr::group_by(country) %>% 
    dplyr::group_map(function(x, y) {
      admin_lev <- ifelse(is.na(x$mean[1]) | is.na(x$variant[1]), 0, x$run_level[1])
      gadm <- st_read(dsn = conn,
                      query = glue::glue("SELECT gid_{admin_lev}, name_{admin_lev}
      FROM admin.gadm_lev{admin_lev}
                         WHERE gid_0 = '{y$country[1]}'"))
      gadm %>% 
        magrittr::set_colnames(c("gid", "name")) 
    }) %>% 
    dplyr::bind_rows() 
  
  DBI::dbDisconnect(conn)
  res
}

#' @title Get Latitudes
#'
#' @description Gets the latitudes of the centroid of countries
#'
#' @param countries
#' @param redo
#'
#' @return a df
#' 
getLatitude <- function(countries = "all",
                        redo = F) {
  
  lat_file <- makeStdResName(country = countries,
                             run_level = "all",
                             model = "all",
                             identifier = "all",
                             suffix = "latitutde",
                             file_type = "csv")
  
  if (file.exists(lat_file) & !redo) {
    latitudes <- readr::read_csv(lat_file,
                                 col_types = readr::cols())
    
  } else {
    
    if (countries == "all") {
      countries <- getAllCountries()
    }
    
    # Get admin 0 shapefiles
    africa.sf_lev0 <- purrr::map_df(countries,
                                    ~ getShapefile(., 0))
    # Extract latitudes
    latitudes <- tibble::tibble(
      country = africa.sf_lev0$gid,
      latitude = sf::st_coordinates(sf::st_centroid(africa.sf_lev0))[, 2],
      longitude = sf::st_coordinates(sf::st_centroid(africa.sf_lev0))[, 1]
    )
    
    # Write result
    readr::write_csv(latitudes, file = lat_file)
  }
  
  return(latitudes)
}


getMeanGadmAreas <- function() {
  outfile <- str_c(getCholeraDirectory(), "generated_data/geodata/mean_areas.rds")
  
  if (!file.exists(outfile)) {
    
    gadm_areas <- purrr::map_df(1:2,
                                function(x) {
                                  getGadmSf(country = getAllCountries(),
                                            gadm_lev = x) %>% 
                                    mutate(gadm_lev = x)})  %>% 
      # Area in sqkm
      mutate(area = sf::st_area(geom) %>% as.numeric()/1e6)
    
    mean_areas <- gadm_areas %>%
      st_drop_geometry() %>% 
      group_by(country, gadm_lev) %>% 
      summarise(median_area = median(area)) %>% 
      group_by(gadm_lev) %>% 
      summarise(mean_area = mean(median_area))
    
    saveRDS(mean_areas, file = outfile)
    
  } else {
    mean_areas <- readRDS(outfile)
  }
  return(mean_areas)
}

#' @title Get shapefile
#'
#' @description Get the shapefile for a given country at a given admin level
#'
#' @param country
#' @param admin_level
#' @param as_df
#' @param source
#' @param source_path
#' 
#' @return an sf object with the country shapefile
#' 
getShapefile <- function(country,
                         admin_level,
                         as_df = FALSE,
                         source = "sql",
                         source_path = paste0(getCholeraDirectory(), "/data")) {
  
  if (!(source %in% c("sql", "csv"))) {
    stop("Source ", source, " not recognized, needs to be either sql or csv")
  }
  
  if (source == "sql") {
    conn <- connectToDB()
    if (!as_df) {
      query <- glue::glue_sql("SELECT gid_{admin_level} as gid,
                              gid_0 as country, ST_SimplifyPreserveTopology(geom, .0075) as geom, area
                              FROM admin.gadm_lev{admin_level}
                              WHERE gid_0 = {country*}",
                              .con = conn)
      
      sf_object <- sf::st_read(conn, query = query) %>% 
        dplyr::mutate(gadm_lev = admin_level)
      
      DBI::dbDisconnect(conn)
      
      return(sf_object)
    } else {
      query <- glue::glue_sql("SELECT gid_{admin_level} as gid,
                              gid_0 as country, area
                              FROM admin.gadm_lev{admin_level}
                              WHERE gid_0 = {country*}",
                              .con = conn)
      
      df <- DBI::dbGetQuery(conn, query) %>% 
        dplyr::mutate(gadm_lev = admin_level)
      
      DBI::dbDisconnect(conn)
      
      return(df)
    }
  } else {
    if (!as_df) {
      stop("Reading sf objects from csv is not allowed")
    } else {
      df <- dir(source_path, full.names = T, pattern = "africa_gadm") %>% 
        stringr::str_subset(paste0(admin_level, ".csv")) %>% 
        purrr::map_df(~readr::read_csv(., col_types = readr::cols()))
      
      df <- df[df$cntry %in% country, ]
      
      return(df)
    }
  }
}

#' Get GADM unit name
#' Pulls the name for a given gadm unit and level from the postgres database
#' 
#' @param gid this can be a vector as well
#' @param gadm_lev this has to be a single level
#'
#' @return a string with the name
#' 
getGadmName <- function(gid, 
                        gadm_lev) {
  conn <- connectToDB()
  DBI::dbGetQuery(conn, 
                  glue::glue_sql("SELECT gid, name
                                       FROM admin.gadm36_lev{gadm_lev}
                                 WHERE gid IN ({gid*})", .con = conn)) %>% 
    dplyr::as_tibble()
}

# Misc --------------------------------------------------------------------

#'  Make period
#' @description Make obervation period
#'
#' @param opt
#' @param time_left
#' @param time_right
#'
#' @return a string
#' 
makePeriod <- function(opt = NULL,
                       time_left = NULL, 
                       time_right = NULL) {
  
  if (!is.null(opt)) {
    str_c(
      ifelse(is.null(opt$start_date), "all", opt$start_date),
      ifelse(is.null(opt$end_date), "all", opt$end_date),
      sep = "-"
    )
  } else {
    paste(ifelse(is.null(time_left), "all", time_left), 
          ifelse(is.null(time_right), "all", time_right), 
          sep = "-")
  }
}
