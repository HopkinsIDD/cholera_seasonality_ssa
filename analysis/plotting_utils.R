# Utils -------------------------------------------------------------------

addSpatialAnnotations <- function(x) {
  x + 
    ggspatial::annotation_north_arrow(location = "bl",
                                      pad_y = unit(.75, "cm"),
                                      pad_x = unit(1.25, "cm"),
                                      height = unit(.8, "cm"),
                                      width = unit(.8, "cm")) +
    ggspatial::annotation_scale(location = "bl",
                                line_width = .7,
                                height = unit(0.15, "cm")) 
}

## This function allows us to specify which facet to annotate
# https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
annotation_custom2 <- function (grob,
                                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                                data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = F, params = list(grob = grob, 
                                       xmin = xmin, xmax = xmax, 
                                       ymin = ymin, ymax = ymax))
}


#' @title Get Climate colors
#'
#' @description Get the color scales used to plot climate variables
#'
#' @return a list with high and low colors for scale_gradient_2
#' 
getClimColors <- function() {
  list(
    Precipitation = list(low = "white",
                         high = "darkblue"),
    Mean_T = list(low = "yellow",
                  high = "darkred")
  )
}


getCovarDict <- function() {
  c("Mean_T" = "Mean air\ntemperature",
    "Max_T" = "Maximum air\ntemperature",
    "Min_T" = "Minimum air\ntemperature",
    "Precipitation" = "Cumulative \nprecipitation",
    "FloodMean" = "Mean flooded \narea",
    "FloodMax" = "Maximum flooded \narea",
    "Runoff" = "Runoff",
    "SoilMoisture" = "Soil moisture\nin top 30cm",
    "Streamflow" = "Streamflow")
}

#' @title Get Cluster colors
#'
#' @description get the color palette to plot seasonality clusters
#'
#' @param k the number of colors
#'
#' @return a vector with colors
#' 
getClusterColors <- function(k) {
  if (k == 5) {
    c("#8F6881", "#B80D6B", "#EE6AA7", "#5DBA00", "#A7F045")
  } else if (k == 4){
    c("#B80D6B", "#EE6AA7", "#5DBA00", "#A7F045")
  } else {
    c("#664A5C", "#8B0A50",  "#458B00")
  }
}

#' @title Get Cluster colors
#'
#' @description get the color palette to plot seasonality clusters
#'
#' @param k the number of colors
#'
#' @return a vector with colors
#' 
getNClusterColors <- function(k = 2) {
  all_cols <- RColorBrewer::brewer.pal(8, "Dark2")  
  all_cols[1:k]
}

#' @title Get excluded countries
#'
#' @description Get the countries that were not run because of lack of data
#'
#' @param identifier
#' @param run_level
#'
#' @return a vector of countries
#' 
getExcludedCountries <- function(identifier = "mean_annual_incidence",
                                 run_level = "best") {
  all_countries <- getAllSSACountries()
  ran_countries <- getRanCountries()
  
  setdiff(all_countries, ran_countries)
}

getRanCountries <- function(identifier = "mean_annual_incidence",
                            run_level = "best") {
  getAllBestModels(countries = "all",
                   run_levels = run_level,
                   identifier = identifier,
                   redo = F,
                   redo_single = F) %>% 
    dplyr::pull(country)
}

#' @title Get no seasonality countries
#'
#' @description Get the countries for which the null was selected
#'
#' @param identifier
#' @param run_level
#'
#' @return a vector of countries
#' 
getNoSeasCountries <- function(identifier = "mean_annual_incidence",
                               run_level = "best") {
  
  ran_countries <- getAllBestModels(countries = "all",
                                    run_levels = run_level,
                                    identifier = identifier,
                                    redo = F,
                                    redo_single = F)
  
  ran_countries$country[ran_countries$best_model == "null"]
}


#' @title Get no offset countries
#'
#' @description Get the countries for which the offset was not selected
#'
#' @param identifier
#' @param run_level
#'
#' @return a vector of countries
#' 
getNoOffsetCountries <- function(identifier,
                                 run_level = "best") {
  
  ran_countries <- getAllBestModels(countries = "all",
                                    run_levels = run_level,
                                    identifier = identifier,
                                    redo = F,
                                    redo_single = F)
  
  ran_countries$country[!stringr::str_detect(ran_countries$best_model, "offset")]
}

#' @title Get figures directory
#'
#' @description Get the directory where to save figures for publication
#'
#' @return a string
#' 
getFigDir <- function() {
  figdir <- "figures/pub_figures"
  
  if (!dir.exists(figdir)) {
    dir.create(figdir)
  }
  return(figdir)
}


#' @title Get colors for seasonality grouping
#'
#' @description Get the color used to plot the two groups in mixture-type models
#'
#' @return a vector with two colors 
#' 
getGroupingColors <- function() {
  c("#5933F2", "#E3B009")
}

#' @title Get Lakes
#'
#' @description Gets large waterbodies in SSA from
#' https://datacatalog.worldbank.org/dataset/africa-water-bodies-2015
#'
#' @param path to data file
#'
#' @return an sf_object
#' 
getLakes <- function(path = "data/geodata/Africa_waterbody.shp") {
  lakes_sf <- sf::st_read(path) %>% 
    dplyr::filter(Shape_area>.5) %>% 
    rmapshaper::ms_simplify(keep = 0.1,
                            keep_shapes = FALSE)
  
  lakes_sf
}

getOrderGroups <- function(betas_grp) {
  # Group A is always the one having the highest values in August
  which_group_A <- betas_grp %>% 
    dplyr::filter(str_detect(param, "8")) %>% 
    dplyr::arrange(param_fam) %>% 
    {which.max(.$mean)}
  
  if (which_group_A == 1) {
    order_groups <- c("A", "B")
  } else {
    order_groups <- c("B", "A")
  }
  order_groups
}


getModelColors <- function() {
  c("green", "red", "purple", "gray", "blue")
}

getMonthGradientCols <- function() {
  c("#FF0000", "#B01CFF", "#00F2FF", "#FFE100", "#FF0000")
}

getRateDict <- function() {
  c("(-Inf,1e-05]" = "< 1",
    "(1e-05,1e-04]" = "1-10",
    "(1e-04,1e-03]" = "10-100",
    "(1e-03,Inf]" = ">100")
}

getRateCuts <- function() {
  c(-Inf, 1e-5, 1e-4, 1e-3, Inf)
}


getSeasCatColors <- function() {
  c("#ffc9b7", "#ff6342", 
    "#dd0001", "#930001")
}


#' @title Get best sf object
#'
#' @description gets the sf objects corresponding to a set of results
#'
#' @param param
#'
#' @return return
#' 
getSfObjects <- function(df) {
  
  sf_objects <- df %>% 
    dplyr::group_by(country) %>% 
    dplyr::group_map(function(x, y) {
      admin_lev <- ifelse(is.na(x$mean[1]) | is.na(x$variant[1]), 0, x$run_level[1])
      getShapefile(y$country[1], admin_lev)
    }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(geometry = geom) 
  
  # Simplify admin units
  sf_objects <- sf_objects %>% 
    rmapshaper::ms_simplify(keep = 0.1,
                            keep_shapes = FALSE)
  
  return(sf_objects)
}

#' @title Lag months
#'
#' @description function to lag months in a cyclic manner
#'
#' @param param
#'
#' @return return
#' 
lagMonths <- function(months, lag) {
  months <- months - lag
  months[months <= 0] <- months[months <= 0]  + 12
  return(months)
}


#' @title Make figure file name
#'
#' @description Make figure file name
#'
#' @param countrie
#' @param run_level
#' @param time_left
#' @param time_right
#' @param prefix
#'
#' @return return
makeFigFilename <- function(countries,
                            run_level,
                            identifier,
                            time_left = NULL,
                            time_right = NULL,
                            prefix = "fig_") {
  
  paste0(getFigDir(), "/", prefix,
         "_",
         stringr::str_c(countries, collapse = "-"),
         "_rl",
         run_level,
         "_",
         identifier,
         "_",
         makePeriod(time_left, time_right),
         ".png")
}


#' @title Tidy best seasonality by group
#'
#' @description description
#'
#' @param best_seas
#' @param all_countries
#'
#' @return return
#' 
tidyBestSeasGrp <- function(best_seas, 
                            all_countries = getAllCountries()) {
  best_seas %>% 
    mutate(param = ifelse(is.na(param_fam), "betas_g1[1]", param)) %>% 
    tidyr::complete(country = all_countries,
                    param = str_c("betas_g1[", 1:12, "]")) %>% 
    mutate(param_fam = ifelse(is.na(param_fam), "betas_g1", param_fam),
           month = str_extract(param, "(?<=\\[)[0-9]+") %>% as.numeric(),
           param = factor(param, levels = c(str_c("betas_g1[", 1:12, "]"),
                                            str_c("betas_g2[", 1:12, "]")))) %>% 
    arrange(country, param)
}


st_shrink <- function(df, scale) {
  areas <- st_area(df) %>% as.numeric()/1e6 %>% log10()
  scales <- scale  * (1 - .015 * (areas-max(areas))/(min(areas)-max(areas)))
  df$s <- scales
  
  df %>% 
    rowwise() %>% 
    group_modify(function(x ,y ){
      dfgeom <- st_geometry(x)
      metadata <- st_drop_geometry(x) 
      
      st_as_sf((dfgeom - st_centroid(dfgeom)) * x$s + st_centroid(dfgeom)) %>% 
        bind_cols(metadata) %>% 
        rename(geometry = x) %>% 
        st_set_geometry("geometry") %>% 
        st_set_crs(4326)
    })
  
}

getCountrySf <- function(all_countries) {
  purrr::map_df(all_countries,
                ~ getShapefile(., 0)) %>% 
    rmapshaper::ms_simplify(keep = 0.05,
                            keep_shapes = FALSE)}

# Internal use ------------------------------------------------------------

#' @title Plot betas
#'
#' @description Plot the seasonality coefficients
#'
#' @param betas_df
#'
#' @return a ggplot
#' 
plotBetas <- function(betas_df = NULL,
                      country = NULL,
                      run_level = NULL,
                      identifier = "all",
                      time_left = NULL, 
                      time_right = NULL,
                      model = "all",
                      gadm_lev = setGadmLev(country, run_level),
                      redo = F) {
  
  runChecks(country = country,
            run_level = run_level,
            identifier = identifier,
            model = model)
  
  if (is.null(betas_df)) {
    if (is.null(country) | is.null(run_level))
      stop("Please provide country, run_level and identifier for model ic extraction")
    
    if (identifier == "all") {
      identifiers <- getAllIdentifiers()
    } else {
      identifiers <- identifier
    }
    
    betas_df <- purrr::map_df(identifiers,
                              ~ getParamPosteriors(country = country,
                                                   run_level = run_level,
                                                   identifier = .,
                                                   time_left = time_left, 
                                                   time_right = time_right,
                                                   gadm_lev = gadm_lev,
                                                   redo = list(main = redo),
                                                   parameter = "betas"))
  } 
  
  betas_df <- betas_df %>% 
    dplyr::mutate(month = as.numeric(str_extract(param, "(?<=\\[)[0-9]+")))
  
  mbreaks <- c(1, 3, 6, 9, 12)
  p <- ggplot2::ggplot(betas_df,
                       ggplot2::aes(x = month, 
                                    y = mean, 
                                    lty = param_fam)) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975, fill = variant),
                         alpha = .1) +
    ggplot2::geom_line(ggplot2::aes(color = variant)) +
    ggplot2::facet_grid(variant~identifier) +
    ggplot2::scale_x_continuous(breaks = mbreaks, 
                                labels = month.abb[mbreaks]) +
    ggplot2::theme_bw()
  
  return(p)
}




plotGroupingWeights <- function(country,
                                identifier,
                                variant,
                                run_level, 
                                time_left = NULL,
                                time_right = NULL, 
                                gadm_lev = setGadmLev(country, run_level),
                                sf_object = NULL,
                                redo = F,
                                redo_single = F,
                                verbose = F) {
  
  if (!str_detect(variant, "mixture")) {
    warning("Selected model does not contain groupin")
    return(NULL)
  }
  
  if (is.null(sf_object)) {
    seas <- getBetas(country = country,
                     run_level = run_level,
                     model = variant,
                     identifier = identifier,
                     time_left = time_left,
                     time_right = time_right)
    
    sf_object <- getSfObjects(seas)
  }
  
  lambdas <- getLambdas(country = country,
                        identifier = identifier,
                        run_level = run_level, 
                        variant = variant,
                        time_left = time_left,
                        time_right = time_right)
  
  sf_object <- lambdas %>% 
    inner_join(sf_object, .)
  
  p_lambdas <- sf_object %>% 
    mutate(group_prob = abs(.5-pmax(lambda, 1-lambda))*2) %>% 
    ggplot(aes(fill = group, alpha = group_prob)) +
    geom_sf(size = .03, color = "gray") +
    scale_fill_manual(values = getGroupingColors()) +
    guides(alpha = "none") +
    ggthemes::theme_map()
  
  p_lambdas
}


plotGrouping <- function(country,
                         identifier,
                         run_level, 
                         time_left = NULL,
                         time_right = NULL, 
                         gadm_lev = setGadmLev(country, run_level),
                         redo = F,
                         verbose = F) {
  
  best_seas <- getBestSeas(countries = country, 
                           identifier = identifier,
                           run_level = run_level,
                           redo = list(main = redo,
                                       best_seas_single = redo,
                                       best_seas = redo))
  
  # Get shapefiles
  sf_object <- getSfObjects(best_seas)
  
  best_model <- best_seas$variant[1]
  best_run_level <- best_seas$run_level[1]
  
  lambdas <- getLambdas(country = country,
                        identifier = identifier,
                        run_level = best_run_level, 
                        variant = best_model,
                        time_left = time_left,
                        time_right = time_right,
                        redo = redo)
  
  best_seas <- best_seas %>% 
    inner_join(lambdas)
  
  p_betas <- best_seas %>% 
    ggplot(aes(x = month, group = gid)) +
    geom_hline(aes(yintercept = 0), lty = 2, size = .1) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = group), alpha = 1/(nrow(sf_object)*1.2)) +
    geom_line(aes(y = mean, color = group), alpha = 1/(nrow(sf_object)/10)) +
    facet_grid(group~.) +
    scale_color_manual(values = getGroupingColors()) +
    scale_fill_manual(values = getGroupingColors()) +
    scale_x_continuous(breaks = seq(1, 12, by = 2), labels = month.abb[seq(1, 12, by = 2)]) +
    guides(fill = "none", color = "none") +
    theme_bw() +
    labs(y = "seaonality coefficient")
  
  p_lambdas <- plotGroupingWeights(country = country,
                                   identifier = identifier,
                                   run_level = best_run_level, 
                                   variant = best_model,
                                   time_left = time_left,
                                   time_right = time_right,
                                   sf_object = sf_object)
  
  p_group <- cowplot::plot_grid(p_lambdas, 
                                p_betas,
                                nrow = 1, 
                                labels = "auto", 
                                rel_widths = c(1, 1.1))
  return(p_group)
  
}


#' @title Plot sample chains
#'
#' @description Plots the samples
#'
#' @param sample_df
#'
#' @return a ggplot object
#' 
plotChains <- function(sample_df) {
  ggplot2::ggplot(sample_df, ggplot2::aes(x = iter, y = sample, 
                                          group = param, color = factor(chain))) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(identifier ~ run_level + variant) +
    ggplot2::theme_bw()
}

#' @title Plot sample chains
#'
#' @description Plots the samples
#'
#' @param sample_df
#'
#' @return a ggplot object
#' 
plotSampleDensities <- function(sample_df) {
  ggplot2::ggplot(sample_df, ggplot2::aes(x = sample, color = param)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(identifier ~ run_level + variant) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = "none")
}

#' @title plot all model status
#'
#' @description makes a tile plot of the model status
#'
#' @param param
#'
#' @return return
plotAllBestModel <- function(countries = "all",
                             identifier,
                             run_level,
                             time_left = NULL,
                             time_right = NULL,
                             redo_single = F,
                             redo = F,
                             save_fig = F,
                             fig_dir = getFigDir()) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Get all best seasonality coefficients
  best_model <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = run_level,
    models = "null",
    fun = getBestModel,
    fun_name = "best_model",
    fun_opts = list(as_df = T,
                    model_ics = NULL),
    redo = list(main = redo),
    error_handling = "stop"
  ) 
  
  p_best_model <- best_model %>% 
    ggplot2::ggplot(ggplot2::aes(x = factor(1), 
                                 y = country, 
                                 fill = best_model)) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = getModelColors())
  
  
  if (save_fig) {
    ggplot2::ggsave(p_best_model, 
                    filename = paste0(getFigDir(), "/model_status_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 7,
                    height = 6,
                    dpi = 300)
  }
  p_best_model
}

#' @title plot all best models
#'
#' @description makes a tile plot of the model status
#'
#' @param param
#'
#' @return return
plotAllModelStatus <- function(countries = "all",
                               identifier,
                               run_level,
                               time_left = NULL,
                               time_right = NULL,
                               redo = F,
                               save_fig = F,
                               fig_dir = getFigDir()) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Get all best seasonality coefficients
  model_status <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = run_level,
    models = "null",
    fun = getModelRunStatus,
    fun_name = "model_run_status",
    fun_opts = list(redo_logs = F),
    redo = redo,
    error_handling = "stop",
  ) 
  
  p_model_status <- model_status %>% 
    dplyr::mutate(fit_avail = ifelse(is.na(fit_avail), "missing", "available")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = model, y = country, fill = fit_avail)) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::scale_fill_manual(values = c("lightblue", "darkgray"))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_model_status, 
                    filename = paste0(getFigDir(), "/model_status_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 7,
                    height = 6,
                    dpi = 300)
  }
  p_model_status
}




#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeas <- function(countries = "all",
                        identifier,
                        run_level,
                        time_left = NULL,
                        time_right = NULL,
                        best_seas = NULL,
                        redo_single = F,
                        redo_data = F,
                        redo = F,
                        save_fig = F,
                        fig_dir = getFigDir(),
                        verbose = T,
                        beta_thresh = 3.5) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas)) {
    best_seas <- getBestSeasByGroup(countries = countries,
                                    identifier = identifier,
                                    run_level = run_level,
                                    redo = list(main = redo))
    
  } else {
    all_countries <- unique(best_seas$country)
  }
  
  best_seas <- tidyBestSeasGrp(best_seas, all_countries) 
  
  # Get Country shapefiles to extract latitude
  latitudes <- getLatitude(countries = "all", redo = F)
  
  # Add latitude to reorder
  best_seas <- best_seas %>% 
    dplyr::inner_join(latitudes) %>% 
    dplyr::mutate(country = reorder(country, -latitude*1.5 + longitude))
  
  # Adjust color scales
  best_seas$mean[best_seas$mean > beta_thresh] <- beta_thresh
  best_seas$mean[best_seas$mean < -beta_thresh] <- -beta_thresh
  
  # Plot heatmap of seasonal coefs
  p_heatmap <- best_seas %>%
    ggplot(aes(x = param_fam, y = 1, fill = mean)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(country ~ ., as.table = T, switch = "y", scales = "free") +
    scale_fill_gradient2(low = "blue", high = "red", na.value = "lightgray") +
    scale_y_continuous(breaks = 0.5:11.5, labels = month.abb, limits = c(0, 12)) +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      panel.spacing = unit(0, "pt"),
      panel.border = element_rect(color = "black", size = .3)
    ) +
    labs(x = "", y = "") +
    guides(fill = guide_colorbar("Seasonality coefficient"))
  
  if (save_fig) {
    ggplot2::ggsave(p_heatmap, 
                    filename = paste0(getFigDir(), "/seas_coef_heatmap_",
                                      stringr::str_c(countries, collapse = "-"),
                                      "_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 7,
                    height = 8,
                    dpi = 300)
  }
  p_heatmap
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasMapsUncertainty <- function(countries = "all",
                                       identifier,
                                       run_level,
                                       time_left = NULL,
                                       time_right = NULL,
                                       best_seas = NULL,
                                       sf_objects = NULL,
                                       africa.sf_lev0 = NULL,
                                       redo_single = F,
                                       redo_data = F,
                                       redo = F,
                                       save_fig = F,
                                       fig_dir = getFigDir(),
                                       verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas)) {
    best_seas <- getBestSeas(countries = countries,
                             identifier = identifier,
                             run_level = run_level,
                             redo = list(main = redo),
                             redo_single = redo_single)
  }
  
  if (is.null(sf_objects)) {
    sf_objects <- getSfObjects(best_seas)
  }
  
  if (is.null(africa.sf_lev0)) {
    africa.sf_lev0 <- getCountrySf(all_countries)
  }
  
  best_seas <- best_seas %>% 
    dplyr::select(gid, month, mean, q2.5, q97.5) %>% 
    dplyr::inner_join(sf_objects, .) %>% 
    dplyr::mutate(width = abs(q97.5-q2.5))
  
  month_dict <- month.abb
  names(month_dict) <- 1:12
  
  
  lakes_sf <- getLakes()
  
  # Get gids for countries for which seasonality was not selected
  sf_objects <- dplyr::bind_rows(
    purrr::map_df(getExcludedCountries(), ~getCountrySf(.)) %>% dplyr::mutate(seas_index = -1),
    purrr::map_df(getNoSeasCountries(), ~getCountrySf(.)) %>% dplyr::mutate(seas_index = 0)
  )
  
  p_seas_maps <- best_seas %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), fill = "gray90", color = "white", size = .3) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), fill = "gray75", color = "white", size = .3) +
    geom_sf(aes(fill = width), color = "lightgray", size = 0) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "lightgray", size = .005, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "darkgray", size = .1, alpha = 0) +
    scale_fill_viridis_c() +
    ggplot2::facet_wrap(~month, ncol = 3, labeller = labeller(month = month_dict)) +
    ggthemes::theme_map() +
    theme(strip.text = element_text(margin = margin(.05, 0, .05, 0, unit = "cm"))) +
    theme(legend.position = "bottom",
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.key.height = unit(.2, unit = "cm"),
          plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = guide_colorbar("Width of 95% CrI",
                                 title.position = "top"))
  
  if (save_fig) {
    ggplot2::ggsave(p_seas_maps, 
                    filename = paste0(getFigDir(), "/seas_coef_maps_uncertainty",
                                      stringr::str_c(countries, collapse = "-"),
                                      "_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 6,
                    height = 7,
                    dpi = 500)
  }
  p_seas_maps
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasMaps <- function(countries = "all",
                            identifier,
                            run_level,
                            time_left = NULL,
                            time_right = NULL,
                            best_seas = NULL,
                            sf_objects = NULL,
                            africa.sf_lev0 = NULL,
                            redo = F,
                            save_fig = F,
                            fig_dir = getFigDir(),
                            verbose = T,
                            beta_thresh = 2.5) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas)) {
    best_seas <- getBestSeas(countries = countries,
                             identifier = identifier,
                             run_level = run_level,
                             redo = list(main = redo))
  }
  
  if (is.null(sf_objects)) {
    sf_objects <- getSfObjects(best_seas)
  }
  
  
  if (is.null(africa.sf_lev0)) {
    africa.sf_lev0 <- getCountrySf(all_countries)
  }
  
  best_seas <- best_seas %>% 
    dplyr::select(gid, month, mean) %>% 
    dplyr::inner_join(sf_objects, .)
  
  month_dict <- month.abb
  names(month_dict) <- 1:12
  
  # Adjust color scales
  best_seas$mean[best_seas$mean > beta_thresh] <- beta_thresh
  best_seas$mean[best_seas$mean < -beta_thresh] <- -beta_thresh
  
  # Lake data from https://datacatalog.worldbank.org/dataset/africa-water-bodies-2015
  lakes_sf <- getLakes()
  
  p_seas_maps <- best_seas %>% 
    dplyr::filter(!is.na(mean)) %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), fill = "gray90", size = 0) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), fill = "gray75", size = 0) +
    geom_sf(aes(fill = exp(mean)), color = "lightgray", size = 0) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "lightgray", size = .005, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "gray40", size = .1, alpha = 0) +
    scale_fill_gradient2(low = "blue", 
                         high = "red", 
                         trans = "log",
                         breaks = c(.1, 1, 10),
                         labels = formatC(c(.1, 1, 10)),
                         na.value = "gray75") +
    ggplot2::facet_wrap(~month, ncol = 3, labeller = labeller(month = month_dict)) +
    ggthemes::theme_map() +
    theme(strip.text = element_text(margin = margin(.05, 0, .05, 0, unit = "cm")))
  
  if (save_fig) {
    ggplot2::ggsave(p_seas_maps, 
                    filename = paste0(getFigDir(), "/seas_coef_maps_",
                                      stringr::str_c(countries, collapse = "-"),
                                      "_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 6,
                    height = 7,
                    dpi = 500)
  }
  p_seas_maps
}

#' @title Plot seasonality combined
#'
#' @description Plots the heatmap along with the monthly maps
#'
#' @param param
#'
#' @return ggplot object
#' 
plotAllSeasCombined <- function(p_seas = NULL,
                                p_maps = NULL,
                                countries = "all",
                                run_level = "best",
                                identifier = "mean_annual_incidence",
                                time_left = NULL,
                                time_right = NULL,
                                redo = F,
                                redo_single = F,
                                save_fig = F,
                                verbose = F) {
  
  if (is.null(p_seas)) {
    # Plot seasonality coefficients
    p_seas <- plotAllSeas(countries = countries,
                          identifier = identifier,
                          run_level = run_level,
                          time_left = time_left,
                          time_right = time_right,
                          best_seas = NULL,
                          redo = redo,
                          save_fig = F)
  }
  
  if (is.null(p_maps)) {
    # Plot of seasonality coefficients as a map
    p_maps <- plotAllSeasMaps(countries = countries,
                              identifier = identifier,
                              run_level = run_level,
                              time_left = time_left,
                              time_right = time_right,
                              redo = redo,
                              redo_single = redo_single,
                              redo_data = F,
                              save_fig = F,
                              verbose = verbose)
    
  }
  
  
  # Combine plots
  p_seas_combined <- cowplot::plot_grid(
    p_seas +
      theme(legend.position = "bottom",
            plot.margin = unit(c(3, 2, 1, 1), "lines")), 
    p_maps +
      guides(fill = "none") +
      theme(legend.position = "bottom",
            plot.margin = unit(c(2, 2, 2, 2), "lines")),
    labels = "auto",
    rel_widths = c(.8, 1)
  )
  
  if (save_fig) {
    ggsave(p_seas_combined,   
           filename = paste0(getFigDir(), "/seas_coef_combine_all",
                             "_rlbest_",
                             identifier,
                             "_",
                             makePeriod(NULL, NULL),
                             ".png"),
           width = 10,
           height = 8,
           dpi = 300)
  }
  
  p_seas_combined
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllClimMaps <- function(countries = "all",
                            covariate = "Precipitation",
                            gadm_lev = 2,
                            save_fig = F,
                            fig_dir = getFigDir(),
                            verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Gadm lev2 objects
  sf_objects <- purrr::map_df(all_countries,
                              ~ getShapefile(., gadm_lev)) %>% 
    rmapshaper::ms_simplify(keep = .1)
  
  africa.sf_lev0 <- getCountrySf(all_countries)
  
  # Pull covariates
  clim <- getClimatology(country = all_countries, 
                         admin_level = gadm_lev, 
                         covar = covariate) %>% 
    dplyr::group_by(month, gid, covar) %>% 
    dplyr::arrange(gid, month, covar, source) %>% 
    dplyr::slice(1) %>% 
    dplyr::select(-source)
  
  # Join with sf
  clim <- dplyr::inner_join(sf_objects, clim)
  
  # Month dictionnary for faceting
  month_dict <- month.abb
  names(month_dict) <- 1:12
  
  p_clim_maps <- ggplot(clim) +
    geom_sf(aes(fill = value), color = "lightgray", size = .01) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .05, alpha = 0) +
    scale_fill_gradient(low = getClimColors()[[covariate]]$low, high = getClimColors()[[covariate]]$high) +
    ggplot2::facet_wrap(~month, ncol = 3, labeller = labeller(month = month_dict)) +
    ggthemes::theme_map() +
    theme(legend.position = "bottom")
  
  if (save_fig) {
    ggplot2::ggsave(p_clim_maps, 
                    filename = makeFigFilename(prefix = paste0("seas_clim_maps_", covariate),
                                               countries = countries,
                                               run_level = "",
                                               identifier = ""),
                    width = 6,
                    height = 7,
                    dpi = 500)
    
  }
  p_clim_maps
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllPrecipPeakLag <- function(countries = "all",
                                 identifier,
                                 run_level,
                                 redo = F,
                                 redo_single = redo,
                                 save_fig = F,
                                 verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  peak_month <- getAllPeakMonths(countries = countries,
                                 identifier = identifier,
                                 run_level = run_level,
                                 redo = list(main = redo_single))
  
  
  # Get sf objects
  sf_objects <- peak_month %>% 
    dplyr::group_by(country) %>% 
    dplyr::group_map(function(x, y) {
      # admin_lev <- ifelse(!stringr::str_detect(x$variant[1], "mixture") | is.na(x$variant[1]), 0, x$run_level[1])
      getShapefile(y$country[1], x$run_level[1])
    }) %>% 
    dplyr::bind_rows() %>% 
    rmapshaper::ms_simplify(keep = .1)
  
  africa.sf_lev0 <- getCountrySf(all_countries)
  
  # Pull covariates
  clim <- peak_month %>% 
    dplyr::group_by(country) %>% 
    dplyr::group_map(function(x, y) {
      getClimatology(country = y$country[1], 
                     admin_level = x$run_level[1], 
                     covar = "Precipitation") %>% 
        dplyr::group_by(month, gid, covar) %>% 
        dplyr::arrange(gid, month, covar, source) %>% 
        dplyr::slice(1) %>% 
        dplyr::select(-source)
    }) %>% 
    dplyr::bind_rows()
  
  precip_peak <- clim %>% 
    dplyr::group_by(gid, gadm_lev) %>% 
    dplyr::summarise(peak_precip_month = which.max(value),
                     peak_precip_val = max(value))
  
  # Join with sf
  peak_data <- dplyr::inner_join(sf_objects, peak_month) %>% 
    dplyr::inner_join(precip_peak)
  
  # Compute lag between precip and seas
  peak_data <- peak_data %>% 
    dplyr::mutate(
      lag_precip = peak_month - peak_precip_month,
      lag_precip = case_when(lag_precip < -6  ~ lag_precip + 12,
                             lag_precip > 6 ~ lag_precip - 12, 
                             T ~ lag_precip)
    )
  
  # Month dictionnary for faceting
  month_dict <- month.abb
  names(month_dict) <- 1:12
  
  p_lag_maps <- ggplot(peak_data) +
    geom_sf(data = africa.sf_lev0, fill = "gray", color = "lightgray", size = .2, alpha = 1) +
    geom_sf(aes(fill = lag_precip), color = "lightgray", size = .1) +
    geom_sf(data = africa.sf_lev0, color = "darkgray", size = .3, alpha = 0) +
    scale_fill_gradient2() +
    ggthemes::theme_map() +
    theme(legend.position = "bottom")
  
  if (save_fig) {
    ggplot2::ggsave(p_lag_maps, 
                    filename = makeFigFilename(prefix = "precip_lag",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 6,
                    height = 5,
                    dpi = 500)
    
  }
  p_lag_maps
}




#' @title plot all data availability
#'
#' @description makes a plot of the availability of data across months and 
#' for each month
#'
#' @param param
#'
#' @return return
plotAllDataMaps <- function(countries = "all",
                            identifier,
                            run_level,
                            time_left = NULL,
                            time_right = NULL,
                            redo_single = F,
                            redo_data = F,
                            redo = F,
                            save_fig = F,
                            verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Get all best seasonality coefficients
  run_data <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = 1,
    models = "best",
    fun = getRunDataWrapper,
    fun_name = "run_data",
    fun_opts = list(redo = redo_single),
    redo = redo
  )
  
  data_stats <- run_data %>% 
    mutate(n_months = map_dbl(months, length),
           month = lubridate::month(month_left),
           n_months_cat = cut(n_months, c(0, 3, 6, 11, 12, Inf))) %>% 
    group_by(gid, n_months_cat) %>% 
    summarise(year_min = min(year),
              year_max = max(year),
              n_years = length(unique(year)))
  
  sf_objects <- map_df(all_countries, function(x) {
    map_df(0:2, function(y) {
      getShapefile(x, y)
    })
  }) %>% 
    rmapshaper::ms_simplify(keep = .1)
  
  africa.sf_lev0 <- getCountrySf(all_countries)
  africa.sf_lev0_facet <- map_df(0:2, function(x) {
    africa.sf_lev0 %>% 
      dplyr::mutate(gadm_lev = x)
  })
  africa.sf_excl <- africa.sf_lev0_facet %>% 
    dplyr::filter(country %in% getExcludedCountries()) 
  
  # africa.sf_lev0 <- purrr::map_df(0:2, ~mutate(africa.sf_lev0, gadm_lev = .))
  
  data_stats <- data_stats %>% 
    dplyr::inner_join(sf_objects, .)
  
  lakes_sf <- getLakes()
  
  obs_length_dict <- c("(0,3]" = "< 3 months",
                       "(3,6]" = "3-6 months",
                       "(6,11]" = "6-11 months",
                       "(11,12]" = "1 year",
                       "(12,Inf]" = "> 1 year")
  
  
  gadm_lev_dict <- c("0" = "GADM0",
                     "1" = "GADM1",
                     "2" = "GADM2")
  
  p_data_maps <- data_stats %>% 
    ggplot() +
    geom_sf(data = africa.sf_excl, fill = "gray75", size = 0) +
    geom_sf(aes(fill = n_years), color = "lightgray", size = .015) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "white", size = .015, alpha = 0) +
    geom_sf(data = africa.sf_lev0_facet, color = "black", size = .2, alpha = 0) +
    scale_fill_viridis_c() +
    ggplot2::facet_grid(n_months_cat ~ gadm_lev, 
                        labeller = labeller(n_months_cat = obs_length_dict,
                                            gadm_lev = gadm_lev_dict)) +
    ggthemes::theme_map() +
    guides(fill = guide_colorbar("Number of years", title.position = "top")) +
    theme(legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = "white"))
  
  if (save_fig) {
    ggplot2::ggsave(p_data_maps, 
                    filename = makeFigFilename(prefix = "data_maps_",
                                               countries = countries,
                                               run_level = "",
                                               identifier = identifier),
                    width = 7,
                    height = 9,
                    dpi = 300)
  }
  p_data_maps
}

plotDefinitionComp <- function(countries = "all",
                               redo = F,
                               redo_single = F,
                               save_fig = T) {
  
  # Get model results for all run levels
  all_model_res <- runAll(countries = countries,
                          models = "best",
                          run_levels = "best",
                          fun = getBestModel,
                          fun_name = "best_model",
                          fun_opts = list(as_df = T,
                                          model_ics = NULL,
                                          redo = redo_single),
                          redo = redo, 
                          error_handling = "remove")
  
  all_model_res <- all_model_res %>% 
    select(-nmodels, -models) %>% 
    pivot_wider(names_from = "identifier",
                values_from = "best_model") %>% 
    pivot_longer(cols = c("occurrence", "cases"),
                 names_to = "identifier") %>% 
    mutate(mean_annual_incidence = factor(mean_annual_incidence, levels = getAllModels()),
           value = factor(value, levels = getAllModels())) %>% 
    filter(!(country %in% getExcludedCountries()))
  
  p_compare <- all_model_res %>% 
    filter(!is.na(value)) %>% 
    count(run_level, mean_annual_incidence, identifier, value) %>% 
    group_by(identifier, run_level) %>% 
    mutate(frac = n/sum(n)) %>% 
    ggplot(aes(x = value, y = mean_annual_incidence)) +
    geom_tile(aes(fill = frac)) +
    theme_bw() +
    facet_grid( ~ identifier,
                labeller = labeller(identifier = c("cases" = "10 cases or more",
                                                   "occurrence" = "1 case or more"))) +
    scale_fill_gradient(low = "gray90", high = "darkgreen") +
    labs(x = "Selected model in alternative cholera occurrence definition", 
         y = "Selected model in main cholera occurrence definition") +
    guides(fill = guide_colorbar("Fraction of\ncountries"))
  
  if (save_fig) {
    
    ggplot2::ggsave(p_compare, 
                    filename = makeFigFilename(
                      prefix = "cholera_def_comp",
                      countries = countries,
                      run_level = "best",
                      identifier = "all"),
                    width = 10,
                    height = 5,
                    dpi = 300)
    
  }
}

#' @title Plot mean 2010-2016 annual cholera incidence rate
#'
#' @description Plot mean 2010-2016 annual cholera incidence rate
#'
#' @param country
#' @param gadm_data
#'
#' @return a plot with the population raster
#' 
plotMeanRate <- function(countries,
                         gadm_data = getGadmData()) {
  
  if (countries == "all") {
    countries <- getAllCountries()
  } 
  
  sf_objects <- purrr::map_df(0:2, function(x) {
    purrr::map_df(countries, function(y) {
      getShapefile(country = y, admin_level = x)
    })
  }) %>% 
    rmapshaper::ms_simplify(keep = .1)
  
  data <- gadm_data %>% 
    dplyr::filter(cntry %in% countries) %>% 
    dplyr::inner_join(sf_objects, .)
  
  ggplot(data, aes(fill = meanrate)) +
    geom_sf(size = .05, color = "white") +
    geom_sf(inherit.aes = F,
            data = purrr::map_df(0:2,
                                 function(x) {
                                   sf_objects %>%
                                     dplyr::filter(gadm_lev == 0) %>% 
                                     dplyr::mutate(gadm_lev = x)}),
            size = .5, color = "white", alpha = 0) +
    facet_wrap(~gadm_lev, labeller = label_both) +
    scale_fill_viridis_c(trans = "log10") +
    ggthemes::theme_map() +
    theme(legend.position = "bottom")
}



#' @title Plot population
#'
#' @description Plot worldpop population
#'
#' @param country
#' @param gadm_data
#'
#' @return a plot with the population raster
#' 
plotPopulation <- function(countries,
                           gadm_data = getGadmData()) {
  
  if (countries == "all") {
    countries <- getAllCountries()
  } 
  
  sf_objects <- purrr::map_df(0:2, function(x) {
    purrr::map_df(countries, function(y) {
      getShapefile(country = y, admin_level = x)
    })
  })
  
  data <- gadm_data %>% 
    dplyr::filter(cntry %in% countries) %>% 
    dplyr::inner_join(sf_objects, .)
  
  map(1:2, function(x) {
    data %>% 
      dplyr::filter(gadm_lev == x) %>% 
      ggplot(aes(fill = popsum)) +
      geom_sf(size = .05, color = "white") +
      geom_sf(inherit.aes = F,
              data = sf_objects %>%
                dplyr::filter(gadm_lev == 0) %>% 
                dplyr::mutate(gadm_lev = x), 
              size = .5, color = "white", alpha = 0) +
      facet_wrap(~gadm_lev, labeller = label_both) +
      scale_fill_viridis_c(trans = "log10") +
      ggthemes::theme_map() +
      theme(legend.position = "bottom")}
  ) %>% 
    cowplot::plot_grid(plotlist = .)
  
}

#' @title plot all seasonality index
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasIndexScatter <- function(countries = "all",
                                    identifier,
                                    run_level,
                                    time_left = NULL,
                                    time_right = NULL,
                                    redo_single = F,
                                    redo = F,
                                    redo_single_probs = F,
                                    save_fig = F,
                                    fig_dir = getFigDir()) {
  
  # Correlations file
  seasind_file <- makeStdResName(country = countries,
                                 run_level = run_level,
                                 identifier = identifier,
                                 model = "best",
                                 gadm_lev = "auto",
                                 suffix = "seas_index",
                                 file_type = "rds")
  
  # Get seasonality index results
  if (file.exists(seasind_file) & ! redo) {
    seas_index <- readRDS(seasind_file)
  } else {
    stop("Need to pre-compute seasonality index results")
  }
  
  # Get population data
  seas_index_long <- seas_index %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-geometry) %>% 
    dplyr::inner_join(getGadmData() %>% 
                        dplyr::mutate(meanrate = cholera_rates_2010_2016,
                                      meanrate = meanrate * 1e5,
                                      popsum = popsum_2020)) %>% 
    dplyr::select(gid, country, seas_index, area, popsum, meanrate) %>% 
    dplyr::mutate(pop_density = popsum/area) %>%
    dplyr::rename(pop = popsum) %>% 
    tidyr::pivot_longer(cols = c("area", "pop", "pop_density", "meanrate"),
                        names_to = "var")
  
  var_dict <- c("area" = "Area [sqkm]",
                "meanrate" = "Mean annual incidence [cases/100'000/year]",
                "pop" = "Population",
                "pop_density" = "Population density [people/sqkm]")
  
  p_scatter <- seas_index_long %>% 
    ggplot(aes(x = log10(value), y = seas_index)) +
    geom_point(alpha = .1) +
    geom_smooth() +
    facet_wrap(~var, scales = "free_x", labeller = labeller(var = var_dict)) +
    theme_bw() +
    labs(x = "Log10 of value", y = "Seasonality index")
  
  if (save_fig) {
    ggplot2::ggsave(p_scatter, 
                    filename = makeFigFilename(prefix = "seas_index_scatter",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 8,
                    height = 6,
                    dpi = 300)
  }
  p_scatter
}

#' @title plot all seasonality index
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasInd2 <- function(countries = "all",
                            identifier,
                            run_level = "best",
                            time_left = NULL,
                            time_right = NULL,
                            redo = F,
                            save_fig = F,
                            fig_dir = getFigDir()) {
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Country shapefiles
  africa.sf_lev0 <- getCountrySf(all_countries)
  
  seas_index <- computeAllSeasIndex(countries = countries,
                                    run_levels = run_level,
                                    identifier,
                                    redo = list(main = redo),
                                    verbose = F,
                                    error_handling = "remove") %>% 
    rmapshaper::ms_simplify(keep = 0.2)
  
  seas_cat_dict <- c(
    "(-2,-1]" = "excluded",
    "(-1,0]" = "non-seasonal",
    "(0,0.3]" = "< 0.3",
    "(0.3,0.5]" = "0.3 - 0.5",
    "(0.5,0.7]" = "0.5 - 0.7",
    "(0.7,1]" = "> 0.7"
  )
  
  # Get gids for countries for which seasonality was not selected
  sf_objects <- dplyr::bind_rows(
    purrr::map_df(getExcludedCountries(), ~getCountrySf(.)) %>% dplyr::mutate(seas_index = -1),
    purrr::map_df(getNoSeasCountries(), ~getCountrySf(.)) %>% dplyr::mutate(seas_index = 0)
  ) 
  
  # Get Gadm data
  pop_data <- getGadmData() %>%
    rename(meanrate = cholera_rates_2010_2016) %>% 
    dplyr::inner_join(seas_index %>% 
                        dplyr::as_tibble() %>% 
                        dplyr::select(gid, seas_index)) %>% 
    bind_rows(., getGadmData() %>%
                rename(meanrate = cholera_rates_2010_2016) %>% 
                inner_join(sf_objects %>% 
                             dplyr::as_tibble() %>% 
                             dplyr::select(gid, seas_index))) %>% 
    dplyr::mutate(meanrate = ifelse(is.na(meanrate), 0, meanrate)) %>% 
    dplyr::mutate(rate_cat = cut(meanrate, getRateCuts())) %>% 
    dplyr::mutate(seas_cat = cut(seas_index, c(-2, -1, 0, .3, .5, .7, 1)) %>% as.character()) %>% 
    {x <-.
    rbind(x %>% dplyr::mutate(var = "Incidence rate"), 
          x %>% dplyr::mutate(rate_cat = "Total", var = "Total"))
    } %>% 
    dplyr::group_by(seas_cat, rate_cat, var) %>% 
    dplyr::summarise(pop = sum(popsum_2020, na.rm = T)) %>%
    dplyr::group_by(rate_cat, var) %>% 
    dplyr::mutate(frac = pop/sum(pop)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(rate_cat = c(getRateDict(), "Total" = "Total")[rate_cat],
                  rate_cat = factor(rate_cat, levels = c(getRateDict(), "Total")),
                  seas_cat = seas_cat_dict[seas_cat],
                  seas_cat = factor(seas_cat, levels = rev(seas_cat_dict)))
  
  # Save output
  saveRDS(pop_data,
          file = makeStdResName(country = countries, 
                                run_level = "best", 
                                identifier = identifier, 
                                time_left = time_left,
                                time_right = time_right, 
                                gadm_lev = "best",
                                model = "best", 
                                suffix = "rate_fraction_seas", 
                                file_type = "rds"))
  
  p_frac <- ggplot(pop_data, aes(x = rate_cat, 
                                 y = frac,
                                 fill = seas_cat)) +
    geom_bar(stat = "identity") +
    facet_grid(~var, scales = "free", space = "free") +
    labs(x = "Cases/100'000/year", 
         y = "Fraction of\npopulation") +
    scale_fill_manual(values = rev(c("gray90", "gray75", getSeasCatColors()))) +
    guides(fill = guide_legend("Seasonality index")) +
    theme_bw()
  
  # Get lakes of SSA
  lakes_sf <- getLakes()
  
  # Plot seasonality indices
  p_seasindex <- seas_index %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), fill = "gray90", color = "white", size = .3) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), fill = "gray75", color = "white", size = .3) +
    geom_sf(aes(fill = seas_index), color = "white", size = .025) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(aes(color = seas_index>.5, size  = seas_index>.5), alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .2, alpha = 0) +
    scale_fill_gradient2(midpoint = .5, 
                         low = "white", 
                         mid = "red",
                         high = "darkred") +
    scale_color_manual(values = c("white", "black")) +
    scale_size_manual(values = c( .035, .15)) +
    guides(color = "none", size = "none", 
           fill = guide_colorbar("Fraction of occurrence\nin three months around peak",
                                 title.position = "top")) +
    ggthemes::theme_map() +
    theme(legend.position = c(.05, .05),
          legend.direction = "horizontal")
  
  # p_seasindex <- addSpatialAnnotations(p_seasindex)
  p_seasindex2 <- cowplot::plot_grid(
    p_seasindex +
      theme(plot.margin = margin(0, 0, 0, 0, unit = "lines"),
            legend.key.height = unit(.18, unit = "cm")),
    p_frac +
      theme(plot.margin = margin(0, 1, 1, 1, unit = "lines"),
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 9),
            plot.background = element_blank(),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 9),
            legend.key.size = unit(.3, units = "cm"),
            strip.text = element_text(size = 8, 
                                      margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"))),
    ncol = 1,
    rel_heights = c(1, .3))
  
  if (save_fig) {
    ggplot2::ggsave(p_seasindex2, 
                    filename = makeFigFilename(prefix = "seas_index_map",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 8,
                    height = 8,
                    dpi = 300)
  }
  
  p_seasindex2
}

getSeasIndSim <- function(country,
                          run_level,
                          model = "best",
                          identifier,
                          time_left = NULL,
                          time_right = NULL,
                          gadm_lev = setGadmLev(country, run_level),
                          redo = getRedoDefaults(),
                          verbose = F,
                          ...) {
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo$data,
                                 verbose = verbose)  
  }
  
  gadm_lev = setGadmLev(country, run_level)
  
  if (verbose) {
    cat("Getting betas for", country, run_level, identifier, "\n")
  }
  
  if (model == "best") {
    model <- getBestModel(country = country,
                          run_level = run_level,
                          identifier = identifier,
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = gadm_lev,
                          redo = redo$best_model)
  }  
  
  # First try to get only betas simulations
  res_file <- try(getResFiles(country = country,
                              run_level = run_level, 
                              identifier = identifier, 
                              time_left = time_left,
                              time_right = time_right,
                              gadm_lev = gadm_lev,
                              what = "seas",
                              paths = getResPaths(what = "outputs_sim"),
                              model = model))
  sims <- readRDS(res_file)
  
  stats <- sims %>% 
    group_by(gid) %>% 
    summarise(seas_index_mean = mean(seas_index),
              seas_index_q025 = quantile(seas_index, .025),
              seas_index_q975 = quantile(seas_index, .975),
              seas_index_postprob = sum(seas_index > .25)/length(seas_index)) %>% 
    bind_cols(parseResFilename(res_file))
  
  return(stats)
}

#' @title plot all seasonality index
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasIndexUncertainty <- function(countries = "all",
                                        identifier,
                                        run_level,
                                        nsim = 500,
                                        time_left = NULL,
                                        time_right = NULL,
                                        redo = F,
                                        redo_single = F,
                                        save_fig = F,
                                        fig_dir = getFigDir(),
                                        error_handling = "remove") {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Correlations file
  seasind_file <- makeStdResName(country = countries,
                                 run_level = run_level,
                                 identifier = identifier,
                                 model = "best",
                                 gadm_lev = "auto",
                                 suffix = "seas_index_uncertainty",
                                 file_type = "rds")
  
  # Country shapefiles
  africa.sf_lev0 <- getCountrySf(all_countries)
  
  if (file.exists(seasind_file) & !redo) {
    all_seas_stats <- readRDS(seasind_file)
  } else {
    
    all_seas_stats <- runAll(countries = countries,
                             run_levels = "best",
                             identifiers = identifier,
                             models = "best",
                             fun = getSeasIndSim,
                             fun_name = "seas_ind_sim",
                             error_handling = error_handling,
                             redo = redo,
                             verbose = verbose)
    
    # Get sf objects
    sf_objects <- getSfObjects(all_seas_stats %>% rename(mean = seas_index_mean))
    
    all_seas_stats <- sf_objects %>% 
      dplyr::inner_join(all_seas_stats)
    
    saveRDS(all_seas_stats, file = seasind_file)
  }
  
  # Get lakes of SSA
  lakes_sf <- getLakes()
  
  # Plot seasonality indices
  p_seasindex <- all_seas_stats %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0, fill = "lightgray", color = "white", size = .3) +
    geom_sf(aes(fill = seas_index_postprob), color = "white", size = .035) +
    geom_sf(data = all_seas_stats %>% 
              dplyr::filter(seas_index_postprob > .95),
            color = "black", size = .1, alpha = 0) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .2, alpha = 0) +
    scale_fill_viridis_c() +
    guides(color = "none", size = "none", 
           fill = guide_colorbar("Posterior probability of seasonality\nstronger than the null",
                                 title.position = "top")) +
    ggthemes::theme_map() +
    theme(legend.position = c(.05, .05),
          legend.direction = "horizontal",
          plot.background = element_rect(fill = "white", color = "white"))
  
  if (save_fig) {
    ggplot2::ggsave(p_seasindex, 
                    filename = makeFigFilename(prefix = "seas_index_map_uncertainty",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 8,
                    height = 8,
                    dpi = 300)
  }
  
  p_seasindex
}

#' @title title
#'
#' @description description
#'
#' @param param
#'
#' @return return


#' @title plot all peak months
#'
#' @description makes a tile plot of the peak month for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllPeakMonth <- function(countries = "all",
                             identifier,
                             run_level,
                             best_seas = NULL,
                             sf_objects = NULL,
                             time_left = NULL,
                             time_right = NULL,
                             redo_single = F,
                             redo = F,
                             save_fig = F,
                             fig_dir = getFigDir()) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  peak_month <- getAllPeakMonths(countries = countries,
                                 identifier = identifier,
                                 run_level = run_level,
                                 best_seas = best_seas,
                                 redo = list(main = redo))
  
  if (is.null(sf_objects)) {
    sf_objects <- getSfObjects(peak_month %>% mutate(mean = peak_month))
  }
  
  africa.sf_lev0 <- getCountrySf(all_countries)
  lakes_sf <- getLakes()
  
  p_peak_month <- sf_objects %>% 
    dplyr::inner_join(peak_month) %>%
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), 
            fill = "gray90", color = "white", size = .3) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), 
            fill = "gray75", color = "white", size = .3) +
    geom_sf(aes(fill = peak_month, alpha = log(peak_val)), color = "white", size = .035) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "white", size = .075, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .2, alpha = 0) +
    scale_fill_gradientn(colors = getMonthGradientCols(),
                         na.value = "lightgray",
                         limits = c(1, 12)) +
    scale_alpha_continuous(range = c(.3, 1), na.value = 1) +
    ggthemes::theme_map()+
    guides(alpha = "none") +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_peak_month, 
                    filename = makeFigFilename(prefix = "peak_month_map",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 8,
                    height = 8,
                    dpi = 300)
  }
  
  p_peak_month
}


#' @title plot all model offsets
#'
#' @description makes a tile plot of the peak month for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllOffsets <- function(countries = "all",
                           identifier,
                           run_level,
                           time_left = NULL,
                           time_right = NULL,
                           redo_single = F,
                           redo = F,
                           save_fig = F,
                           fig_dir = getFigDir()) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  # Get all best seasonality coefficients
  best_offsets <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = run_level,
    models = "best",
    fun = getOffsets,
    fun_name = "best_offsets_admin",
    fun_opts = list(redo = list(main = redo_single)),
    redo = redo, 
    error_handling = "remove"
  ) 
  
  offset_df <- best_offsets %>% 
    dplyr::group_by(country, gid, variant) %>% 
    dplyr::summarise(offset_month = which.max(prob)) %>% 
    dplyr::ungroup() %>% 
    # Expand to have all combinations
    tidyr::complete(country = all_countries) %>% 
    dplyr::mutate(gid = dplyr::case_when(
      is.na(gid) ~ country,
      T ~ gid
    ))
  
  sf_objects <- best_offsets %>% 
    dplyr::group_by(country) %>% 
    dplyr::group_map(function(x, y) {
      admin_lev <- ifelse(is.na(x$gid[1]) | x$gid[1] == y$country[1], 0, x$run_level[1])
      getShapefile(y$country[1], admin_lev)
    }) %>% 
    dplyr::bind_rows() %>% 
    rmapshaper::ms_simplify(keep = .1)
  
  africa.sf_lev0 <- getCountrySf(all_countries)
  lakes_sf <- getLakes()
  
  p_offset <- sf_objects %>% 
    dplyr::inner_join(offset_df) %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), 
            fill = "gray90", color = "white", size = .3)  +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoOffsetCountries(identifier = identifier)), 
            fill = "gray75", color = "white", size = .3) +
    geom_sf(aes(fill = offset_month), color = "white", size = .035) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "white", size = .075, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .2, alpha = 0) +
    scale_fill_gradientn(colors = getMonthGradientCols(),
                         na.value = "gray75",
                         limits = c(1, 12)) +
    ggthemes::theme_map() +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_offset, 
                    filename = makeFigFilename(prefix = "offset_map",
                                               countries = countries,
                                               run_level = run_level,
                                               identifier = identifier),
                    width = 8,
                    height = 8,
                    dpi = 300)
  }
  
  p_offset
  
}

#' @title Plot data coverage
#'
#' @description Plots the data coverage in SSA
#'
#' @param data_path
#'
#' @return return
#' 
plotDataCoverage <- function(data_path = getMonthlyDataPath(),
                             save_fig = F,
                             fig_dir = getFigDir()) {
  
  all_data <- readr::read_csv(data_path, 
                              col_types = cols(), 
                              progress = F) %>% 
    mutate(type = case_when(period_type %in% c("day", "week",
                                               "biweek", "Oweek", "month") ~ "sub-monthly",
                            period_type == "year" ~ "year", 
                            T ~ "other")) %>% 
    mutate(cholera = pmax(suspected_cases, confirmed_cases)) %>% 
    dplyr::filter(country %in% getAllCountries() &
                    ((in_upper & in_lower) | (!in_lower & cholera > 0) | (!in_upper & cholera == 0)))
  
  
  # Figure of data coverage
  all_sf_objects <- map(
    dir("data/", pattern = "shp", full.names = T) %>% 
      str_subset("gadm_lev") ,
    ~ st_read(.) %>% 
      select(id, contains("gid"))
  ) %>% 
    bind_rows() %>% 
    mutate(gid_0 = as.character(gid_0),
           gid_1 = as.character(gid_1),
           gid_2 = as.character(gid_2),
           gid_0 = case_when(!is.na(gid_1) | !is.na(gid_2) ~ NA_character_,
                             T ~ gid_0),
           gid_1 = case_when(!is.na(gid_2) ~ NA_character_,
                             T ~ gid_1)) %>% 
    pivot_longer(cols = c("gid_0", "gid_1", "gid_2"), 
                 names_to = "gadm_lev",
                 values_to = "gid_chr") %>% 
    dplyr::filter(!is.na(gid_chr)) %>% 
    mutate(gadm_lev = str_extract(gadm_lev, "[0-2]") %>% as.numeric()) %>% 
    rename(gadm_id = id)
  
  other_shp <- map(
    dir("data/", patter = "shp", full.names = T) %>% 
      str_subset("africa") ,
    ~ st_read(.) %>% 
      select(cntry, contains("gid"))
  ) %>% 
    bind_rows() %>% 
    mutate(gid_lev0 = case_when(is.na(gid_lev1) & is.na(gid_lev2) ~ cntry,
                                T ~ cntry)) %>% 
    mutate(gid_lev0 = as.character(gid_lev0),
           gid_lev1 = as.character(gid_lev1),
           gid_lev2 = as.character(gid_lev2),
           gid_lev0 = case_when(!is.na(gid_lev1) | !is.na(gid_lev2) ~ NA_character_,
                                T ~ gid_lev0),
           gid_lev1 = case_when(!is.na(gid_lev2) ~ NA_character_,
                                T ~ gid_lev1)) %>% 
    pivot_longer(cols = c("gid_lev0", "gid_lev1", "gid_lev2"), 
                 names_to = "gadm_lev",
                 values_to = "gid_chr") %>% 
    mutate(gadm_lev = str_extract(gadm_lev, "[0-2]") %>% as.numeric()) %>% 
    dplyr::filter(!is.na(gid_chr))
  
  all_sf_objects <- all_sf_objects %>% 
    as_tibble() %>% 
    select(-geometry) %>% 
    inner_join(other_shp)
  
  # gadm units
  distinct_admin_units <- all_data %>% 
    dplyr::filter(type == "sub-monthly") %>% 
    mutate(`observation length` = "\nmonthly or sub-monthly") %>% 
    distinct(gadm_id, gadm_lev, `observation length`) %>% 
    bind_rows(all_data %>% 
                mutate(`observation length` = "\nall") %>% 
                distinct(gadm_id, gadm_lev, `observation length`))
  
  distinct_units_sf <- distinct_admin_units %>% 
    inner_join(all_sf_objects)
  
  p_coverage <- distinct_units_sf %>% 
    mutate(`Administrative level` = gadm_lev) %>% 
    ggplot(aes(geometry = geometry)) +
    geom_sf(data = all_sf_objects %>%
              dplyr::filter(gadm_lev == 0) %>% 
              select(-gadm_lev),
            fill = "white", size = .1) +
    geom_sf(fill = "darkgrey", color = "white", size = .03) +
    geom_sf(data = all_sf_objects %>%
              dplyr::filter(gadm_lev == 0) %>% 
              select(-gadm_lev), alpha = 0, size = .1) +
    ggthemes::theme_map() +
    facet_grid(`observation length` ~ `Administrative level`, labeller = label_both)
  
  if (save_fig) {
    ggplot2::ggsave(p_coverage, 
                    filename = paste0(getFigDir(), "/", "data_coverage.png"),
                    width = 8,
                    height = 6,
                    dpi = 300)
  }
  
  p_coverage
  
}

#' Plot all posterior coverage
#' Plot of the coverage of observations for different CrI widths
#'
#' @param countries 
#' @param identifier 
#' @param time_left 
#' @param time_right 
#' @param redo 
#' @param save_fig 
#' @param fig_dir 
#' @param verbose 
#'
#' @return return
#' 
plotAllCoverage <- function(countries = "all",
                            identifier,
                            time_left = NULL,
                            time_right = NULL,
                            redo = F,
                            save_fig = F,
                            fig_dir = getFigDir(),
                            verbose = T) { 
  
  all_sims <- runAll(
    countries = countries,
    identifiers = identifier,
    models = "best",
    run_levels = "best",
    fun = getSims,
    fun_name = "sim_stats",
    redo = redo
  )
  
  all_sims <- all_sims %>% mutate(lev = as.character(lev)) %>% 
    bind_rows(all_sims %>% mutate(lev = "all")) %>% 
    group_by(country) %>% 
    mutate(gadm_lev = case_when(run_level == 1 & lev == 1 ~ "GADM1",
                                run_level == 1 & lev == 2 ~ "GADM0",
                                run_level == 2 & lev == 1 ~ "GADM2",
                                run_level == 2 & lev == 2 ~ "GADM1",
                                run_level == 2 & lev == 3 ~ "GADM0",
                                lev == "all" ~ "all",
                                T ~ NA_character_),
           gadm_lev = factor(gadm_lev, levels = c("GADM0", "GADM1", "GADM2", "all")))
  
  widths <- c(seq(.2, .9, by = .1), .95)
  
  all_coverage <- map_df(
    widths,
    function(x) {
      bounds <- widthToInt(x)
      qleft <- str_c("q", str_trim(formatC(bounds[1]*100, format = "g", digits  = 2)))
      qright <- str_c("q", str_trim(formatC(bounds[2]*100, format = "g", 
                                            digits  = ifelse(bounds[2] == 0.975, 3, 2))))
      
      all_sims %>% 
        dplyr::mutate(obs_low = !!sym(qleft),
                      obs_high = !!sym(qright),
                      width = x) %>% 
        dplyr::select(obs, gadm_lev, obs_low, obs_high, width)
    }) %>% 
    group_by(country, gadm_lev, width) %>% 
    summarise(coverage = sum(obs >= obs_low & obs <= obs_high)/n()) 
  
  cov_labels <- all_coverage %>% 
    dplyr::filter(width == max(width),
                  gadm_lev == "all") %>% 
    mutate(label = formatC(coverage, digits = 2) %>% str_trim())
  
  p_coverage <- all_coverage %>% 
    ggplot(aes(x = width, y = coverage, color = gadm_lev)) +
    geom_line(aes(lty = gadm_lev)) +
    geom_point(data = cov_labels) +
    ggrepel::geom_label_repel(
      data = cov_labels, aes(label = label), 
      size = 3, box.padding = 0.5, segment.linetype = 4,
      max.overlaps = Inf) +
    facet_wrap(~country) +
    ggthemes::scale_color_colorblind() +
    scale_linetype_manual(values = c(2, 2, 2, 1)) +
    guides(color = guide_legend("GADM level"),
           lty = guide_legend("GADM level")) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(.05, .25, .5, .75, .95)) +
    theme_bw() +
    labs(x = "Credible interval width", y = "Observation coverage fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom")
  
  
  
  if (save_fig) {
    ggplot2::ggsave(p_coverage, 
                    filename =  makeFigFilename(prefix = "ppc_coverage",
                                                countries = countries,
                                                run_level = "best",
                                                identifier = identifier),
                    width = 8,
                    height = 9,
                    dpi = 300)
  }
  
  p_coverage
}

plotCorrStats <- function(countries = "all",
                          identifier,
                          run_level,
                          covariates = c("Precipitation"),
                          redo_single = redo,
                          method = "pearson",
                          save_fig = F,
                          verbose = T) {
  
  
  # Correlations file
  cor_file <- makeStdResName(country = countries,
                             run_level = run_level,
                             identifier = identifier,
                             model = "best",
                             gadm_lev = "auto",
                             suffix = paste0("corr_", 
                                             paste0(covariates, collapse = "-"),
                                             "_", method),
                             file_type = "rds")
  
  
  if (file.exists(cor_file)) {
    clim_corr <- readRDS(cor_file) %>% 
      tibble::as_tibble() %>% 
      dplyr::left_join(getGadmData())
    
  } else {
    stop("Compute correlations first")
  }
  
  cors <- seq(.5, .9, by = .05)
  
  cor_stats <-  clim_corr %>% 
    mutate(popsum = popsum_2020) %>% 
    select(gid, area, lag, covar, popsum, corr) %>%
    mutate(ind = 1,
           covar = getCovarDict()[covar]) %>% 
    pivot_longer(cols = c("area", "popsum", "ind"),
                 names_to = "var") %>% 
    group_by(covar, lag, var) %>% 
    group_modify(function(x, y) {
      map_df(cors, function(z) {
        x %>% 
          summarise(frac = sum(value[abs(corr)>=z], na.rm=T)/sum(value[!is.na(corr)])) %>% 
          mutate(corr_thresh = z)
      })
    })
  
  var_dict <- c("area" = "Area",
                "ind" = "Administrative\nunits",
                "popsum" = "Population")
  
  p_corr_stats <- cor_stats %>% 
    ggplot(aes(x = corr_thresh, y = frac)) +
    geom_line(aes(color = covar, lty = covar)) +
    facet_grid(var ~ lag, labeller = labeller(lag = label_both,
                                              var = var_dict)) +
    theme_bw() +
    labs(x = "Absolute value of Spearman correlation", 
         y = "Fraction with correlation above threshold") +
    scale_color_manual(values = c("blue", "red", "green"))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_corr_stats, 
                    filename =  makeFigFilename(prefix = "corr_stats",
                                                countries = countries,
                                                run_level = run_level,
                                                identifier = identifier),
                    width = 7,
                    height = 4,
                    dpi = 300)
  }
  
}

#' Get prior or posterior selection
#' Get the prior or posterior statistics for a set of parameters
#' 
#' @param country 
#' @param run_level 
#' @param model 
#' @param identifier 
#' @param time_left 
#' @param time_right 
#' @param gadm_lev 
#' @param redo 
#' @param verbose 
#' @param what whether to select the prior or the posterior
#' @param ... 
#'
#' @return
#' 
getParamStats <- function(country,
                          run_level,
                          model = "best",
                          identifier,
                          time_left = NULL,
                          time_right = NULL,
                          gadm_lev = setGadmLev(country, run_level),
                          redo = getRedoDefaults(),
                          verbose = F,
                          what = "prior",
                          ...) {
  
  redo <- completeRedo(redo)
  
  if (what == "prior") {
    suffix <- "_prior"
  } else if (what == "posterior") {
    suffix <- ""  
  } else {
    stop("Set ", what, " unkown. Needs to be one of prior/posterior")
  }
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo$data,
                                 verbose = verbose)  
  }
  
  if (verbose) {
    cat("Getting betas for", country, run_level, identifier, "\n")
  }
  
  if (model == "best") {
    model <- getBestModel(country = country,
                          run_level = run_level,
                          identifier = identifier,
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = gadm_lev,
                          redo = redo$best_model)
  }  
  
  # Get result file
  res_file <- try(getResFiles(country = country,
                              run_level = run_level, 
                              identifier = identifier, 
                              time_left = time_left,
                              time_right = time_right,
                              gadm_lev = gadm_lev,
                              what = str_glue("output{suffix}"),
                              paths = getResPaths(what = str_glue("outputs{suffix}")),
                              model = model))
  
  output <- readRDS(res_file)
  
  param_stats <- getOutputStats(output) %>% 
    bind_cols(parseResFilename(res_file)) %>% 
    mutate(set = what)
  
  return(param_stats)
}

#' Get output statistics
#' Gets the stats for selected model parameters from the output objects
#' @param output 
#'
#' @return a tibble
#'
getOutputStats <- function(output) {
  bind_rows(
    output$par_mod_df %>%
      janitor::clean_names() %>%
      magrittr::set_colnames(colnames(.) %>% 
                               str_replace("x2_", "x02_") %>%
                               str_remove_all("_percent|_") %>%
                               str_replace_all("x", "q")
      ),
    output$par_conv_df %>%
      janitor::clean_names() %>%
      group_by(param, param_fam) %>%
      summarise(mean = mean(sample),
                sd = sd(sample),
                q025 = quantile(sample, .025),
                q25 = quantile(sample, .25),
                q50 = quantile(sample, .50),
                q75 = quantile(sample, .75),
                q975 = quantile(sample, .975)) %>%
      ungroup()
  )}

#' Plot prior shrinkage
#' Plot the shrinkage betwee the prior and posterior of parameters. This figure only shows
#' seasonality-related parameters
#'
#' @param countries 
#' @param identifier 
#' @param redo 
#' @param redo_single 
#' @param save_fig 
#' @param verbose 
#'
#' @return
#'
plotAllPriorShrinkage <- function(countries = "all",
                                  identifier,
                                  redo = F,
                                  redo_single = redo,
                                  save_fig = F,
                                  verbose = T) {
  
  best_priors <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = "best",
    models = "best",
    fun = getParamStats,
    fun_name = "best_priors",
    fun_opts = list(redo = list(main = redo_single),
                    what = "prior"),
    redo = redo
  )
  
  best_posteriors <- runAll(
    countries = countries,
    identifiers = identifier,
    run_levels = "best",
    models = "best",
    fun = getParamStats,
    fun_name = "best_posteriors",
    fun_opts = list(redo = list(main = redo_single),
                    what = "posterior"),
    redo = redo
  )
  
  # Combine sets
  param_stats <- dplyr::bind_rows(best_priors, 
                                  best_posteriors)
  
  # Plot shrinkage for seasonality-related parameters
  p_shrinkage <- param_stats %>%
    dplyr::select(param, param_fam, country, sd, set) %>%
    tidyr::pivot_wider(values_from = "sd",
                       names_from = "set") %>%
    dplyr::mutate(param_fam = dplyr::case_when(
      stringr::str_detect(param, "\\[") ~ stringr::str_extract(param, ".*(?=\\[)"),
      T ~ as.character(param))) %>% 
    dplyr::filter(str_detect(param_fam, "beta|tau_betas"),
                  param_fam != "tau_lambdas",
                  !(country %in% getNoSeasCountries()),
                  !(country %in% getExcludedCountries())) %>% 
    dplyr::mutate(param_fam = factor(param_fam,
                                     levels = c("tau_betas", "betas_g1", "betas_g2"),
                                     labels = c("Precision", "Seasonality first group", "Seasonality second group"))) %>%
    ggplot() +
    geom_histogram(aes(x = 1-posterior^2/prior^2, fill = param_fam)) +
    geom_vline(aes(xintercept = 0), col = "gray20", size = .5, lty = 2) +
    theme_bw() +
    labs(x = "Parameter shrinkage") +
    guides(fill = guide_legend("Parameters")) +
    facet_wrap(~country, scales = "free_y") +
    scale_fill_manual(values = c("gray", getGroupingColors()))
  
  if (save_fig) {
    ggsave(p_shrinkage,   
           filename = paste0(getFigDir(), "/prior_shrinkage",
                             "_rlbest_",
                             identifier,
                             "_",
                             makePeriod(NULL, NULL),
                             ".png"),
           width = 10,
           height = 6,
           dpi = 300)
  }
}

plotClusteringModelComp <- function(identifier,
                                    suffix = "sigma1_tau0.5",
                                    save_fig = F) {
  
  # Get clustering outputs and compute groups
  fit_files <- dir("generated_data/output_clustering/", 
                   pattern = "kmeans_output", 
                   full.names = T) %>% 
    stringr::str_subset(pattern = suffix)
  
  group_loo_comp <- map(fit_files, function(x) {readRDS(x) %>% .$loo}) %>%
    loo::loo_compare()
  
  group_loo_comp <- group_loo_comp %>% 
    as_tibble() %>% 
    mutate(model = rownames(group_loo_comp) %>%
             str_extract("[0-9]+") %>% 
             as.numeric()) %>% 
    mutate(k = map_dbl(model, ~ str_extract(fit_files[.], "[0-9]+") %>% as.numeric))
  
  # Plot
  p_loo <- ggplot(group_loo_comp, aes(x = k, y = elpd_diff)) +
    geom_hline(aes(yintercept = 0), lty = 2, size = .3) +
    geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = elpd_diff - 2 * se_diff, 
                      ymax = elpd_diff + 2 * se_diff),
                  size = .5,
                  width = 0) +
    theme_bw() +
    labs(x = "Number of groups", y = "Difference in LOO-CV estimate\nwith respect to maximum")
  
  if (save_fig) {
    
    ggsave(p_loo, 
           filename = paste0(getFigDir(), "/seas_group_model_comp",
                             "_",
                             identifier,
                             ".png"),
           width = 6,
           height = 4)
  }
  
  return(p_loo)
}

# Publication figures -----------------------------------------------------



#' @title Plot seasonality combined
#'
#' @description Plots the heatmap along with the monthly maps
#'
#' @param param
#'
#' @return ggplot object
#' 
plotAllSeasCombined2 <- function(p_seas_index = NULL,
                                 p_maps = NULL,
                                 countries = "all",
                                 run_level = "best",
                                 identifier = "mean_annual_incidence",
                                 time_left = NULL,
                                 time_right = NULL,
                                 redo = F,
                                 redo_single = F,
                                 save_fig = F,
                                 verbose = F) {
  
  if (is.null(p_seas_index)) {
    # Plot seasonality coefficients
    p_seas_index <- plotAllSeasInd2(countries = countries,
                                    identifier = identifier,
                                    run_level = run_level,
                                    time_left = time_left,
                                    time_right = time_right,
                                    redo = redo,
                                    save_fig = F)
  }
  
  if (is.null(p_maps)) {
    # Plot of seasonality coefficients as a map
    p_maps <- plotAllSeasMaps(countries = countries,
                              identifier = identifier,
                              run_level = run_level,
                              time_left = time_left,
                              time_right = time_right,
                              redo = redo,
                              redo_single = redo_single,
                              redo_data = F,
                              save_fig = F,
                              verbose = verbose)
    
  }
  
  
  # Combine plots
  p_seas_combined <- cowplot::plot_grid(
    p_maps +
      theme(legend.position = "bottom",
            plot.margin = unit(c(1, 1, 1, 1), "lines"),
            legend.key.height = unit(.2, unit = "cm")) +
      guides(fill = guide_colorbar("Odds ratio",
                                   title.position = "top")),
    p_seas_index +
      theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
    labels = "auto",
    rel_widths = c(.9, 1)
  ) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  if (save_fig) {
    ggsave(p_seas_combined,   
           filename = paste0(getFigDir(), "/seas_coef_combine_all_v2",
                             "_rlbest_",
                             identifier,
                             "_",
                             makePeriod(NULL, NULL),
                             ".png"),
           width = 9,
           height = 6,
           dpi = 300)
  }
  
  p_seas_combined
}



# A function to plot the inset map
plotMap <- function(df, 
                    identifier) {
  if (str_detect(df$variant[1], "mixture")) {
    plotGroupingWeights(country = df$country[1],
                        variant = df$variant[1],
                        run_level = df$run_level[1],
                        identifier = identifier,
                        sf_object = dplyr::filter(sf_objects, country == df$country[1])) +
      guides(fill = "none") +
      theme(plot.background = element_rect(fill = "white", colour = "white"),
            plot.margin = unit(c(0, 0, 0, 0), unit = "lines"))
  } else {
    df <- data.frame()
    ggplot(df) + geom_point()+
      theme(plot.background = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), unit = "lines"))
  }
}

#' @title Plot all seasonality time series
#'
#' @description Plots the seasonal coefficients``
#'
#' @param param
#'
#' @return 
#' 
plotAllSeasTS <- function(countries = "all",
                          identifier,
                          run_level,
                          time_left = NULL,
                          time_right = NULL,
                          best_seas_grp = NULL,
                          sf_objects = NULL,
                          redo_single = F,
                          redo_data = F,
                          redo = F,
                          save_fig = F,
                          fig_dir = getFigDir(),
                          verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas_grp)) {
    best_seas_grp <- getBestSeasByGroup(countries = countries,
                                        identifier = identifier,
                                        run_level = run_level,
                                        redo = list(main = redo))
    
  } else {
    all_countries <- unique(best_seas_grp$country)
  }
  
  
  if (is.null(sf_objects)) {
    sf_objects <- getSfObjects(best_seas_grp)
  }
  
  # Tidy with month number
  best_seas_grp <- tidyBestSeasGrp(best_seas_grp, all_countries)
  
  # Make inset maps of the prob of being in each group
  # Based on solution in:
  # https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
  
  insets <- best_seas_grp %>% 
    dplyr::filter(!is.na(mean)) %>% 
    split(f = .$country) %>%
    purrr::map(~annotation_custom2(
      grob = ggplotGrob(plotMap(., identifier = identifier)),
      data = data.frame(country = unique(.$country)),
      xmin = 0, xmax = 4,
      ymin = 2.5, ymax = ifelse(unique(.$country) != "SDN", 5.7, 9)
    ))
  
  # Add fake data to set axis limits
  fake_data <- best_seas_grp %>% 
    filter(!is.na(mean)) %>% 
    group_by(country) %>% 
    summarise(x = 0,
              ymax = ifelse(country[1] == "SDN", exp(9), exp(5.7)),
              ymin = ifelse(country[1] == "SDN", exp(-12.5), exp(-5)))
  
  p_best_seas_ts <- best_seas_grp %>%
    dplyr::filter(!is.na(mean)) %>%
    group_by(country) %>% 
    group_modify(function(x, y) {
      og <- getOrderGroups(x)
      x %>% 
        mutate(group = ifelse(param_fam == "betas_g1", og[1], og[2]))
    }) %>% 
    ggplot() +
    geom_hline(aes(yintercept = exp(0)), lty = 5, size = .4, color = "darkgray") +
    geom_ribbon(aes(x = month, ymin = exp(q025), ymax = exp(q975), fill = group), alpha = .2) +
    geom_line(aes(x = month, y = exp(mean), color = group)) +
    geom_point(data = fake_data, 
               inherit.aes = F,
               aes(x = x, y = ymax), 
               size = 0, col = "white", alpha = 0) +
    geom_point(data = fake_data, 
               inherit.aes = F,
               aes(x = x, y = ymin), 
               size = 0, col = "white", alpha = 0) +
    facet_wrap(~country, ncol = 4, scales = "free_y") +
    labs(x = "Month", y = "Odds ratio") +
    scale_color_manual(values = getGroupingColors()) +
    scale_fill_manual(values = getGroupingColors()) +
    scale_x_continuous(breaks = seq(1, 12, by = 2), labels = month.abb[seq(1, 12, by = 2)]) +
    scale_y_continuous(trans = "log", breaks = 10^(-3:3), labels = formatC(10^(-3:3))) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank()) +
    guides(lty = "none", color = "none", fill = "none")
  
  
  if (save_fig) {
    ggsave(p_best_seas_ts + insets,   
           filename = paste0(getFigDir(), "/seas_coef_ts_all",
                             "_rlbest_",
                             identifier,
                             "_",
                             makePeriod(NULL, NULL),
                             ".png"),
           width = 9,
           height = 12,
           dpi = 600)
  }
  
  p_best_seas_ts + insets
}



#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasGroups <- function(countries = "all",
                              identifier,
                              run_level,
                              time_left = NULL,
                              time_right = NULL,
                              best_seas = NULL,
                              sf_objects = NULL,
                              africa.sf_lev0 = NULL,
                              n_groups = "best",
                              redo = F,
                              save_fig = F,
                              fig_dir = getFigDir(),
                              verbose = T,
                              label = "auto") {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas)) {
    best_seas <- getBestSeas(countries = countries,
                             identifier = identifier,
                             run_level = run_level,
                             redo = list(main = redo))
  }
  
  # if (is.null(sf_objects)) {
  #   sf_objects <- getSfObjects(best_seas)
  # }
  
  if (is.null(africa.sf_lev0)) {
    africa.sf_lev0 <- getCountrySf(all_countries)
  }
  
  # Get clustering outputs
  fit_files <- dir("generated_data/output_clustering/", 
                   pattern = "kmeans_output", 
                   full.names = T) %>% 
    stringr::str_subset(pattern = "sigma1_tau0.5")
  
  cluster <- fit_files %>% 
    stringr::str_subset(glue::glue("_{n_groups}_")) %>% 
    readRDS()
  
  cluster$sf_object <- cluster$sf_object %>% 
    rmapshaper::ms_simplify(keep = .2)
  
  best_seas_data <- inner_join(best_seas ,
                               cluster$sf_object %>% 
                                 as_tibble() %>% 
                                 dplyr::select(gid, cluster)) 
  
  # Get gids for countries for which seasonality was not selected
  sf_objects <- dplyr::bind_rows(
    purrr::map_df(getExcludedCountries(), ~getCountrySf(.)) %>% 
      dplyr::mutate(cluster = "excluded"),
    purrr::map_df(getNoSeasCountries(), ~getCountrySf(.)) %>% 
      dplyr::mutate(cluster = "non-seasonal")
  )
  
  p_seas_groups <- best_seas_data %>% 
    mutate(cluster = str_c("Cluster ", cluster)) %>% 
    ggplot(aes(x = month, y = exp(mean), color = cluster)) +
    geom_hline(aes(yintercept = 1), lty = 2, size = .1) +
    geom_line(aes(group = gid), alpha = .01) +
    geom_smooth(size = 1, se = F) +
    facet_wrap(~cluster, ncol = ifelse(cluster$k == 2 , 1, 2)) +
    theme_bw() +
    scale_color_manual(values = getNClusterColors(k = as.numeric(n_groups))) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_y_continuous(trans = "log", breaks = c(.1, 1, 10), labels = formatC(c(.1, 1, 10))) +
    labs(y = "Odds ratio", x = "Month") +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          panel.grid.minor = element_blank()) +
    guides(color = "none") +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  
  # Lake data from 
  lakes_sf <- getLakes()
  
  p_groups_map <- cluster$sf_object %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(!(country %in% cluster$sf_object$country)), 
            fill = "lightgray", color = "darkgray", size = .3) +
    geom_sf(aes(fill = factor(cluster), alpha = max_prob_cat), color = "white", size = .02) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "white", size = .02, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .25, alpha = 0) +
    scale_alpha_discrete(range = c(.7, 1)) +
    scale_fill_manual(values = getNClusterColors(k = as.numeric(n_groups))) +
    ggthemes::theme_map()  +
    guides(fill = "none",
           alpha = "none") +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  
  p_grouping <- cowplot::plot_grid(
    p_groups_map, 
    p_seas_groups +
      theme(plot.margin = unit(c(2, 1, 2, 1), units = "lines")),
    rel_widths = c(1.4, ifelse(cluster$k == 2, .7, 1.3)),
    labels = label) +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_grouping, 
                    filename = paste0(getFigDir(), "/seas_group_map_",
                                      stringr::str_c(countries, collapse = "-"),
                                      "_k",
                                      n_groups,
                                      "_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 9,
                    height = 6,
                    dpi = 500)
  }
  p_grouping
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllSeasGroups2 <- function(countries = "all",
                               identifier,
                               run_level,
                               time_left = NULL,
                               time_right = NULL,
                               best_seas = NULL,
                               sf_objects = NULL,
                               africa.sf_lev0 = NULL,
                               redo_single = F,
                               redo_data = F,
                               n_groups = "best",
                               redo = F,
                               save_fig = F,
                               fig_dir = getFigDir(),
                               verbose = T,
                               label = "auto") {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(best_seas)) {
    best_seas <- getBestSeas(countries = countries,
                             identifier = identifier,
                             run_level = run_level,
                             redo = list(main = redo))
  }
  
  if (is.null(sf_objects)) {
    sf_objects <- getSfObjects(best_seas)
  }
  
  if (is.null(africa.sf_lev0)) {
    africa.sf_lev0 <- getCountrySf(all_countries)
  }
  
  # Get clustering outputs and comput groups
  fit_files <- dir("generated_data/output_clustering/", 
                   pattern = "kmeans_output", 
                   full.names = T) %>% 
    stringr::str_subset(pattern = "sigma1_tau0.5")
  
  cluster_dict_2g <- c("Cluster 1" = "SDN",
                       "Cluster 2" = "Late peak - high amplitude",
                       "Cluster 3" = "Early peak - high amplitude")
  
  cluster_2g <- fit_files %>% 
    stringr::str_subset("_3") %>% 
    readRDS() 
  
  cluster_2g$sf_object <- cluster_2g$sf_object %>% 
    mutate(cluster = str_c("Cluster ", cluster),
           cluster = cluster_dict_2g[cluster],
           cluster = factor(cluster, levels = cluster_dict_2g))
  
  cluster_dict_5g <- c("Cluster 1" = "SDN",
                       "Cluster 3" = "Late peak - high amplitude",
                       "Cluster 5" = "Late peak - low amplitude",
                       "Cluster 4" = "Early peak - high amplitude",
                       "Cluster 2" = "Early peak - low amplitude",
                       "non-seasonal" = "non-seasonal",
                       "excluded" = "excluded")
  
  cluster_dict_5g_order <- c(1, 2, 3, 4, 5, 6, 7)
  
  cluster_5g <- fit_files %>% 
    stringr::str_subset("_5") %>% 
    readRDS() 
  
  cluster_5g$sf_object <- cluster_5g$sf_object %>% 
    # rmapshaper::ms_simplify(keep = .1) %>%
    dplyr::mutate(cluster = str_c("Cluster ", cluster),
                  cluster = cluster_dict_5g[cluster],
                  cluster = factor(cluster, levels = cluster_dict_5g[cluster_dict_5g_order]))
  
  best_seas_data_2g <- dplyr::inner_join(best_seas ,
                                         cluster_2g$sf_object %>% 
                                           dplyr::as_tibble() %>% 
                                           dplyr::select(gid, cluster))
  esri_crs <- st_crs()
  esri_crs$wkt <- 'PROJCS["Africa_Sinusoidal",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",15],UNIT["Meter",1]]'    
  
  # Make union of 2g to make outline
  union_2g <- dplyr::inner_join(cluster_2g$sf_object %>% 
                                  dplyr::mutate(
                                    cluster = as.character(cluster),
                                    cluster = dplyr::case_when(country == "SDN" ~ cluster_dict_2g[1],
                                                               country == "SLE" ~ cluster_dict_2g[2],
                                                               country == "CMR" ~ cluster_dict_2g[2],
                                                               country == "GIN" ~ cluster_dict_2g[2],
                                                               country == "NGA" ~ cluster_dict_2g[2],
                                                               country == "SSD" ~ cluster_dict_2g[3],
                                                               T ~ cluster),
                                    cluster = factor(cluster, levels = cluster_dict_2g)
                                  ) %>% 
                                  sf::st_transform(crs = esri_crs) %>%
                                  sf::st_make_valid() %>% 
                                  lwgeom::st_snap_to_grid(2000) %>% 
                                  sf::st_buffer(0) %>% 
                                  dplyr::select(gid, cluster),
                                best_seas) %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(n = n())
  
  # Remove wholes
  union_2g <- nngeo::st_remove_holes(union_2g) %>% 
    sf::st_transform(crs = 4326) 
  
  # To parts to remove small regions
  union_2g_parts <- sf::st_cast(union_2g, "POLYGON")
  
  # Remove small and fill wholes
  union_2g_filled <- union_2g_parts %>% 
    dplyr::filter(as.numeric(st_area(geometry)) > 1e11) %>% 
    nngeo::st_remove_holes() %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(n = n()) %>% 
    rmapshaper::ms_simplify(keep = .2)
  
  best_seas_data_5g <- dplyr::inner_join(best_seas ,
                                         cluster_5g$sf_object %>% 
                                           dplyr::as_tibble() %>% 
                                           dplyr::select(gid, cluster)) 
  
  # Get gids for countries for which seasonality was not selected
  sf_objects <- dplyr::bind_rows(
    purrr::map_df(getExcludedCountries(), ~getCountrySf(.)) %>% 
      dplyr::mutate(cluster = "excluded"),
    purrr::map_df(getNoSeasCountries(), ~getCountrySf(.)) %>% 
      dplyr::mutate(cluster = "non-seasonal")
  )
  
  # Compute statistics on the number of people
  cluster_data <-  best_seas_data_5g %>% 
    dplyr::filter(month == 1) %>% 
    dplyr::inner_join(getGadmData()) %>% 
    dplyr::bind_rows(
      dplyr::inner_join(sf_objects %>% 
                          dplyr::as_tibble() %>% 
                          dplyr::select(country, gid, cluster),
                        getGadmData())
    ) %>% 
    dplyr::mutate(cholera_rates_2010_2016 = ifelse(is.na(cholera_rates_2010_2016), 0, cholera_rates_2010_2016),
                  rate_cat = cut(cholera_rates_2010_2016, getRateCuts()),
                  rate_cat = getRateDict()[rate_cat],
                  rate_cat = factor(rate_cat, levels = getRateDict()),
                  cluster = ifelse(cluster == "SDN", cluster_dict_5g[2], cluster),
                  cluster = factor(cluster, levels = cluster_dict_5g[-1])) 
  
  pop_frac <- cluster_data %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(pop = sum(popsum_2020)) %>% 
    dplyr::mutate(frac = pop/sum(pop),
                  var = "Total", 
                  rate_cat = "Total") %>% 
    dplyr::ungroup()
  
  
  rate_frac <- cluster_data %>% 
    dplyr::group_by(rate_cat, cluster) %>% 
    dplyr::summarise(pop = sum(popsum_2020)) %>% 
    dplyr::group_by(rate_cat) %>% 
    dplyr::mutate(frac = pop/sum(pop),
                  var = "Incidence rate") %>% 
    dplyr::ungroup()
  
  
  # Save output
  saveRDS(rbind(rate_frac, pop_frac),
          file = makeStdResName(country = countries, 
                                run_level = "best", 
                                identifier = identifier, 
                                time_left = time_left,
                                time_right = time_right, 
                                gadm_lev = "best",
                                model = "best", 
                                suffix = "rate_fraction_clusters", 
                                file_type = "rds"))
  
  p_frac <- rbind(rate_frac, 
                  pop_frac) %>%
    ggplot(aes(x = rate_cat, y = frac, fill = cluster)) +
    geom_bar(stat = "identity") +
    facet_grid(~ var, scales = "free", space = "free") +
    theme_bw() +
    labs(x = "Cases/100'000/year", 
         y = "Fraction of\npopulation") +
    scale_fill_manual(values = c(getClusterColors(k = 4), "gray75", "gray90")) +
    guides(fill = guide_legend(title = "Seasonality cluster"))
  
  p_seas_groups <- best_seas_data_5g  %>% 
    mutate(cluster = ifelse(cluster == "SDN", cluster_dict_5g[2], as.character(cluster)),
           cluster = factor(cluster, levels = cluster_dict_5g)) %>%
    filter(country != "SDN") %>% 
    ggplot(aes(x = month, y = exp(mean), color = cluster)) +
    geom_hline(aes(yintercept = 1), lty = 2, size = .1) +
    geom_line(aes(group = gid), alpha = .01) +
    geom_smooth(size = 1, se = F) +
    facet_wrap(~cluster, ncol = 1, dir = "v") +
    theme_bw() +
    scale_color_manual(values = getClusterColors(k = 5)[-1]) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    coord_cartesian(ylim = c(exp(-3.5), exp(3.5))) +
    scale_y_continuous(trans = "log", 
                       breaks = c(.1, 1, 10), 
                       labels = formatC(c(.1, 1, 10))) +
    labs(y = "Odds ratio", x = "Month") +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          panel.grid.minor = element_blank()) +
    guides(color = "none")
  
  p_seas_SDN <- best_seas_data_5g  %>% 
    filter(country == "SDN") %>% 
    ggplot(aes(x = month, y = exp(mean), color = cluster)) +
    geom_hline(aes(yintercept = 1), lty = 2, size = .1) +
    geom_line(aes(group = gid), alpha = .01) +
    geom_smooth(size = .3, se = F) +
    facet_wrap(~cluster, ncol = 1, dir = "v") +
    theme_bw() +
    scale_color_manual(values = getClusterColors(k = 5)) +
    scale_x_continuous(breaks = c(1, 7, 12), labels = month.abb[c(1, 7, 12)]) +
    scale_y_continuous(trans = "log", 
                       breaks = c(.1, 10), 
                       labels = str_pad(formatC(c(.1,10)), 3)) +
    labs(y = "Odds ratio", x = "Month") +
    theme(axis.text = element_text(size = 4),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(vjust = 4),
          axis.text.y = element_text(hjust = 4, margin = margin(0, 0, 0 ,0)),
          strip.text = element_text(size = 5, margin = margin(0.3, 0.3, 0.3, 0.3)),
          plot.margin = margin(0, 0, 0, 0)) +
    guides(color = "none")
  
  # Lake data from 
  lakes_sf <- getLakes()
  
  # Shrink to fit lines 
  union_2g_filled_shrunk <- sf::st_cast(union_2g_filled, "MULTIPOLYGON") %>% 
    sf::st_cast("POLYGON", do_split = T) %>% 
    st_shrink(.9975) 
  
  # union_2g_filled_shrunk <- st_buffer(union_2g_filled, 10)
  
  p_groups_map <- cluster_5g$sf_object  %>% 
    dplyr::filter(!(country %in% getNoSeasCountries())) %>% 
    ggplot() +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), 
            fill = "gray90", color = "gray40", size = 0) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), 
            fill = "gray75", color = "gray40", size = 0) +
    geom_sf(aes(fill = cluster), color = "white", size = .02) +
    geom_sf(data = lakes_sf, fill = "#aad3df", color = "#7fb4c4", size = .06) +
    geom_sf(color = "white", size = .02, alpha = 0) +
    geom_sf(data = africa.sf_lev0, color = "white", size = .25, alpha = 0) +
    # geom_sf(data = union_2g_filled,
    #         aes(color = cluster), size = 1, alpha = 0) +
    geom_sf(data = union_2g_filled_shrunk,
            aes(color = cluster), size = .8, alpha = 0) +
    scale_alpha_discrete(range = c(.7, 1)) +
    scale_fill_manual(values = getClusterColors(k = cluster_5g$k)) +
    scale_color_manual(values = getClusterColors(k = cluster_2g$k)) +
    ggthemes::theme_map()  +
    guides(fill = "none",
           color = "none",
           alpha = "none") +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  p_groups_map2 <- cowplot::plot_grid(
    p_groups_map +
      theme(plot.margin = margin(0, 0, 0, 0, unit = "lines"),
            legend.key.height = unit(.18, unit = "cm")),
    p_frac +
      theme(plot.margin = margin(.2, 1, 1, 1, unit = "lines"),
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 9),
            # plot.background = element_blank(),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 9),
            legend.key.size = unit(.3, units = "cm"),
            strip.text = element_text(size = 8, 
                                      margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"))),
    ncol = 1,
    rel_heights = c(1, .275),
    # rel_widths = c(1,1),
    align = "v")
  
  
  p_grouping <- cowplot::plot_grid(
    p_groups_map2, 
    p_seas_groups +
      theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines")) +
      annotation_custom2(grob=ggplotGrob(p_seas_SDN), 
                         data = data.frame(cluster="Late peak - high amplitude"),
                         ymin = 0.1, ymax = 3.5, 
                         xmin = 0.75, xmax = 5.25),
    rel_widths = c(1, .5),
    labels = "auto")  +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  
  if (save_fig) {
    ggplot2::ggsave(p_grouping, 
                    filename = paste0(getFigDir(), "/seas_group_map_",
                                      stringr::str_c(countries, collapse = "-"),
                                      "_rl",
                                      run_level,
                                      "_",
                                      identifier,
                                      "_",
                                      makePeriod(time_left, time_right),
                                      ".png"),
                    width = 7.5,
                    height = 6,
                    dpi = 500)
  }
  p_grouping
}


#' @title plot all seasonality
#'
#' @description makes a tile plot of the seasonality coefficients for the best
#' model for a given identifier at a given run level
#'
#' @param param
#'
#' @return return
plotAllCovarCorr <- function(countries = "all",
                             identifier,
                             run_level,
                             covariates = c("Precipitation"),
                             best_seas = NULL,
                             sf_objects = NULL,
                             africa.sf_lev0 = NULL,
                             lags = 0, 
                             redo = F,
                             method = "pearson",
                             save_fig = F,
                             verbose = T) {
  
  if (countries == "all") {
    all_countries <- getAllSSACountries()
  } else {
    all_countries <- countries
  }
  
  if (is.null(africa.sf_lev0)) {
    africa.sf_lev0 <- getCountrySf(all_countries)
  }
  
  
  # Correlations file
  cor_file <- makeStdResName(country = countries,
                             run_level = run_level,
                             identifier = identifier,
                             model = "best",
                             gadm_lev = "auto",
                             suffix = paste0("corr_", 
                                             paste0(covariates, collapse = "-"),
                                             "_", method),
                             file_type = "rds")
  
  
  if (file.exists(cor_file) & !redo) {
    clim_corr <- readRDS(cor_file)
  } else {
    
    if (is.null(best_seas)) {
      best_seas <- getBestSeas(countries = countries,
                               identifier = identifier,
                               run_level = run_level,
                               redo = list(main = redo))
    }
    
    if (is.null(sf_objects)) {
      sf_objects <- getSfObjects(best_seas)
    }
    
    covar_file <- makeStdResName(country = countries,
                                 run_level = run_level,
                                 identifier = identifier,
                                 model = "best",
                                 gadm_lev = "auto",
                                 suffix = str_c(str_c(covariates, collapse = "-"), "_", str_c("l", str_c(lags, collapse = "-"))),
                                 file_type = "rds")
    
    # Pull covariates
    if (file.exists(covar_file) & !redo) {
      clim <- readRDS(covar_file) 
    } else {
      # Pull covariates
      clim <- best_seas %>% 
        dplyr::group_by(country) %>% 
        dplyr::group_modify(function(x, y) {
          getClimatology(country = y$country[1], 
                         admin_level = x$run_level[1], 
                         covar = covariates) %>% 
            dplyr::group_by(month, gid, covar) %>% 
            dplyr::arrange(gid, month, covar, source) %>% 
            dplyr::slice(1) %>% 
            dplyr::select(-source)
        }) 
      
      saveRDS(clim, file = covar_file)
    }
    
    # Join with sf
    clim_data <- dplyr::inner_join(clim, best_seas)
    
    # Compute lag between precip and seas
    clim_corr <- clim_data %>%
      dplyr::group_by(gid, covar) %>% 
      dplyr::group_modify(function(x, y) {
        # Compute correlations at different lags
        map_df(lags, function(l) {
          corr <- cor(x$mean, x$value[lagMonths(x$month, l)], method = method)
          corr_pval <- cor.test(x$mean, x$value[lagMonths(x$month, l)], method = method)$p.value
          tibble::tibble(
            lag = l,
            corr = corr,
            corr_pval = corr_pval
          )
        }) 
      }) %>% 
      dplyr::inner_join(sf_objects, .) 
    
    saveRDS(clim_corr, file = cor_file)
  }
  
  # Labellers
  covar_dict <- getCovarDict() 
  
  lag_dict <- str_c("lag: ", 
                    c(
                      "0 months",
                      "1 month",
                      str_c(2:5, " months")))
  
  names(lag_dict) <- 0:5
  
  # Plot
  p_corr <- ggplot(clim_corr) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getExcludedCountries()), 
            fill = "gray90", color = "gray", size = .1, alpha = 1) +
    geom_sf(data = africa.sf_lev0 %>% 
              dplyr::filter(country %in% getNoSeasCountries()), fill = "gray75", 
            color = "gray", size = .1, alpha = 1) +
    geom_sf(data = clim_corr %>%
              dplyr::filter(corr_pval > .05),
            aes(fill = corr, alpha = corr_pval < .05),
            color = "gray", size = .01) +
    geom_sf(data = africa.sf_lev0, color = "gray40", size = .1, alpha = 0) +
    geom_sf(data = clim_corr %>%
              dplyr::filter(corr_pval < .05),
            aes(fill = corr, alpha = corr_pval < .05),
            color = "black", size = .03) +
    facet_grid(covar ~ lag, labeller = labeller(covar = covar_dict, 
                                                lag = lag_dict)) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    scale_alpha_manual(values = c(.3, 1)) +
    ggthemes::theme_map() +
    theme(legend.position = "right",
          plot.margin = unit(c(1, 1, 1, 1), "lines")) +
    guides(alpha = "none",
           fill = guide_colorbar("Correlation")) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  if (save_fig) {
    n_covar <- length(unique(clim_corr$covar))
    
    ggplot2::ggsave(p_corr, 
                    filename = makeFigFilename(
                      prefix = paste0("covar_corr_", paste(covariates, collapse = "-"),
                                      "_lags", paste(lags, collapse = "-"),
                                      "_", method),
                      countries = countries,
                      run_level = run_level,
                      identifier = identifier),
                    width = 7,
                    height = min(2*n_covar, 11),
                    dpi = 600)
    
  }
  p_corr
}
