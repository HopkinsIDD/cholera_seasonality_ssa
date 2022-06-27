# Result pulling ----------------------------------------------------------

#' Extract parameters of interest
#' 
#' @description extract the parameter statistics from the outputs object
#'  
#' @param model_fit 
#' @param model_type 
#' @param parset 
#' @param trace 
#'
#' @details 
#' 
#' @return
#'
extractPars <- function(model_fit,
                        model_type = "null",
                        parset = "modeling",
                        trace = T) {
  
  if (parset == "convergence") {
    # Parameters for modeling and inference
    model_pars <- c("sigma", "sigma_etas", "rho", "xis")
    
    if (!str_detect(model_type, "null")) {
      model_pars <- c(model_pars, "tau_betas")
    }
    if (str_detect(model_type, "mixture")) {
      model_pars <- c(model_pars, "rho_lambdas", "tau_lambdas")
    }
  } else if (parset == "modeling") {
    # Parameter to check convergence
    #model_pars <- c("etas", "xis", "theta", "phi", "convolved_re")
    model_pars <- c("etas", "xis", "convolved_re")
    
    if (!str_detect(model_type, "null")) {
      model_pars <- c(model_pars, "betas_g1")
    }
    
    if (str_detect(model_type, "mixture")) {
      model_pars <- c(model_pars, "lambdas", "betas_g2")
    }
  } else {
    stop("Parset not known")
  }
  
  if (trace) {
    par_df <- map_df(model_pars, 
                     function(x) {
                       out <- try(rstan::extract(model_fit, pars = x, permute = FALSE))
                       if (!inherits(out, "try-error")) {
                         out  %>% 
                           as_tibble() %>% 
                           mutate(sim = row_number()) %>% 
                           gather(chain, sample, -sim) %>% 
                           mutate(param = map_chr(str_split(chain, "\\."), ~.[2]),
                                  chain = str_extract(chain, "(?<=chain:)[0-9]{1}"))
                       }
                     })
  } else {
    par_df <- rstan::summary(model_fit, pars = model_pars)$summary %>% 
      {mutate(as_tibble(.), param = row.names(.))}
  }
  
  par_df <- par_df %>% 
    mutate(param_fam = map_chr(str_split(param, "\\["), ~.[1]),
           param = factor(param, levels = unique(param)))
  
  return(par_df)
}

#' Extract offset
#'
#' @param model_type 
#' @param chol_stanfit 
#'
#' @return
#'
extractOffset <- function(model_type, chol_stanfit) {
  has_offset <- str_detect(model_type, "offset")
  if (has_offset) {
    offset_probs <- rstan::extract(chol_stanfit, 
                                   pars = ifelse(str_detect(model_type, "mixture"),
                                                 "ll_z2",
                                                 "ll_z1"))[[1]]
    
    offset_probs <- colMeans(offset_probs[!is.nan(offset_probs[,1]), ])
    if (sum(offset_probs) < 1) {
      offset_probs <- offset_probs * 12
    }
  } else {
    offset_probs <- NULL
  }
  return(offset_probs)
}

#' @title Extract to tibble
#'
#' @description Converts a stan sample extraction without permutation to a 
#' tibble with iteration and chain columns
#'
#' @param mat the 2 or 3 dimensional array of the extraction 
#'
#' @return a tibble
#' 
extractToTibble <- function(mat) {
  nchain <- dim(mat)[2]
  pars <- dimnames(mat)$parameters
  npars <- length(pars)
  
  samples <- purrr::map_df(1:nchain, function(i) {
    {
      if (length(pars) > 1) {
        mat[, i, ] %>% 
          tibble::as_tibble() 
      } else {
        mat[, i, 1] %>% 
          tibble::as_tibble() %>% 
          magrittr::set_colnames(pars)
      }
    } %>% 
      tidyr::pivot_longer(cols = dplyr::everything(),
                          values_to = "sample",
                          names_to = "param") %>% 
      dplyr::group_by(param) %>% 
      dplyr::mutate(iter = dplyr::row_number(),
                    chain = i) %>% 
      dplyr::ungroup()})
  
  return(samples)
}

#' @title Get unit Betas
#'
#' @description Get the admin-unit level betas
#'
#' @param country
#' @param run_level
#' @param model
#' @param identifier
#'
#' @return df
#' 
getBetas <- function(country,
                     run_level,
                     model = "best",
                     identifier,
                     time_left = NULL,
                     time_right = NULL,
                     gadm_lev = setGadmLev(country, run_level),
                     redo = getRedoDefaults(),
                     verbose = F,
                     ...) {
  
  redo <- completeRedo(redo)
  
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
  
  # First try to get only betas simulations
  res_file <- try(getResFiles(country = country,
                              run_level = run_level, 
                              identifier = identifier, 
                              time_left = time_left,
                              time_right = time_right,
                              gadm_lev = gadm_lev,
                              what = "output_sim",
                              paths = getResPaths(what = "outputs_sim"),
                              model = model,
                              keep = "betas"))
  
  if (inherits(res_file, "try-error")) {
    # The try object with all simulations
    res_file <- try(getResFiles(country = country,
                                run_level = run_level, 
                                identifier = identifier, 
                                time_left = time_left,
                                time_right = time_right,
                                gadm_lev = gadm_lev,
                                what = "output_sim",
                                paths = getResPaths(what = "outputs_sim"),
                                model = model,
                                ignore = "ts_probs"))
    
  }
  
  betas <- readRDS(res_file)$beta_stats
  
  if (is.null(betas)) {
    print("hello")
    betas <- tibble(variable = NA, mean = NA, q2.5 = NA, q25 = NA, q75 = NA, q97.5 = NA, unit = NA, month = NA)
  }
  
  betas <- betas %>% 
    bind_cols(parseResFilename(res_file)) %>% 
    rename(gid = unit)
  
  return(betas)
}


#' @description Get result files from a set of paths
#' 
#' @param country what country to get
#' @param run_level what run levels to get
#' @param identifier
#' @param paths what paths to look into
#' @param what one of output or fit
#' 
#' @return a vector of absolute file paths
getResFiles <- function(country, 
                        run_level,
                        identifier,
                        model = "all",
                        time_left = NULL,
                        time_right = NULL,
                        gadm_lev = setGadmLev(country, run_level),
                        paths = getResPaths(),
                        what = "output",
                        verbose = F,
                        ignore = NULL,
                        keep = NULL) {
  
  # Make period
  period <- makePeriod(opt = NULL, time_left, time_right)
  
  res_files <- dir(paths, full.names = T, pattern = what) %>% 
    stringr::str_subset(country) %>% 
    {
      if (identifier == "occurrence") {
        stringr::str_subset(., "cases|mean", negate = T)
      } else {
        stringr::str_subset(., identifier)
      }
    } %>% 
    {
      if (!is.null(keep)) {
        stringr::str_subset(., keep)
      } else {
        .
      }
    } %>% 
    {
      if (!is.null(ignore)) {
        stringr::str_subset(., ignore, negate = T)
      } else {
        .
      }
    } %>% 
    stringr::str_subset(paste0("runlev", run_level)) %>% 
    stringr::str_subset(period) %>% 
    stringr::str_subset(paste0("l", gadm_lev)) %>% 
    stringr::str_subset("rds") 
  
  if (length(res_files) == 0)
    stop("Couldn't find any results for ", country, " ", identifier, " runlev ", 
         run_level, " period ", period, " gadm_lev ", gadm_lev, " of type ", what, " in paths ", paste(paths, collapse = ", "))
  
  # Handle the case for multiple runs by selecting the most recent result
  res_files_df <- purrr::map_df(res_files, function(x) {
    parseResFilename(x) %>% 
      mutate(dir = stringr::str_split(x, "/")[[1]][2],
             f = x,
             time = file.info(x)$ctime) 
  }) %>% 
    dplyr::group_by(variant, identifier) %>%
    dplyr::arrange(variant, identifier, time) %>% 
    dplyr::slice_tail(n = 1)
  
  if (model != "all") {
    res_files <- res_files_df %>% 
      dplyr::filter(variant == model) %>% 
      .$f
  } else {
    res_files <- res_files_df$f
  }
  
  if (verbose)
    cat("Found", length(res_files), "result files for ", country, identifier,
        "runlev", run_level, "period", period, "\n")
  
  return(res_files)
}

#' @title Get best run level
#'
#' @description gets the best run level for final results
#'
#' @param country what country to get
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param time_left
#' @param time_right
#' @param redo whether to redo the extraction
#' @param verbose 
#' 
#' @details Takes in either the dataframe of model ics or the country, run_level
#' and path arguments to get the results
#' 
#' @return a tibble with the best model name for the given country, run_level and identifier
#' 
getBestRunLevel <- function(country,
                            identifier = "occurrence",
                            time_left = "all",
                            time_right = "all",
                            redo = FALSE,
                            verbose = F,
                            ...) {
  if (verbose) {
    cat("-- Selecting best run level for ", country, identifier, time_left, time_right, "\n")
  }
  
  run_data <- getRunDataWrapper(country = country,
                                identifier = identifier,
                                time_left = time_left,
                                time_right = time_right)
  
  run_level <- runLevel(country_iso3 = country,
                        run_data = run_data,
                        cholera_directory = getCholeraDirectory(),
                        verbose = verbose)
  # Enforced levels
  # run_level <- dplyr::case_when(country == "BEN" ~ 1,
  #                               country == "CIV" ~ 1,
  #                               T ~ run_level)
  
  if (verbose) {
    cat(" |---> Selected level ", run_level, "\n")
  }
  
  return(run_level)
}

#' @title Get parameter samples
#'
#' @description Gets the samples of a given set of parameter 
#' names. The function first checks if parameter have not already been extracted.
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param model what model to process; one of null, base, mixture, offset, mixture_offset
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param paths what paths to look into
#' @param parameter what parameter to extract
#' @param redo whether to redo the extraction if it has alread been done
#' @param write whether to write the extraction to file for further use

#' @return return
getParamSamples <- function(country,
                            run_level,
                            model = "all",
                            identifier = "occurrence",
                            time_left = NULL,
                            time_right = NULL,
                            gadm_lev = setGadmLev(country, run_level),
                            paths = getResPaths(),
                            parameter = NULL,
                            redo = F,
                            write = T,
                            verbose = F) {
  
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo_data,
                                 verbose = verbose) 
  }
  
  runChecks(country = country,
            run_level = run_level,
            identifier = identifier,
            model = model)
  
  if (is.null(parameter))
    stop("Please provide a paramter to extract")
  
  param_filename <- makeStdResName(country = country,
                                   run_level = run_level,
                                   identifier = identifier,
                                   model = model,
                                   time_left = time_left,
                                   time_right = time_right,
                                   gadm_lev = gadm_lev,
                                   suffix = paste0(makeParamFilename(parameter), "_samples"))
  
  if (file.exists(param_filename) & !redo) {
    return(readr::read_csv(param_filename, col_types = readr::cols()))
  } else {
    
    # Result files to parse
    res_files <- getResFiles(country = country,
                             run_level = run_level, 
                             identifier = identifier, 
                             time_left = time_left,
                             time_right = time_right,
                             gadm_lev = gadm_lev,
                             what = "fit",
                             paths = paths,
                             model = model)
    
    param_df <- purrr::map_df(res_files, function(x) {
      
      # Load parsed output file
      res <- readRDS(x)
      stan_pars <- try(names(rstan::get_inits(res$chol_stanfit)[[1]]), 
                       silent = T)
      
      if (inherits(stan_pars, "try-error")) {
        stan_pars <- getParNames(res$chol_stanfit)
      }
      
      purrr::map_df(parameter, function(y) {
        par <- ifelse(stringr::str_detect(y, "\\["), stringr::str_extract(y, "(.)*(?=\\[)"), y)
        
        if (any(stringr::str_detect(stan_pars, paste0("^", par)))) {
          
          # Extract parameter
          res$chol_stanfit %>%
            rstan::extract(par = y, permuted = F) %>%
            extractToTibble(.) %>% 
            dplyr::bind_cols(parseResFilename(x))
        }
      })
    })
    
    if (write) {
      readr::write_csv(param_df, file = param_filename)
    }
  }
  
  return(param_df)
}

#' @title Get parameters
#'
#' @description Gets the posterior mean and 95% CrI of a given set of parameter 
#' names. The function first checks if parameter have not already been extracted.
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param model what model to process; one of null, base, mixture, offset, mixture_offset
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param paths what paths to look into
#' @param parameter what parameter to extract
#' @param redo whether to redo the extraction if it has alread been done
#' @param write whether to write the extraction to file for further use

#' @return return
getParamPosteriors <- function(country,
                               run_level,
                               model = "all",
                               identifier = "occurrence",
                               time_left = NULL,
                               time_right = NULL,
                               gadm_lev = setGadmLev(country, run_level),
                               paths = getResPaths(),
                               parameter = NULL,
                               redo = getRedoDefaults(),
                               write = T,
                               verbose = T) {
  
  redo <- completeRedo(redo)
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo$data,
                                 verbose = verbose)  
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
  
  if (verbose) {
    cat("Getting ", parameter, "from", country, run_level, model, gadm_lev, "\n")
  }
  
  runChecks(country = country,
            run_level = run_level,
            identifier = identifier,
            model = model)
  
  if (is.null(parameter))
    stop("Please provide a paramter to extract")
  
  param_filename <- makeStdResName(country = country,
                                   run_level = run_level,
                                   identifier = identifier,
                                   model = model,
                                   time_left = time_left,
                                   time_righ = time_right,
                                   gadm_lev = gadm_lev,
                                   suffix = paste(parameter, collapse = "-"))
  
  if (file.exists(param_filename) & !redo$main) {
    
    res <- readr::read_csv(param_filename, col_types = readr::cols())
    
    if (nrow(res) > 0){
      return(res)
    }
  } else {
    
    # Result files to parse
    res_files <- getResFiles(country = country,
                             run_level = run_level, 
                             identifier = identifier, 
                             paths = paths,
                             model = model,
                             gadm_lev = gadm_lev,
                             time_left = time_left,
                             time_righ = time_right)
    
    param_df <- purrr::map_df(res_files, function(x) {
      # Load parsed output file
      res <- readRDS(x)
      
      # Extract parameter
      purrr::map_df(parameter, function(x) {
        res$par_mod_df %>%
          tibble::as_tibble() %>%
          dplyr::filter(., stringr::str_detect(param, stringr::str_c("^", x))) %>% 
          dplyr::select(mean, `2.5%`, `97.5%`, Rhat, param, param_fam)
      }) %>% 
        dplyr::bind_cols(parseResFilename(x)) %>% 
        dplyr::rename(q025 = `2.5%`,
                      q975 = `97.5%`) %>% 
        dplyr::mutate(source = x)
    })
    
    if (write) {
      readr::write_csv(param_df, file = param_filename)
    }
  }
  
  return(param_df)
}

getParNames <- function(stanfit_obj) {
  names(stanfit_obj) %>% 
    stringr::str_split(., "\\[") %>% 
    purrr::map_chr(~ .[[1]][1]) %>% 
    unique()
}

#' @title Get time series of predicted probabilities of occurrene
#'
#' @param country
#' @param run_level
#' @param model
#' @param identifier
#'
#' @return df
#' 
getProbs <- function(country,
                     run_level,
                     model = "best",
                     identifier,
                     time_left = NULL,
                     time_right = NULL,
                     gadm_lev = setGadmLev(country, run_level),
                     redo = getRedoDefaults(),
                     verbose = F,
                     ...) {
  
  redo <- completeRedo(redo)
  
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
  
  # Result files to parse
  res_file <- getResFiles(country = country,
                          run_level = run_level, 
                          identifier = identifier, 
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = gadm_lev,
                          what = "output_sim",
                          paths = getResPaths(what = "outputs_sim"),
                          model = model, 
                          ignore = "betas")
  
  ts_probs <- readRDS(res_file)$ts_prob_stats %>% 
    bind_cols(parseResFilename(res_file))
  
  return(ts_probs)
}

#' @title Get posterior retrodictive simulated obsevations
#'
#' @param country
#' @param run_level
#' @param model
#' @param identifier
#'
#' @return df
#' 
getSims <- function(country,
                    run_level,
                    model = "best",
                    identifier,
                    time_left = NULL,
                    time_right = NULL,
                    gadm_lev = setGadmLev(country, run_level),
                    redo = getRedoDefaults(),
                    verbose = F,
                    ...) {
  
  redo <- completeRedo(redo)
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo$data,
                                 verbose = verbose)  
  }
  
  gadm_lev <- setGadmLev(country, run_level)
  
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
  
  data_mapping <- getDataMapping(country = country,
                                 variant = model,
                                 run_lev = run_level,
                                 time_left = time_left,
                                 time_right = time_right,
                                 identifier = identifier)
  
  # Result files to parse
  res_file <- getResFiles(country = country,
                          run_level = run_level, 
                          identifier = identifier, 
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = gadm_lev,
                          what = "output_sim",
                          paths = getResPaths(what = "outputs_sim"),
                          model = model, 
                          ignore = "betas|seas|ts_prob")
  
  sim_stats <- readRDS(res_file)$sim_stats %>% 
    bind_cols(parseResFilename(res_file)) %>% 
    mutate(lev = data_mapping$map_2_lev)
  
  return(sim_stats)
}

#' @title Get time series of predicted probabilities of occurrene
#'
#' @param country
#' @param run_level
#' @param model
#' @param identifier
#'
#' @return df
#' 
getProbSamples <- function(country,
                           run_level,
                           model = "best",
                           identifier,
                           time_left = NULL,
                           time_right = NULL,
                           gadm_lev = setGadmLev(country, run_level),
                           redo = getRedoDefaults(),
                           verbose = F,
                           ...) {
  
  redo <- completeRedo(redo)
  
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
  
  # Result files to parse
  res_file <- getResFiles(country = country,
                          run_level = run_level, 
                          identifier = identifier, 
                          time_left = time_left,
                          time_right = time_right,
                          gadm_lev = gadm_lev,
                          what = "sim",
                          paths = getResPaths(what = "cmdstanr_sim"),
                          model = model)
  
  fit <- readRDS(res_file)
  
  ts_probs <- fit$draws(variable = "ts_prob") 
  
  return(ts_probs)
}

getOffsets <- function(country,
                       run_level,
                       model = "best",
                       identifier,
                       time_left = NULL,
                       time_right = NULL,
                       gadm_lev = setGadmLev(country, run_level),
                       redo = getRedoDefaults(),
                       verbose = F,
                       ...) {
  
  redo <- completeRedo(redo)
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo$data,
                                 verbose = verbose)  
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
  
  
  if (stringr::str_detect(model, "offset")) {
    
    offset_filename <- makeStdResName(country = country,
                                      run_level = run_level,
                                      identifier = identifier,
                                      model = model,
                                      time_left = time_left,
                                      time_righ = time_right,
                                      gadm_lev = gadm_lev,
                                      suffix = "offset")
    
    if (file.exists(offset_filename) & !redo$main) {
      
      res <- readr::read_csv(offset_filename, col_types = readr::cols())
      
      if (nrow(res) > 0){
        return(res)
      }
    } else {
      
      # Result files to parse
      res_files <- getResFiles(country = country,
                               run_level = run_level, 
                               identifier = identifier, 
                               model = model,
                               time_left = time_left,
                               time_righ = time_right,
                               gadm_lev = gadm_lev)
      
      offset_df <- purrr::map_df(res_files, function(x) {
        # Load parsed output file
        res <- readRDS(x)
        
        offset_probs <- tibble::tibble(month = 1:12, 
                                       prob = res$offset_probs)
        
        if (stringr::str_detect(x, "mixture")) {
          
          data_mapping <- getDataMapping(country = country, 
                                         run_lev = run_level, 
                                         identifier = identifier, 
                                         variant = model,
                                         time_left = time_left, 
                                         time_right = time_right)
          
          lambdas <- getParamPosteriors(country = country,
                                        run_level = run_level,
                                        identifier = identifier,
                                        time_left = time_left, 
                                        time_right = time_right,
                                        model = model,
                                        redo = redo,
                                        parameter = "lambdas") %>% 
            dplyr::select(mean, param) %>% 
            dplyr::mutate(unit = extractNum(param),
                          gid = data_mapping$u_unit[unit]) 
          
          offset_df <- lambdas %>% 
            dplyr::mutate(group = ifelse(mean > .5, 1, 2)) %>% 
            dplyr::inner_join(
              rbind(offset_probs %>% dplyr::mutate(group = 2), 
                    offset_probs %>% dplyr::mutate(prob = ifelse(month == 1, 1, 0),
                                                   group = 1))
            ) %>% 
            dplyr::select(-group) %>% 
            dplyr::select(-mean, -param, -unit) %>% 
            dplyr::bind_cols(parseResFilename(x))
          
        } else {
          offset_df <- offset_probs  %>% 
            dplyr::mutate(gid = country) %>% 
            dplyr::bind_cols(parseResFilename(x))
        }
      })
    }
    
  } else {
    offset_df <- tibble::tibble(
      month = NA,
      prob = NA,
      country = country,
      gid = country,
      run_level = run_level,
      variant = NA,
      identifier = identifier,
      time_left = time_left,
      time_right = time_right
    )
  }
  return(offset_df)
}

#' @title Get run data
#'
#' @description gets the data used to run the model
#'
#' @param country what country to get
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param time_left
#' @param time_right
#' @param paths what paths to look into
#' 
#' @details Takes in either the dataframe of model ics or the country, run_level
#' and path arguments to get the results
#' 
#' @return a tibble with the best model name for the given country, run_level and identifier
#' 
getRunDataFile <- function(country,
                           run_level,
                           identifier,
                           time_left = NULL,
                           time_right = NULL,
                           gadm_lev = setGadmLev(country, run_level),
                           paths = getDataPaths()){
  
  data_file <- makeRunDataFile(
    opt = list(cholera_dir = getCholeraDirectory()),
    config = list(
      data = "./",
      country_iso3 = country,
      run_level = run_level,
      thresh = identifier,
      case_thresh = 10,
      start_date = time_left,
      end_date = time_right,
      keep = TRUE,
      gadm_levels = gadm_lev))
  
  return(data_file)
}

makeParamFilename <- function(parameter) {
  upar <- c(stringr::str_subset(parameter, "\\[", negate = T),
            stringr::str_subset(parameter, "\\[") %>% 
              stringr::str_extract("(.)*(?=\\[)") %>% unique())
  
  out <- "_"
  for (par in upar) {
    inds <- stringr::str_subset(parameter , paste0("^", par, "\\[")) %>% 
      stringr::str_extract("(?<=\\[)[0-9]+(?=\\])") %>% 
      stringr::str_c(collapse = "-")
    
    out <- paste0(out, par, "-", inds)
  }
  return(out)
}

#' @title Parse result file name
#'
#' @description extract the country, run level, identifier from result file name
#'
#' @param res_filename
#'
#' @return a tibble with the extraction
#' 
parseResFilename <- function(res_file) {
  country <- stringr::str_extract(res_file, "[A-Z]{3}")
  
  run_level <- stringr::str_extract(res_file, "(?<=runlev)[1-2]{1}") %>% as.numeric()
  
  # Model variant
  variant <- dplyr::case_when(
    stringr::str_detect(res_file, "base") ~ "base",
    stringr::str_detect(res_file, "mixture_offset") ~ "mixture_offset",
    stringr::str_detect(res_file, "offset") ~ "offset",
    stringr::str_detect(res_file, "mixture") ~ "mixture",
    T ~ "null"
  )
  
  # Period
  period <- stringr::str_extract(res_file, paste0("(?<=", variant, "_)(.)*(?=_k)"))
  
  left <- ifelse(stringr::str_detect(period, "^all"), 
                 "all", 
                 stringr::str_sub(period, 1, 10))
  
  right <- ifelse(stringr::str_detect(period, "all$"), 
                  "all", 
                  stringr::str_remove(period, paste0(left, "-")))
  
  # Run identifier
  identifier <- stringr::str_extract(res_file, 
                                     paste0("(?<=", country, "_)",
                                            "(.)*",
                                            "(?=_", variant, ")")) %>% 
    ifelse(is.na(.), "occurrence", .)
  
  identifier <- ifelse(identifier == "10cases", "cases", identifier)
  
  res <- tibble::tibble(country = country,
                        run_level = run_level,
                        variant = variant,
                        identifier = identifier,
                        time_left = left,
                        time_right = right)
  
  return(res)
}


# Postprocessing ----------------------------------------------------------

#' @title Calculate seasonality index 2
#'
#' @description description
#'
#' @param param
#'
#' @return return
#' 
calcSeasIndex2 <- function(probs, 
                           months,
                           peak_month,
                           window_width = 3) {
  
  window_half <- (window_width-1)/2 # half the window width
  
  # The window defining the peak season
  window <- c()
  for (i in -window_half:window_half) {
    m <- peak_month + i
    if (m<1) {
      m <- m + 12
    }
    if (m > 12) {
      m <- m - 12
    }
    window <- c(window, m)
  } 
  
  prob_ratio <- sum(probs[months %in% window])/sum(probs)
  
  return(prob_ratio)
}

#' #' @title Get admin unit mapping
#' #'
#' #' @description Get the mapping from the unique admin units to the GADM units
#' #'
#' #' @param country
#' #' @param admin_lev
#' #'
#' #' @return df
#' #' 
#' getAdminUnitMapping <- function(country,
#'                                 admin_lev, 
#'                                 verbose = F,
#'                                 ...) {
#'   
#'   sf_object <- getShapefile(country = country,
#'                             admin_lev = admin_lev)
#'   
#'   #  Read in generated data
#'   data_mapping_file <- try(getResFiles(country = country,
#'                                        run_level = admin_lev,
#'                                        identifier = "cases",
#'                                        model = "base",
#'                                        paths = "generated_data",
#'                                        what = "data_mapping",
#'                                        verbose = verbose))
#'   
#'   
#'   if (inherits(data_mapping_file, "try-error") | length(data_mapping_file) == 0) {
#'     stop("Failed getting admin mapping because of missing file ", data_mapping_file)
#'   }
#'   
#'   generated_data <- readRDS(data_mapping_file)
#'   
#'   
#'   #  Read in GADM ids
#'   conn <- connectToDB()
#'   id_mapping <- DBI::dbGetQuery(conn, 
#'                                 glue::glue_sql(
#'                                   "SELECT id, gid_{admin_lev} as gid
#'                                 FROM admin.gadm_lev{admin_lev}
#'                                 WHERE id IN ({generated_data$u_unit*});"
#'                                 )) %>% 
#'     dplyr::inner_join(
#'       tibble::tibble(id = generated_data$u_unit,
#'                      serial_id = 1:length(generated_data$u_unit)))
#'   
#'   return(id_mapping)
#' }

#' @title Get best model
#'
#' @description gets the best model based on the ICs
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param paths what paths to look into
#' @param redo whether to redo the extraction
#' 
#' @details Takes in either the dataframe of model ics or the country, run_level
#' and path arguments to get the results
#' 
#' @return a tibble with the best model name for the given country, run_level and identifier
#' 
getBestModel <- function(model_ics = NULL,
                         country = NULL,
                         run_level = NULL,
                         time_left = NULL,
                         time_right = NULL,
                         gadm_lev = setGadmLev(country, run_level),
                         identifier = "occurrence",
                         redo = FALSE,
                         as_df = F,
                         se_diff_thresh = 2,
                         set_diff_tol = .1,
                         ...) {
  # cat(country, run_level, time_left, time_right, identifier, "\n")
  
  if (is.null(model_ics)) {
    if (is.null(country) | is.null(run_level) | is.null(identifier))
      stop("Please provide country, run_level and identifier for model ic extraction")
    
    model_ics <- getModelIC(country = country,
                            run_level = run_level,
                            identifier = identifier,
                            time_left = time_left,
                            time_right = time_right,
                            gadm_lev = gadm_lev,
                            redo = redo)
  }
  
  model_complexity <- 1:5
  names(model_complexity) <- c("null", "base", "offset", "mixture", "mixture_offset")
  
  if (nrow(model_ics) != length(model_complexity))
    warning("Not all models are available, choosing among: ", 
            paste(model_ics$variant, collapse = "; "))
  
  # Get the models that are not significanlty different than the best
  # best_model <- model_ics$variant[model_ics$loo_rank == 1]
  best_model_set <- dplyr::filter(model_ics, elpd_diff == 0 | abs(elpd_diff) < se_diff * (se_diff_thresh - set_diff_tol))
  
  # If multiple models fall within the SE thresh limit choose the simplest
  if (nrow(best_model_set) == 1){
    best_model <- best_model_set$variant[1]
  } else {
    # Get the simplest model
    best_model <- best_model_set %>% 
      # dplyr::filter(variant %in% c(best_model_set$variant[1], "null")) %>%
      dplyr::mutate(complexity = model_complexity[variant]) %>% 
      dplyr::arrange(complexity) %>% 
      dplyr::slice(1) %>% 
      dplyr::pull(variant)
  }
  
  if (as_df) {
    best_model <- tibble::tibble(
      country = country,
      run_level = run_level,
      identifier = identifier,
      period = makePeriod(opt = NULL, time_left, time_right),
      best_model = best_model,
      nmodels = nrow(model_ics),
      models = list(model_ics$variant)
    )
  }
  
  return(best_model)
}



getRedoDefaults <- function() {
  list(main = F,
       best_seas = F,
       best_seas_single = F,
       best_model = F,
       data = F,
       probs = F)
}

#' @title Get best seasonality coefficients
#'
#' @description get best seasonality coefficients for a given run
#'
#' @param param
#'
#' @return return
#' 
getBestSeas <- function(countries,
                        identifier,
                        run_level,
                        redo = getRedoDefaults()) {
  
  redo <- completeRedo(redo)
  
  res_file <- makeStdResName(country = countries,
                             run_level = run_level,
                             identifier = identifier,
                             model = "best",
                             gadm_lev = "auto",
                             suffix = "best_seas")
  
  if (file.exists(res_file) & !redo$main) {
    best_seas <- readr::read_csv(res_file, col_types = readr::cols())
    
  } else {
    if (countries == "all") {
      all_countries <- getAllCountries()
    } else {
      all_countries <- countries
    }
    
    # Get all best seasonality coefficients
    best_seas <- runAll(
      countries = countries,
      identifiers = identifier,
      run_levels = run_level,
      models = "best",
      fun = getBetas,
      fun_name = "best_seas_coefs_admin",
      fun_opts = list(redo = list(main = redo$best_seas_single,
                                  best_model = redo$best_model,
                                  data = redo$data)),
      redo = redo$best_seas
    )
    
    readr::write_csv(best_seas, file = res_file)
  }
  
  return(best_seas)
}

#' @title Get best seasonality coefficients
#'
#' @description get best seasonality coefficients for a given run
#'
#' @param param
#'
#' @return return
#' 
getBestSeasByGroup <- function(countries = "all",
                               identifier = "mean_annual_incidence",
                               run_level = "best",
                               redo = getRedoDefaults()) {
  
  redo <- completeRedo(redo)
  
  res_file <- makeStdResName(country = countries,
                             run_level = run_level,
                             identifier = identifier,
                             model = "best",
                             gadm_lev = "auto",
                             suffix = "best_seas_group")
  
  if (file.exists(res_file) & !redo$main) {
    best_seas <- readr::read_csv(res_file, col_types = readr::cols())
  } else {
    
    if (countries == "all") {
      all_countries <- getAllCountries()
    } else {
      all_countries <- countries
    }
    
    # Get all best seasonality coefficients
    best_seas <- runAll(
      countries = countries,
      identifiers = identifier,
      run_levels = run_level,
      models = "best",
      fun = getParamPosteriors,
      fun_name = "best_seas_coefs_group",
      fun_opts = list(parameter = "betas",
                      redo = list(main = redo$best_seas_single,
                                  data = redo$data,
                                  best_model = redo$best_model)),
      redo = redo$best_seas
    ) %>% 
      # Expand to have all combinations
      tidyr::complete(country = all_countries)
    
    readr::write_csv(best_seas, file = res_file)
  }
  
  return(best_seas)
}

#' @title Get the model information criteria
#'
#' @description Get the model comparison criteria for a set of models
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param model what model to process; one of null, base, mixture, offset, mixture_offset
#' @param identifier what data identifier to process; one of occurrence, cases, mean_annual_incidence
#' @param paths what paths to look into
#' @param redo whether to redo the extraction
#' @param write whether to write the extraction to file for further use
#' 
#' @return a tibble with the information criteria and model weights
getModelIC <- function(country,
                       run_level,
                       identifier = "occurrence",
                       time_left = NULL,
                       time_right = NULL,
                       gadm_lev = setGadmLev(country, run_level),
                       paths = getResPaths(),
                       redo = F,
                       write = T,
                       verbose = F,
                       ...) {
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo_data,
                                 verbose = verbose) 
    
    gadm_lev <- setGadmLev(country, run_level)
  }
  
  ic_filename <- makeStdResName(country = country,
                                run_level = run_level,
                                identifier = identifier,
                                time_left = time_left,
                                time_right = time_right,
                                gadm_lev = gadm_lev,
                                suffix = "model_ics_loo")
  
  if (file.exists(ic_filename) & !redo) {
    return(readr::read_csv(ic_filename, col_types = readr::cols()))
  } else {
    
    # Result files to parse
    res_files <- getResFiles(country = country,
                             run_level = run_level, 
                             identifier = identifier, 
                             time_left = time_left,
                             time_right = time_right,
                             gadm_lev = gadm_lev,
                             paths = paths)
    
    # Get raw loo output and use loo_compare
    loo_objects <- purrr::map(res_files, function(x) {
      res <- readRDS(x)
      res$model_IC$loo}
    )
    
    # Extract info to name loo objects
    loo_models <- purrr::map_df(res_files, ~parseResFilename(.))
    
    # Compare models
    loo_compare_df <- loo::loo_compare(loo_objects) %>% 
      { x <- .
      as.data.frame(x) %>%
        mutate(id = rownames(x) %>% str_extract("[0-9]") %>% as.numeric())
      } %>% 
      tibble::as_tibble() %>% 
      dplyr::inner_join(loo_models %>% mutate(id = row_number())) %>% 
      dplyr::select(-id)
    
    if (write) {
      # Write to file
      readr::write_csv(loo_compare_df, file = ic_filename)
    }
  }
  
  return(loo_compare_df)
}


# RunAll calls ------------------------------------------------------------

getAllPeakMonths <- function(countries = "all",
                             identifier,
                             run_level,
                             best_seas = NULL,
                             redo = getRedoDefaults(),
                             redo_single = F) {
  
  redo <- completeRedo(redo)
  
  if (countries == "all") {
    all_countries <- getAllCountries()
  } else {
    all_countries <- countries
  }
  
  res_file <- makeStdResName(country = paste(countries, collapse = "-"),
                             run_level = run_level,
                             identifier = identifier,
                             model = "best",
                             time_left = "all",
                             time_right = "all",
                             suffix = "peak_months",
                             file_type = "csv")
  
  if (!file.exists(res_file) || redo$main) {
    
    
    if (is.null(best_seas)) {
      best_seas <- getBestSeas(countries = countries,
                               identifier = identifier,
                               run_level = run_level,
                               redo = list(main = redo$best_seas,
                                           best_seas_single = redo$best_seas_single))
    }
    
    peak_month <- best_seas %>% 
      dplyr::group_by(country, gid, variant, run_level) %>% 
      dplyr::summarise(peak_month = which.max(mean),
                       peak_val = max(mean)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(country) %>% 
      dplyr::group_map(function(x, y){
        if (is.na(x$variant[1])) {
          x <- x %>% 
            dplyr::mutate(gid = y$country[1], 
                          country = y$country[1],
                          variant = x$variant[1]) %>% 
            dplyr::distinct()
        } else {
          x <- x %>% 
            dplyr::mutate(country = y$country[1],
                          variant = x$variant[1])
        }
        return(x)
      }) %>% 
      dplyr::bind_rows()
    
    readr::write_csv(peak_month, file = res_file)
    
  } else {
    peak_month <- readr::read_csv(file = res_file)
  }
  
  return(peak_month)
}

getAllBestModels <- function(countries = "all",
                             run_levels = "best",
                             identifier,
                             redo = F,
                             redo_single = F,
                             verbose = F,
                             error_handling = "remove"){
  
  runAll(countries = countries,
         run_levels = run_levels,
         identifiers = identifier,
         models = "all",
         fun = getBestModel,
         fun_name = "best_model",
         fun_opts = list(as_df = T,
                         model_ics = NULL,
                         redo = redo_single),
         error_handling = error_handling,
         redo = redo,
         verbose = verbose)
}


computeAllSeasIndex <- function(countries = "all",
                                run_levels = "best",
                                identifier,
                                redo = getRedoDefaults(),
                                verbose = F,
                                error_handling = "remove") {
  
  redo <- completeRedo(redo)
  
  seasind_file <- makeStdResName(country = countries,
                                 run_level = run_levels,
                                 identifier = identifier,
                                 model = "best",
                                 gadm_lev = "auto",
                                 suffix = "seas_index",
                                 file_type = "rds")
  
  if (file.exists(seasind_file) & ! redo$main) {
    seas_index <- readRDS(seasind_file)
  } else {
    
    # Get time series of probabilities
    probs <- runAll(countries = countries,
                    run_levels = run_levels,
                    identifiers = identifier,
                    models = "best",
                    fun = getProbs,
                    fun_name = "probs",
                    redo = redo$probs,
                    fun_opts = list(
                      redo = list(
                        best_model = redo$best_model)
                    ),
                    error_handling = error_handling) %>% 
      filter(variant != "null", !is.na(mean)) %>% 
      mutate(year = factor(year)) %>% 
      rename(gid = unit)
    
    # Get peak month for each admin unit
    peak_month <- getAllPeakMonths(countries = "all",
                                   identifier = identifier,
                                   run_level = run_levels,
                                   redo = list(main = redo$peak_month,
                                               best_seas = F,
                                               best_seas_single = F))
    # Join 
    probs <-inner_join(probs, peak_month)
    
    # Compute seasonality indices
    seas_index <- probs %>% 
      dplyr::group_by(country, gid, variant) %>% 
      dplyr::summarise(seas_index = calcSeasIndex2(mean, month, peak_month)) %>% 
      dplyr::ungroup()
    
    # Get sf objects
    sf_objects <- getSfObjects(probs)
    
    seas_index <- sf_objects %>% 
      dplyr::inner_join(seas_index)
    
    saveRDS(seas_index, file = seasind_file)
  }
  seas_index
}


#' @title Run all
#'
#' @description Runs a function over a set of combinations of country,
#' run level, and, identifier
#'
#' @param countries
#' @param run_levels
#' @param identifiers
#' @param models
#' @param times_left
#' @param times_right
#' @param fun function to run, must take in arguments country, run_level,
#' identifier, model, time_left and time_right
#' @param redo
#'
#' @return a dataframe
#' 
runAll <- function(countries = "all",
                   run_levels = "all",
                   identifiers = "all",
                   models = "all",
                   times_left = "all",
                   times_right = "all",
                   gadm_lev = "auto",
                   fun,
                   fun_name,
                   fun_opts = NULL,
                   postfun = NULL,
                   error_handling = "remove",
                   redo = F,
                   ...) {
  
  res_file <- makeStdResName(country = paste(countries, collapse = "-"),
                             run_level = paste(run_levels, collapse = "-"),
                             identifier = paste(identifiers, collapse = "-"),
                             model = paste(models, collapse = "-"),
                             time_left = paste(times_left, collapse = "-"),
                             time_right = paste(times_right, collapse = "-"),
                             gadm_lev = gadm_lev,
                             suffix = fun_name,
                             file_type = "rds")
  
  if (file.exists(res_file) & !redo) {
    all_res <- readRDS(res_file)
  }  else {
    if (countries == "all") {
      countries <- getAllCountries()
    }
    
    if (run_levels == "all") {
      run_levels <- getAllRunLevels()
    }
    
    if (identifiers == "all") {
      identifiers <- getAllIdentifiers()
    }
    
    # Parallel setup
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    
    export_packages <- c("tidyverse", "magrittr", "foreach", "rstan",
                         "lubridate", "sf")
    
    all_res <- foreach(cntry = countries,
                       .combine = bind_rows,
                       .errorhandling = error_handling,
                       .packages = export_packages) %:% 
      foreach(id = identifiers,
              .combine = bind_rows,
              .errorhandling = error_handling,
              .packages = export_packages) %:%
      foreach(rl = run_levels,
              .combine = bind_rows,
              .errorhandling = error_handling,
              .packages = export_packages) %:%
      foreach(mdl = models,
              .combine = bind_rows,
              .errorhandling = error_handling,
              .packages = export_packages) %:%
      foreach(tl = times_left,
              tr = times_right,
              .combine = bind_rows,
              .errorhandling = error_handling,
              .packages = export_packages) %do% { 
                cat("Running", cntry, id, rl, mdl, "\n")
                
                source(str_c(getCholeraDirectory(), "/analysis/utils.R"))
                source(str_c(getCholeraDirectory(), "analysis/modeling_utils.R"))
                source(str_c(getCholeraDirectory(), "analysis/postprocessing_utils.R"))
                
                args <- c(
                  list(
                    country = cntry,
                    run_level = rl,
                    model = mdl,
                    identifier = id,
                    time_left = tl,
                    time_right = tr),
                  fun_opts)
                
                res <- do.call(fun, args)
                
                # res <- try(do.call(fun, args))
                # 
                # if (!inherits(res, "try-error")) {
                #   
                #   if (!is.null(postfun)) {
                #     res <- postfun(res)
                #   }
                #   return(res)  
                # }
              }
    
    parallel::stopCluster(cl)
    
    saveRDS(all_res, file = res_file)
  }
  return(all_res)
}

# Covariates --------------------------------------------------------------

#' @title Get climatology
#'
#' @description Get the climatology data for a given country at a given admin level
#'
#' @param country
#' @param admin_level
#' 
#' @return a tibble with the climatology data
#' 
getClimatology <- function(country,
                           admin_level,
                           covar = "all") {
  
  conn <- connectToDB()
  
  all_covar <- DBI::dbGetQuery(conn, "SELECT DISTINCT covar FROM climatology;")$covar
  
  if (covar == "all") {
    covar <- all_covar
  } else {
    if (!any(purrr::map_lgl(covar, ~ . %in% all_covar))) {
      stop("Covariate ", covar, " not found")
    }
  }
  
  if (country == "all") {
    country <- getAllCountries()
  }
  
  clim_data <- DBI::dbGetQuery(conn, glue::glue_sql(
    "SELECT a.gadm_lev, a.gid, a.gadm_id, a.covar, a.month, a.value, a.source
    FROM climatology a
    JOIN 
    (SELECT id, gid_0 as gid, 0 as gadm_lev 
    FROM admin.gadm_lev0
    WHERE gid_0 IN ({country*}) 
    UNION
    SELECT id, gid_1 as gid, 1 as gadm_lev 
    FROM admin.gadm_lev1
    WHERE gid_0 IN ({country*}) 
    UNION
    SELECT id, gid_2 as gid, 2 as gadm_lev 
    FROM admin.gadm_lev2
    WHERE gid_0 IN ({country*})) b
    ON a.gadm_id = b.id AND a.gadm_lev = b.gadm_lev
    WHERE a.gadm_lev IN ({admin_level*})
    AND a.covar IN ({covar*})",
    .con = conn
  )) %>% 
    tibble::as_tibble() 
  
  DBI::dbDisconnect(conn)
  
  return(clim_data)
}


#' @title Get time covar
#'
#' @description Get the time varying covariate data for a given country at a given admin level
#'
#' @param country
#' @param admin_level
#' 
#' @return a tibble with the time varying covariate data
#' 
getTimeCovar <- function(country,
                         admin_level,
                         covar = "all") {
  
  conn <- connectToDB()
  
  all_covar <- getAllCovar(type = "time_covar")$covar  
  
  if (covar == "all") {
    covar <- all_covar
  } else {
    if (!any(purrr::map_lgl(covar, ~ . %in% all_covar))) {
      stop("Covariate ", covar, " not found")
    }
  }
  
  if (country == "all") {
    country <- getAllCountries()
  }
  
  tc_data <- DBI::dbGetQuery(conn, glue::glue_sql(
    "SELECT a.gadm_lev, a.gid, a.gadm_id, a.covar, a.month, a.year, a.value
    FROM time_covar a
    JOIN 
    (SELECT id, gid_0 as gid, 0 as gadm_lev 
    FROM admin.gadm_lev0
    WHERE gid_0 IN ({country*}) 
    UNION
    SELECT id, gid_1 as gid, 1 as gadm_lev 
    FROM admin.gadm_lev1
    WHERE gid_0 IN ({country*}) 
    UNION
    SELECT id, gid_2 as gid, 2 as gadm_lev 
    FROM admin.gadm_lev2
    WHERE gid_0 IN ({country*})) b
    ON a.gadm_id = b.id AND a.gadm_lev = b.gadm_lev
    WHERE a.gadm_lev IN ({admin_level*})
    AND a.covar IN ({covar*})",
    .con = conn
  )) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(date = as.Date(stringr::str_c(year, month, 1, sep = "-")))
  
  DBI::dbDisconnect(conn)
  
  return(tc_data)
}

#' @title Ingest Climatology
#'
#' @description Ingests climatology data
#'
#' @param clim_files climatology files to ingest
#' @param overwrite
#' @return NULL
#' 
ingestClimatology <- function(files = c("data/GADM0_monthlyclimatology_19902018.nc",
                                        "data/GADM1_monthlyclimatology_19902018.nc",
                                        "data/GADM2_monthlyclimatology_19902018.nc"),
                              overwrite = T) {
  # Connect to db
  conn <- connectToDB()
  
  if (overwrite) {
    DBI::dbSendStatement(conn, "DROP TABLE IF EXISTS climatology;")
  }
  
  for (f in files) {
    # Get the gadm level
    lev <- stringr::str_extract(f, "[0-2]{1}") %>% as.numeric()
    # Get the corresponding lookup table
    lookup <- readr::read_csv(paste0("data/lookup_gadm", lev, ".csv"), 
                              col_types = readr::cols())
    
    dat <- tidync::tidync(f) %>%
      tidync::hyper_tibble() %>% 
      dplyr::rename(uid = GADM) %>% 
      dplyr::inner_join(lookup %>% 
                          magrittr::set_colnames(c("uid", "count", "gid")), 
                        by = "uid") %>% 
      dplyr::mutate(gadm_lev = lev) %>% 
      dplyr::select(-count, -uid) %>% 
      tidyr::pivot_longer(cols = !c("month", "gid", "gadm_lev"),
                          names_to = "covar",
                          values_to = "value") %>% 
      dplyr::mutate(source = f) %>% 
      dplyr::inner_join(DBI::dbGetQuery(
        conn, 
        glue::glue_sql("SELECT id as gadm_id, gid_{lev} as gid
                       FROM admin.gadm_lev{lev};",
                       .conn = conn)))
    
    # Append file to climatology table
    DBI::dbWriteTable(conn, name = "climatology", value = dat, append = T)
  }
  
  DBI::dbSendStatement(conn, "CREATE INDEX tc_covar ON time_covar (covar);")
  DBI::dbSendStatement(conn, "CREATE INDEX tc_gadmid ON time_covar (gadm_id);")
  DBI::dbSendStatement(conn, "CREATE INDEX tc_gadmlev ON time_covar (gadm_lev);")
  
  DBI::dbSendStatement(conn, "VACUUM time_covar;")
  
  # Disconnect
  DBI::dbDisconnect(conn)
}

#' @title Ingest time covar
#'
#' @description Ingest time varying environmental covariate
#'
#' @param files
#' @param overwrite
#'
#' @return NULL
#' 
ingestTimeCovar <- function(files = c("data/GADM0_FloodScan_20002020.nc",
                                      "data/GADM1_FloodScan_20002020.nc"),
                            overwrite = T) {
  # Connect to db
  conn <- connectToDB()
  
  if (overwrite) {
    DBI::dbSendStatement(conn, "DROP TABLE IF EXISTS time_covar;")
  }
  
  for (f in files) {
    
    # Get the gadm level
    lev <- stringr::str_extract(f, "[0-2]{1}") %>% as.numeric()
    # Get the corresponding lookup table
    lookup <- readr::read_csv(paste0("data/lookup_gadm", lev, ".csv"), 
                              col_types = readr::cols())
    
    cat("-- Ingesting", f, "\n")
    
    dat <- tidync::tidync(f) %>%
      tidync::hyper_tibble() %>% 
      {
        if ("ID" %in% colnames(.)) {
          dplyr::rename(., uid = ID)
        } else if ("GADM" %in% colnames(.)){
          dplyr::rename(., uid = GADM)
        }
      } %>% 
      dplyr::inner_join(lookup %>% 
                          magrittr::set_colnames(c("uid", "count", "gid")), 
                        by = "uid") %>% 
      dplyr::mutate(gadm_lev = lev) %>% 
      dplyr::select(-count, -uid) %>% 
      tidyr::pivot_longer(cols = !c("month", "year", "gid", "gadm_lev"),
                          names_to = "covar",
                          values_to = "value") %>% 
      dplyr::mutate(source = f) %>% 
      dplyr::inner_join(DBI::dbGetQuery(
        conn, 
        glue::glue_sql("SELECT id as gadm_id, gid_{lev} as gid
                       FROM admin.gadm_lev{lev};",
                       .conn = conn)))
    
    # Append file to time_covar table
    DBI::dbWriteTable(conn, 
                      name = "time_covar", 
                      value = dat, 
                      append = T)
    
    # Append monthly means to climatology table
    dat %>% 
      dplyr::group_by(month, gid, gadm_id, gadm_lev, covar) %>%
      dplyr::summarise(value = mean(value)) %>% 
      mutate(source = f) %>% 
      DBI::dbWriteTable(conn, 
                        name = "climatology", 
                        value = ., 
                        append = T)
    
    cat("-- Done ingesting", f, "\n")
  }
  
  # Disconnect
  DBI::dbDisconnect(conn)
}

#' @title Make Metadata
#'
#' @description makes metadata table
#'
#' @param redo
#'
#' @return none
makeMetadata <- function(overwrite = F) {
  conn <- connectToDB()
  
  if (overwrite) {
    DBI::dbSendStatement(conn, "DROP TABLE IF EXISTS covar_metadata;")
  }
  
  DBI::dbSendStatement(conn, "CREATE TABLE covar_metadata AS (
  WITH tmp AS (SELECT covar, gadm_lev, source,
  MIN(year) as year_left, 
  MAX(year) as year_right,
  'time_covar' as type
  FROM time_covar
  GROUP BY covar, source, gadm_lev
  UNION
  SELECT DISTINCT covar, gadm_lev, source, 1990 as year_left, 2018 as year_right,
  'climatology' as type
  FROM climatology WHERE source LIKE '%climatology%')
                  SELECT * FROM tmp ORDER BY covar, source, gadm_lev;)")
  
  
  DBI::dbDisconnect(conn)
}

# Wrappers ----------------------------------------------------------------

# Functions in this section depend on utility functions in other scripts

getAdminUnitMappingWrapper <-  function(country,
                                        identifier,
                                        run_level,
                                        time_left = NULL,
                                        time_right = NULL,
                                        verbose = F,
                                        ...) {
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 verbose = verbose)  
  }
  
  getAdminUnitMapping(country = country,
                      admin_lev = run_level, 
                      verbose = verbose) 
}


#' @title Get run data wrapper
#'
#' @description wrapper around the call to getRunData
#'
#' @param country
#' @param identifier
#' @param time_left
#' @param time_right
#'
#' @return data frame with run data
#' 
getRunDataWrapper <- function(country, 
                              identifier, 
                              time_left = NULL, 
                              time_right = NULL,
                              redo = F,
                              verbose = T,
                              ...) {
  
  # Get data files
  data_file <-  getRunDataFile(country = country,
                               run_level = 2, 
                               identifier = identifier,
                               time_left = time_left,
                               time_right = time_right)
  
  if (file.exists(data_file) & !redo) {
    run_data <- readRDS(data_file)
  } else {
    if (verbose) {
      cat("-- Re-pulling data from file ",  getMonthlyDataPath(), "\n")
    }
    run_data <- getRunData(path_to_monthly_data = getMonthlyDataPath(),
                           country_iso3 = country,
                           start_date = time_left,
                           end_date = time_right,
                           thresh = identifier)
    
    saveRDS(run_data, file = data_file)
  }
  
  return(run_data)
}

#' @title Wrapper get time covar
#'
#' @description Wrapper to be able to call runall
#'
#' @param param
#'
#' @return return
#' 
getTimeCovarWrapper <- function(country,
                                identifier,
                                run_level,
                                time_left = NULL,
                                time_right = NULL,
                                redo_data = F,
                                verbose = F,
                                covar = "all",
                                ...) {
  
  if (run_level == "best") {
    run_level <- getBestRunLevel(country = country,
                                 identifier = identifier,
                                 time_left = time_left,
                                 time_right = time_right,
                                 redo = redo_data,
                                 verbose = verbose)  
  }
  
  getTimeCovar(country = country,
               admin_level = run_level,
               covar = covar)
}

# Misc --------------------------------------------------------------------

#' @title Check all models done
#'
#' @description Checks whether all models were computed
#'
#' @param country
#' @param run_level
#' @param identifier  
#'
#' @return a logical whether all models were done

checkAllDone <- function(country,
                         run_level,
                         identifier,
                         time_left = NULL,
                         time_right = NULL,
                         suffix,
                         file_type) {
  
  sample_files <- purrr::map_chr(getAllModels(),
                                 ~makeStdResName(country = country, 
                                                 run_level = run_level,
                                                 identifier = identifier, 
                                                 model = .,
                                                 time_left = time_left,
                                                 time_right = time_right,
                                                 suffix = suffix,
                                                 file_type = file_type))
  
  res_files <- purrr::map_chr(getAllModels(),
                              function(x) {
                                res <- getResFiles(country = country, 
                                                   run_level = run_level,
                                                   identifier = identifier, 
                                                   model = x,
                                                   time_left = time_left,
                                                   time_right = time_right)
                                if (length(res) == 0){
                                  NA
                                } else {
                                  res
                                }
                              })
  
  all_done <- all(purrr::map2_lgl(sample_files, res_files, ~ (checkFile(.x) | is.na(.y))))
  
  return(all_done)
}

checkFile <- function(x) {
  if (file.exists(x)) {
    res <- readRDS(x)
    if (!is.null(res)) {
      TRUE
    } else {
      FALSE
    }
  } else {
    FALSE
  }
}


#' @title Compute yearly data avail
#'
#' @description Computes the number of years with available monthly-submonthly dat
#'
#' @param obs_len_thresh
#'
#' @return tibble
#' 
computeNumYearValid <-  function(df = NULL,
                                 country = NULL,
                                 identifier, 
                                 time_left = NULL, 
                                 time_right = NULL,
                                 redo = F,
                                 obs_len_thresh = 1,
                                 ...) {
  
  if (is.null(df)) {
    df <- getRunDataWrapper(country = country,
                            identifier = identifier,
                            time_left = time_left,
                            time_right = time_right)
  }
  
  df %>% 
    dplyr::filter(map_dbl(months, length) <= obs_len_thresh) %>% 
    count(year, gadm_lev) %>% 
    group_by(gadm_lev) %>% 
    summarise(n = n()) %>% 
    mutate(country = country)
}

completeRedo <- function(redo) {
  rdefault <- getRedoDefaults()
  
  for (i in names(rdefault)){
    if (is.null(redo[[i]])) {
      redo[[i]] <- rdefault[[i]]
    }
  }
  redo
}

#' @title Extract num
#'
#' @description Extract a number from the parameter name
#'
#' @param str
#'
#' @return double
#' 
extractNum <- function(str) {
  as.numeric(stringr::str_extract(str, "(?<=\\[)[0-9]+"))
}

#' @title Make standardized file
#'
#' @description Makes a standardized file path for result extraction
#'
#' @param country what country to get
#' @param run_level what run levels to get
#' @param suffix additional optional suffix
#' @param res_path the path to the results, defaults to generated_data/model_outputs/
#' @param verbose say where things are saved
#' 
#' @return a string with the results file name
#' 
makeStdResName <- function(country,
                           run_level,
                           identifier,
                           model = "all",
                           time_left = NULL,
                           time_right = NULL,
                           gadm_lev = setGadmLev(country, run_level),
                           suffix,
                           file_type = "csv",
                           res_path = str_c(getCholeraDirectory(), "/generated_data/postprocessing/"),
                           verbose = F) {
  
  if (!dir.exists(res_path)) {
    dir.create(res_path)
  }
  
  period <- makePeriod(opt = NULL, time_left, time_right)
  
  if (period == "all-all") {
    filename <- paste0(res_path, 
                       paste(country, model, identifier, 
                             paste0("runlev", run_level), gadm_lev, suffix, sep = "_"), ".", file_type)
  } else {
    filename <- paste0(res_path, 
                       paste(country, model, identifier, 
                             paste0("runlev", run_level), 
                             period, gadm_lev, suffix, sep = "_"), ".", file_type)
  }
  
  return(filename)
}

#' @title Parse columns
#'
#' @description Parse columns name
#'
#' @param col
#'
#' @return string
parseCol <- function(col) {
  str_c("('", str_c(col, collapse = "' '"), "')")
}


#'
#' @description checks whether requested parameters are allowed
#'
#' @param country
#' @param run_level
#' @param identifier
#' @param model
#'
#' @return none
runChecks <- function(country,
                      run_level,
                      identifier,
                      model) {
  
  if (!(country %in% getAllCountries())) {
    stop("Country ", country, " not allowed, needs to be one of {'",
         paste(getAllCountries(), collapse = "', '"), "'}")
  }
  
  if (!(run_level %in% c(1, 2))) {
    stop("Run level", run_level, " not allowed, needs to be one of {'",
         paste(getAllRunLevels(), collapse = "', '"), "'}")
  }
  
  if (!(model %in% getAllModels() | model == "all")) {
    stop("Model ", model, " not known, needs to be one of {'", 
         paste(getAllModels(), collapse = "', '"), "'}")
  }
  
  if (!(identifier %in% getAllIdentifiers() | identifier == "all")) {
    stop("Identifier ", identifier, " not known, needs to be one of {'",
         paste(getAllIdentifiers(), collapse = "', '"), "'}")
  }
  
}

setGadmLev <- function(country, run_level) {
  ifelse(country %in% c("COD", "NGA") & run_level == 2, "12", "012")
}

