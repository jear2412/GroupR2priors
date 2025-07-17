# Misc  ----------------------------------------------------

#used in brms 
bind <- function(...) cbind(...)

# Plot distributions ----------------------------------------------------

#------ Plot distributions
plot_dist <- function(dist, bounds, pars, xtype = c("c", "d"), xname="", 
                      prefix = c("d", "p", "q"), parnames = NULL, 
                      package = NULL, ...) {
  xtype <- match.arg(xtype)
  prefix <- match.arg(prefix)
  pos <- -1
  if (!is.null(package)) {
    pos <- asNamespace(package)
  }
  dist_fun <- get(paste0(prefix, dist), pos = pos, mode = "function")
  if (xtype == "c") {
    # continuous
    df <- data.frame(x = seq(bounds[1], bounds[2], 0.001))
  } else if (xtype == "d") {
    # discrete
    df <- data.frame(x = bounds[1]:bounds[2])
  }
  if (!is.null(parnames)) {
    parnames <- paste0(parnames, " = ")
  }
  cnames <- rep(NA, length(pars))
  for (i in seq_along(pars)) {
    tmp <- do.call(dist_fun, c(list(df$x), pars[[i]], list(...)))
    cnames[i] <- paste0("$", parnames, pars[[i]], "$", collapse = ", ")
    df[paste0(parnames, pars[[i]], collapse = ", ")] <- tmp
  }
  
  df <- df %>%
    gather("pars", "dens", -x) %>%
    mutate(pars = factor(pars, unique(pars)))
  
  gg <- ggplot(df, aes(x, dens, color = pars))
  if (xtype == "c") {
    gg <- gg + geom_line(size=1.5)
  } else if (xtype == "d") {
    gg <- gg + 
      geom_linerange(aes(ymin=0, ymax=dens), size = 1) +
      geom_line(size = 0.8, linetype = "dotted", alpha = 0.8)
  }
  
  gg <- gg + 
    #scale_colour_manual(values=my_colors,labels = unname(latex2exp::TeX(cnames)))+
    scale_colour_viridis_d(labels = unname(latex2exp::TeX(cnames)))+
    labs(x = xname, y = "", color = "") + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "bottom",
      axis.text.x=element_text(size=20),
      legend.text = element_text(size = 20)
    )
  if (prefix == "p") {
    gg <- gg +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  } else if (prefix == "q") {
    gg <- gg +
      scale_y_continuous() +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  }
  gg
}

mytdist <- function( x, nu=1, mean=0, sd=1, log=FALSE){
  return(dnorm(x, mean= mean, sd=sd)^nu)
}

trim_vec <- function(x, lo, hi = lo){
  if(lo + hi > length(x)){
    stop("cannot remove more than 'length(x)' elements, 'lo + hi > length(x)'")
  }
  i <- order(x)
  x <- x[i]
  x <- x[-seq.int(to = length(x), length.out = hi)]
  x <- x[-seq.int(to = lo)]
  x <- x[order(i)]
  x[!is.na(x)]
}


# DivergenceFunctions ----------------------------------------------------

kl.diri.logitnormal <- function(alpha, ref = 1){
  
  p <- length(alpha)
  mu <- digamma(alpha)-digamma(alpha[ref])
  #mu <- mu[-ref]
  Sigma <- matrix( trigamma(alpha[ref]),p,p)
  diag(Sigma) <- trigamma(alpha)+trigamma(alpha[ref])
  
  list(alpha=alpha, mu=mu, Sigma=Sigma)  
  
}

# Matrix Functions ----------------------------------------------------

tridiag <- function(upper, lower, main){
  out <- matrix(0,length(main),length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx+1,indx)] <- lower
  out[cbind(indx,indx+1)] <- upper
  return(out)
}


# Projpred functions ----

predfun <- function(fit, transform = FALSE, newdata = NULL) {
  if (!is.null(newdata)) {
    X <- newdata$X
  } else {
    X <- fit$data$X
  }
  draws <- posterior::as_draws_matrix(fit$fit$draws(c("b_Intercept", "beta")))
  b_Intercept <- as.numeric(draws[, 1])
  b <- draws[, 2:ncol(draws)]
  t(b_Intercept + b %*% t(X))
  
}


# Simulation Functions ----------------------------------------------------


dataset_cond_sim <- function(sim_cond,
                             sim_params,
                             seed= NULL,
                             path= NULL,
                             ncores= 1,
                             exp_name= NULL){
  
  
  # Set seed for reproducibility.
  # global seed to control the simulations
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  seed_list <- sample(1000000000:.Machine$integer.max,
                      size = sim_params$nsims) #seeds used for each simulation
  
  if (file.exists(paste0(paste(path, paste0(exp_name,"_",sim_cond$id), sep = "/"), ".RDS"))) {
    return(readRDS(paste0(paste(path, paste0(exp_name,"_",sim_cond$id), sep = "/"), ".RDS")))
  } else {
    if (ncores > 1) {
      # Multiprocessing setup
      
      cluster <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cluster)
      
      parallel::clusterExport(cl = cluster, c("sim_params"))
      
      parallel::clusterEvalQ(cl = cluster, {
        library(brms)
        library(tidybayes)
        library(cmdstanr)
        library(projpred)
        library(posterior)
        library(dplyr)
        library(tibble)
        library(dplyr)
        library(mvtnorm)
        library(LaplacesDemon)
        library(parallel)
        library(Matrix)
        library(ppcor)
        library(gtools)
        
        
        # scripts
        source("scripts/aux_functions/all_auxfunctions.R") #script containing auxiliary functions
        source("scripts/aux_functions/coef_functions.R") #script containing coefficients auxiliary functions
        source("scripts/aux_functions/cov_functions.R") #script containing covariance auxiliary functions
        source("scripts/aux_functions/dgp_functions.R") #script containing dgp auxiliary functions
        source("scripts/aux_functions/summary_functions.R") #script containing summary auxiliary functions
        source("scripts/aux_functions/fit_params.R") #mcmc related functions 
        source("scripts/aux_functions/fit_params_data.R") #function to create mcmc params data
        source("scripts/aux_functions/stan_fits.R") #stan fits 
        source("scripts/aux_functions/full_sim_comparison.R") #special function to run the simulation
        
        source("scripts/aux_functions/R2_alpha_gen.R") #how to generate concentration vector from r2d2

        cmdstanr::set_cmdstan_path(sim_params$cmdstan_path)
        #options(matrixStats.vars.formula.freq = 2e6 ) # Weird bug ?
        
    
      })
      
      `%dopar%` <- foreach::`%dopar%`
      
      # Multiprocessing run
      # Parallelizing over number of simulations (represented by number of seeds)
      results <- foreach::foreach(
        par_seed = seed_list
      ) %dopar% {
        cond_sim(
          sim_params = sim_params,
          sim_cond = sim_cond,
          seed = par_seed
        )
      }
      
      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    }else { # Single processing
      results <- vector("list", length = length(seed_list))
      par_seed = seed_list
      for (i in seq_along(results)) {
        results[[i]] <- cond_sim(
          sim_params = sim_params,
          sim_cond = sim_cond,
          seed = par_seed[i]
        )
      }
    } 
    
    final_result <- do.call(rbind, results)
    
    if (!is.null(path)) {
      saveRDS(final_result, paste0(paste(path, 
                                         paste0(exp_name,"_",sim_cond$id), 
                                         sep = "/"), ".RDS"))
    }
    return(final_result) 
  }
}

#Simulate data and summarise a fit given a simulation condition
cond_sim <-  function(sim_params, sim_cond, seed = NULL, verbose = TRUE){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #------ Generate data
  # Include seed for the data-generating process (dgp)
  
  sim_cond$seed <- seed  #need the seed to control the dgp!

  dgp_list <- tryCatch({
    make_dgp(list(sim_cond = sim_cond))
  }, error = function(e) {
    stop("Error generating data: ", e)
  })
  
  #------  Fit different models
  # Generate parameters of fits
  fits_params <-  get_sim_fit_params_data(fit_params = sim_params$fit_params,
                                           sim_cond = sim_cond)
  
  fit_options <- list( nchains = sim_params$nchains, 
                       iter_warmup = sim_params$iter_warmup, 
                       iter_sample = sim_params$iter_sample, 
                       prior_only = sim_params$prior_only)
  
  
  summary_list <- vector(mode = "list", length = length(fits_params$ids_list)  )
  
  for(i in 1: length(fits_params$ids_list) ){
    
    if (verbose) message("------- Fitting model: ", fits_params$ids_list[i], " -------")
   
    tryCatch({
      
      fit_mcmc_params <- make_fit_params(fits_params$names_list[i], 
                                       fits_params$fit_params[[i]], 
                                       dgp_list)
    
      params_list <- c(dgp_list,
                     fit_mcmc_params, 
                     fit_options)
    
      # Fit the model using the corresponding fitting function
      fitfn <- get(paste0(fits_params$fits_list[i],"fit")) # Get fit function
      fit <- fitfn(params_list) # Run stan fit
      
      #Summarize
      if(!fit$fit.error){
        
        fit_summary_params <- list(
          fit= fit$fit,
          real_theta= dgp_list$real_theta,
          standat= params_list,
          voi= sim_params$smqoi$voi,
          moi= sim_params$smqoi$moi,
          probsoi= sim_params$smqoi$probsoi, 
          proj_pred_flag = sim_cond$proj_pred_flag, 
          seed= seed)
        
        summary_list[[i]] <- fit_summary(fit_summary_params)
      }else{
        message("Warning: Fitting error for model ", fits_params$ids_list[i])
      }
    }, error = function(e){
      message("Error occurred for model ", fits_params$ids_list[i], ": ", e$message)
  
    })
  }
  
  names(summary_list) <- fits_params$ids_list
  return(summary_list)
  
}

cond_sim_nsims <-  function(sim_params, sim_cond, seed = NULL){
  
  # simulation condition over number of simulations
  # data is pregenerated
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #------ Read data
  sim_cond$seed <- seed #need the seed to control the process
  dataset_path <- paste0(sim_params$sim_cond_directory, "/", sim_params$exp_name, "_", sim_params$global_seed, "_", sim_cond$id
                         , "_", seed )
  
  dgp_list <- readRDS(paste0(dataset_path, ".RDS"))
  #------  Fit different models
  # Generate parameters of fits
  
  fits_params <-  get_sim_mcmc_params_data(mcmc_params = sim_params$mcmc_params,
                                           gen_coef = sim_cond$gencoef_function,
                                           sim_cond = sim_cond)
  
  nfits <- fits_params$nfits #number of fits
  fits_list <- fits_params$fits_list #fits to considered
  
  ids_list <- fits_params$ids_list # ids assigned to different models 
  names_list <- fits_params$names_list # names of the stan fits
  nids <- length(ids_list) 
  
  summary_list <- vector(mode = "list", length = nids)
  params_list <- vector(mode = "list", length = nids)
  
  #Cycle through models
  for(i in 1:nids){
    
    # MCMC params of a specific fit
    fit_mcmc_params <- make_mcmc_params(names_list[i], 
                                        fits_params$mcmc_params[[i]], 
                                        dgp_list)
    
    fit_mcmc_params$nchains <- sim_params$nchains
    fit_mcmc_params$iter_warmup <- sim_params$iter_warmup
    fit_mcmc_params$iter_sample <- sim_params$iter_sample
    
    params_list[[i]] <- append(mcmc_params,
                               dgp_list)
    
    
    fitfn <- get(paste0(fits_list[i],"fit")) # Get fit function
    fit <- fitfn(params_list[[i]]) # Run stan fit
    
    #Summarize
    #TODO: change
    if(!fit$fit.error){
      fit_summary_params <- list(fit= fit$fit,
                                 real_theta= dgp_list$real_theta,
                                 standat= params_list[[i]],
                                 voi= sim_params$smqoi$voi,
                                 moi= sim_params$smqoi$moi,
                                 probsoi= sim_params$smqoi$probsoi, 
                                 proj_pred_flag = sim_cond$proj_pred_flag, 
                                 seed= seed)
      
      summary_list[[i]] <- fit_summary(fit_summary_params)
    }
    
    
    
  }
  
  names(summary_list) <- ids_list
  return(summary_list)
  
}



# SBC Functions ----------------------------------------------------

# 31.03.2025 needs refactoring 


dataset_cond_sbc <- function(sim_cond,
                             sim_params,
                             seed=NULL,
                             path=NULL,
                             ncores=1,
                             exp_name=NULL){
  
  # Set seed for reproducibility.
  # global seed to control the simulations
  if (!is.null(seed)) {
    set.seed(seed)
  }
  seed_list <- sample(1000000000:.Machine$integer.max,
                      size = sim_params$nsims)
  
  if (file.exists(paste0(paste(path, paste0(exp_name,"_",sim_cond$id), sep = "/"), ".RDS"))) {
    return(readRDS(paste0(paste(path, paste0(exp_name,"_",sim_cond$id), sep = "/"), ".RDS")))
  } else {
    if (ncores > 1) {
      # Multiprocessing setup
      
      cluster <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cluster)
      
      parallel::clusterExport(cl= cluster, c('sim_params'))
      parallel::clusterEvalQ(cl = cluster, {
        library(brms)
        library(tidybayes)
        library(cmdstanr)
        library(MASS)
        library(extraDistr)
        library(tibble)
        library(dplyr)
        library(mvtnorm)
        library(LaplacesDemon)
        library(parallel)
        library(Matrix)
        library(gtools)
        library(ppcor)
        library(matrixStats)
        
        
        # scripts
        source("scripts/aux_functions/all_auxfunctions.R") #script containing auxiliary functions
        source("scripts/aux_functions/coef_functions.R") #script containing coefficients auxiliary functions
        source("scripts/aux_functions/cov_functions.R") #script containing covariance auxiliary functions
        source("scripts/aux_functions/dgp_functions.R") #script containing dgp auxiliary functions
        source("scripts/aux_functions/summary_functions.R") #script containing summary auxiliary functions
        source("scripts/aux_functions/mcmc_params.R") #mcmc related functions 
        source("scripts/aux_functions/mcmc_params_data.R") #function to create mcmc params data
        source("scripts/aux_functions/stan_fits.R") #stan fits 
        
        source("scripts/sim/sim_comparison.R") #special function to run the simulation
        source("scripts/aux_functions/R2D2_alpha_gen.R") #how to generate concentration vector from r2d2
        
        
        source("scripts/sbc/sbc_full.R") 
        source("scripts/sbc/sbc_rng.R") 
        
        cmdstanr::set_cmdstan_path(sim_params$cmdstan_path)
        options(matrixStats.vars.formula.freq = 2e6 ) # Weird bug 
        
      })
      
      `%dopar%` <- foreach::`%dopar%`
      
      # Multiprocessing run
      results <- foreach::foreach(
        par_seed = seed_list
      ) %dopar% {
        cond_sbc(
          sim_params = sim_params,
          sim_cond = sim_cond,
          seed = par_seed
        )
        
        
      }
      
      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    }else {
      results <- vector("list", length = length(seed_list))
      par_seed = seed_list
      for (i in seq_along(results)) {
        print(par_seed[i])
        cond_sbc(
          sim_params = sim_params,
          sim_cond = sim_cond,
          seed = par_seed[i]
        )
      }
    } 
    
    final_result <- do.call(rbind, results)
    #final_result$data_config_seed <- seed
    
    if (!is.null(path)) {
      saveRDS(final_result, paste0(paste(path, paste0(exp_name,"_",sim_cond$id), sep = "/"), ".RDS"))
    }
    return(final_result) 
  }
}

#Simulate data and summarise a fit given a simulation condition
#TODO: Update
cond_sbc <-  function(sim_params, sim_cond, seed = NULL){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  sim_cond$seed <- seed 
  
  #-----Cycle through the models 
  
  #------  Fit models
  # Generate parameters of fits
  gen_coef <-   sim_cond$gencoef_function
  fits_params <-  get_sim_mcmc_params_data(sim_params$mcmc_params,
                                           gen_coef, 
                                           sim_cond)
  
  fits_list <- fits_params$fits_list #fits to considered
  
  ids_list <- fits_params$ids_list
  names_list <- fits_params$names_list #names of the mcmcs
  nids <- length(ids_list)  #number of different models
  
  sbc_list <- vector(mode = "list", length = nids)
  params_list <- vector(mode = "list", length = nids)
  sbc_result_list <- vector(mode = "list", length = nids)
  
  for( i in 1:nids){
    #----- Generate prior draw
    # Generate prior draws. Should return a list.
    real_theta <- make_rgn(append( list( name= names_list[i]),
                                   append(fits_params$mcmc_params[[i]],
                                  sim_cond)))
    #------ Generate data
    # Use the prior draw to generate data.
    # dgp should use real value of theta and sim_cond
   
    dgp_list <- make_dgp(list(sim_cond= sim_cond,
                              real_theta = real_theta,
                              sim_params = sim_params))
    
    real_theta <- dgp_list$real_theta #update the real theta vector adding the intercept
    #------ MCMC
    
    mcmc_params <- make_mcmc_params(names_list[i], 
                                    fits_params$mcmc_params[[i]], 
                                    dgp_list)
    
    mcmc_params$iter <-  fits_params$mcmc_params[[i]]$iter
    
    params_list[[i]] <- append(mcmc_params,
                               dgp_list)
    
    fitfn <- get(paste0(fits_list[i],"fit")) # Get fit function
    fit <- fitfn(params_list[[i]]) # Run stan fit
    
    #Summarize
    if(fit$fit.error){
      sbc_result_list[[i]] <- list(  name= names_list[i], 
                                     id = ids_list[i],
                                     seed = seed, 
                                     real_values= unlist(real_theta),
                                smdf= NULL,
                                time=NULL )
    }else{
      
      
      params_names <- c(names(real_theta), 
                  "mu_tilde", 
                  "mu_tilde_test", 
                  "log_lik",
                  "log_lik_test")
      
      pdraws <- as.data.frame(fit$fit$draws( params_names, 
                                         format = "df"))
      
      pdraws <- pdraws[, 1:(dim(pdraws)[2]-3) ] #remove last columns
  
      real_values <- c(unlist(real_theta),
                       dgp_list$mu_tilde,
                       dgp_list$mu_tilde_test,
                       dgp_list$log_lik,
                       dgp_list$log_lik_test)
                       #get real values as a vector
      ranks <- colSums(sweep(pdraws, 2, real_values   , "<")) #calculate ranks
      
      smdf <- fit$fit$summary( params_names, 
                       posterior::default_summary_measures()[1:4],
                       quantiles = ~ posterior::quantile2(., probs = sim_params$smqoi$probsoi),
                       posterior::default_convergence_measures(),
                       prob_gt_0 = ~ mean(. > 0),
                       c("min", "max")  ) %>% 
              mutate( real_theta = real_values, .after= 1 ) %>% 
              mutate( ranks = ranks, .after= 6)  %>%
              mutate( zscore = (mean-real_theta)/sd, .after = 7  )
       
      
      sbc_result_list[[i]] <- list( name= names_list[i], 
                                    id = ids_list[i],
                                    seed = seed, 
                                    real_values= unlist(real_theta),
                                    smdf= smdf,
                                    time=fit$fit$time()$total)
      
    }      
    
  }
  
  return(sbc_result_list)
  
}




# Data Functions ----------------------------------------------------


dataset_dataexp_sim <- function(exp_cond,
                                exp_params,
                                seed=NULL,
                                path=NULL,
                                ncores=1,
                                exp_name=NULL){
  
  # Set seed for reproducibility.
  # global seed to control the simulations
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  if( file.exists(paste0(paste(path, paste0(exp_name,"_",exp_cond$id), sep = "/"), ".RDS")) ) {
    return(readRDS(paste0(paste(path, paste0(exp_name,"_",exp_cond$id), sep = "/"), ".RDS")))
  } else {
    if(ncores > 1) {
      # Multiprocessing setup
      
      cluster <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cluster)
      
      parallel::clusterExport(cl= cluster, c('exp_params'))
      
      parallel::clusterEvalQ(cl = cluster, {
        library(brms)
        library(tidybayes)
        library(cmdstanr)
        
        library(tibble)
        library(dplyr)
        library(mvtnorm)
        library(LaplacesDemon)
        library(parallel)
        library(Matrix)
        library(gtools)
        library(ppcor)
        
        
        source("scripts/aux_functions/all_auxfunctions.R") #script containing auxiliary functions
        source("scripts/aux_functions/full_data_experiment.R") #special function to run the simulation
        
        source("scripts/aux_functions/fit_params.R") #mcmc related functions 
        source("scripts/aux_functions/fit_params_data.R") #function to create mcmc params data
        
        source("scripts/aux_functions/stan_fits.R") #stan fits 
        source("scripts/aux_functions/compile_stan_models.R") #compile and save stan models
        source("scripts/aux_functions/summary_functions.R") # summary functions 
        
        source("scripts/dasp/data/data_mcmc_params_data_lists.R") #function to create mcmc params data
        source("scripts/aux_functions/cov_functions.R") #script containing covariance auxiliary functions
        
        source("scripts/aux_functions/R2_alpha_gen.R") #script containing concentration vector gen functions
        
        
        cmdstanr::set_cmdstan_path(exp_params$cmdstan_path)
        
      })
      
      `%dopar%` <- foreach::`%dopar%`
      
      # Multiprocessing run
      
      if(!exp_cond$loocv_flag){
        # K fold cross validation
        
        seed_list <- sample(1000000000:.Machine$integer.max,
                            size = exp_cond$nreps)
        
        message(paste0("Running K fold cross validation with K= ", exp_cond$nreps) ) 
        
        results <- foreach::foreach(
          par_seed = seed_list
        ) %dopar% {
          dataexp_sim(
            exp_params = exp_params,
            exp_cond = exp_cond,
            seed = par_seed, 
            loo_id = NULL) 
        }
      }else{
        # loo cv
  
        message(paste0("Running LOOCV with n= ", exp_params$list_df$n) ) 
        
        seed_list <- sample(1000000000:.Machine$integer.max,
                            size = exp_params$list_df$n)
        
        loo_ids <- seq_len(exp_params$list_df$n)
        
        results <- foreach::foreach(
          loo_id = loo_ids
        ) %dopar% {
          dataexp_sim(
            exp_params = exp_params,
            exp_cond = exp_cond,
            seed = seed_list[loo_id], 
            loo_id = loo_id)
        }
        
        
      }
  
      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    } else {
      results <- vector("list", length = length(seed_list))
      for (i in seq_along(results)) {
        dataexp_sim(
          exp_params = exp_params,
          exp_cond = exp_cond,
          seed = seed_list[[i]]
        )
      }
    }
   
    final_result <- do.call(rbind, results)
  
    if (!is.null(path)) {
      saveRDS(final_result, paste0(paste(path, paste0(exp_name,"_", exp_cond$id), sep = "/"), ".RDS"))
    }
    return(final_result) 
  }
}

#Simulate data and summarise a fit given a simulation condition
dataexp_sim <-  function(exp_params, exp_cond,  seed = NULL, loo_id = NULL, verbose = TRUE){
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  reps <- exp_cond$reps
  train_size <- exp_cond$train_size
  
  #--- Extract general data info
  
  list_df <- exp_params$list_df
  
  X <- list_df$X
  n <- list_df$n
  p <- list_df$p
  
  df <- list_df$df
  
  #---- Data partition
  
  df$id <- seq_len(nrow(df))
  
  if( !exp_cond$loocv_flag ){
    train <- df %>% dplyr::sample_frac(train_size)
    test <- df[setdiff(df$id, train$id), ]
  }else{
    # loo 
    train <- df[-loo_id, ]
    test <- df[loo_id, ]
  }  
  
  Xtrain <- train[, 1:p]
  Xtest <- test[, 1:p]
  y <- train$y
  ytest <- test$y
  
  cv_id <-  test$id
  
  ntrain <- nrow(Xtrain)
  ntest <- nrow(Xtest)
  
  #--- Collect data params
  
  data_params <- list(n=ntrain,
                      X=Xtrain,
                      ntest= ntest,
                      Xtest= Xtest,
                      p= p, 
                      y= y, 
                      ytest= ytest, 
                      omega_flag  = exp_cond$omega_flag,
                      seed= seed)
  
  data_params$Omega <- calculate_omega(data_params)
  
      
  #------  Fit different models
  # Generate parameters of fits
  fits_params <-  get_sim_fit_params_data(fit_params = exp_params$fit_params,
                                          sim_cond = exp_cond)
  
  fit_options <- list( nchains     = exp_params$nchains, 
                       iter_warmup = exp_params$iter_warmup, 
                       iter_sample = exp_params$iter_sample, 
                       adapt_delta = exp_params$adapt_delta, 
                       prior_only  = 0
                       )
  
  summary_list <- vector(mode = "list", length = length(fits_params$ids_list)  )
  
  for(i in 1: length(fits_params$ids_list) ){
    
    if (verbose) message("------- Fitting model: ", fits_params$ids_list[i], " -------")
    
    tryCatch({
      
      fit_mcmc_params <- make_fit_params(fits_params$names_list[i], 
                                         fits_params$fit_params[[i]], 
                                         data_params)
      
      
      params_list <- c(data_params,
                       fit_mcmc_params, 
                       fit_options)
      
      # Fit the model using the corresponding fitting function
      fitfn <- get(paste0(fits_params$fits_list[i],"fit")) # Get fit function
      fit <- fitfn(params_list) # Run stan fit
      
      if(!fit$fit.error){
      
        fit_summary_params <- list(
          fit= fit$fit,
          standat= params_list,
          voi= exp_params$exp_qoi$voi,
          moi= exp_params$exp_qoi$moi,
          probsoi= exp_params$exp_qoi$probsoi,
          proj_pred_flag = FALSE, 
          seed= seed, 
          exp_cond = exp_cond, 
          cv_id = cv_id)
        
        summary_list[[i]] <- dataexp_summary(fit_summary_params)
        
      }else{
        message("Warning: Fitting error for model ", fits_params$ids_list[i])
      }
    }, error = function(e){
      message("Error occurred for model ", fits_params$ids_list[i], ": ", e$message)
      
    })
  }
  
  names(summary_list) <- fits_params$ids_list
  return(summary_list)
  
  
}










