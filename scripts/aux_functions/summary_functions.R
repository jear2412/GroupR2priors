#Simulate data and summarise a fit given a simulation condition

#---- User summary functions ----

#Percentile estimator
invq <- function(theta,draws){
  #estimator of rqs, such that given theta
  # P( X <= theta)= rqs
  rqs= c()
  for(j in 1:length(theta)){
    rqs[j]= mean(draws[,,j]<theta[j])
  }
  return(rqs) 
}

prmse <-function(standat, fit){
  #Calculate posterior rmse observations
  #posterior rmse
  # Extract the observed data
  y=standat$y
  ytest=standat$ytest
  N=standat$n
  Ntest=standat$ntest
  
  # Extract the posterior draws for y_tilde and mu_tilde
  tildes <-  fit$draws(c("y_tilde","mu_tilde",
                      "y_tilde_test" ,"mu_tilde_test"),
                    format="draws_matrix")
  
  niter <- dim(tildes)[1]
  
  # Split the posterior draws into their respective components
  y_tilde <- tildes[, 1:N]
  mu_tilde <- tildes[, (N+1):(2*N)]
  y_tilde_test <- tildes[, (2*N+1):(2*N+Ntest)]
  mu_tilde_test <- tildes[, (2*N+1+Ntest):(2*N+2*Ntest)]
  
  # Train MSEs
  ytemp=matrix(y, byrow = TRUE, ncol = N, nrow = niter) 
  mse1_train=mean((y-colMeans(y_tilde))^2)
  mse2_train=mean((ytemp-y_tilde)^2)
  mse3_train=mean((y-colMeans(mu_tilde))^2)
  mse4_train=mean((ytemp-mu_tilde)^2)
  
  # Test MSEs
  ytesttemp=matrix(ytest, byrow = TRUE, ncol = Ntest, nrow = niter) 
  mse1_test=mean((ytest-colMeans(y_tilde_test))^2)
  mse2_test=mean((ytesttemp-y_tilde_test)^2)
  mse3_test=mean((ytest-colMeans(mu_tilde_test))^2)
  mse4_test=mean((ytesttemp-mu_tilde_test)^2)
  
  # Combine the results
  mses_train= c(mse1_train,mse2_train,mse3_train,mse4_train)
  mses_test= c(mse1_test,mse2_test, mse3_test, mse4_test)
  
  # Return the square root of the MSEs for both train and test datasets
  return(sqrt(c(mses_train, mses_test)))
  
}  

prmse_theta <- function(fit, theta_name, real_theta, subset_type = "all")  {
  # Calculate posterior RMSE for theta based on the specified subset type.
  #
  # Args:
  #   fit: A fitted Stan model object.
  #   theta_name: The name of the parameter to evaluate.
  #   real_theta: A vector of true values for the parameter.
  #   subset_type: Subset to compute RMSE for. One of "all", "null", or "non-null".
  #
  # Returns:
  #   The RMSE for the specified subset, or NA if no matching subset exists.
  
  theta_hat = as_draws_matrix(fit$draws(variables = theta_name))
  
  # Determine the subset of real_theta to use
  subset_index <- switch(
    subset_type,
    "all" = rep(TRUE, length(real_theta)),
    "null" = real_theta == 0,
    "non-null" = real_theta != 0,
    stop("Invalid subset_type. Must be 'all', 'null', or 'non-null'.")
  )
  
  # If no matching subset, return NA
  if (sum(subset_index) == 0) {
    return(NA)
  }
  
  # Extract relevant subset
  theta_hat_subset <- theta_hat[, subset_index, drop = FALSE]
  real_theta_subset <- real_theta[subset_index]
  
  # Replicate real_theta for RMSE computation
  theta_real_matrix <- matrix(
    real_theta_subset,
    byrow = TRUE,
    nrow = nrow(theta_hat_subset),
    ncol = ncol(theta_hat_subset)
  )
  
  # Compute RMSE
  residuals <- theta_hat_subset - theta_real_matrix
  sqrt(mean( rowSums(residuals^2) ))
  
}

prmse_per_parameter <- function(fit, theta_name, real_theta, subset_type = "all") {
  # Compute per-parameter RMSE and return the mean RMSE across parameters.
  #
  # Args:
  #   fit: A fitted Stan model object.
  #   theta_name: The name of the parameter to evaluate.
  #   real_theta: A vector of true values for the parameter.
  #   subset_type: Subset to compute RMSE for. One of "all", "null", or "non-null".
  #
  # Returns:
  #   The mean of per-parameter RMSEs, or NA if no matching subset exists.
  
  theta_hat <- as_draws_matrix(fit$draws(variables = theta_name))
  
  # Determine the subset of real_theta to use
  subset_index <- switch(
    subset_type,
    "all" = rep(TRUE, length(real_theta)),
    "null" = real_theta == 0,
    "non-null" = real_theta != 0,
    stop("Invalid subset_type. Must be 'all', 'null', or 'non-null'.")
  )
  
  # If no matching subset, return NA
  if (sum(subset_index) == 0) {
    return(NA)
  }
  # Extract relevant subset
  theta_hat_subset <- theta_hat[, subset_index, drop = FALSE]
  real_theta_subset <- real_theta[subset_index]
  
  per_param_rmse <- numeric(ncol(theta_hat_subset))  # Initialize a vector to store RMSEs
  
  for (j in seq_len(ncol(theta_hat_subset))) {
    # Extract the posterior draws for the jth parameter
    draws <- theta_hat_subset[, j]
    
    # Compute the RMSE for the jth parameter
    per_param_rmse[j] <- sqrt(mean((draws - real_theta_subset[j])^2))
  }
  
  # Return mean RMSE across parameters
  mean(per_param_rmse)
  
}



#---- Summaries ####

fit_summary <- function(fit_summary_params){
  # Extract input parameters
  
  fit <-  fit_summary_params$fit  # stan fit
  real_theta <-  fit_summary_params$real_theta #real values
  standat <-  fit_summary_params$standat #stan data used in fit
  voi <-  fit_summary_params$voi #variables of interest
  moi <-  fit_summary_params$moi # metrics of interest
  probsoi <-  fit_summary_params$probsoi #percentiles of interest
  proj_pred_flag <-  fit_summary_params$proj_pred_flag #Run projpred or not?
  seed <-  fit_summary_params$seed #seed
   
  # Define parameter names and real values
  
  params_names <- c(names(real_theta), 
                    "mu_tilde","log_lik",
                    "mu_tilde_test",
                    "log_lik_test")
  
  real_values <- c(unlist(real_theta), 
                   standat$mu_tilde,
                   standat$log_lik, 
                   standat$mu_tilde_test, 
                   standat$log_lik_test)
  
  # Summarize fit results
  model_summary <-  fit$summary(params_names , 
                      moi,
                      quantiles= ~posterior::quantile2(.,probs=probsoi), 
                      prob_gt_0 = ~ mean(. > 0))  %>% 
                      mutate( real_values = real_values, .after = 1) %>% 
                      mutate( zscore = (mean-real_values)/sd) %>% 
                      mutate( time = as.numeric(fit$time()$total) ) %>% 
                      mutate( seed = seed )
                    
  # Compute performance metrics
  
  model_performance <- list( lpd_train =  sum( fit$summary(c("log_lik"), mean)$mean ),
                     lpd_test  = sum(fit$summary(c("log_lik_test"), mean)$mean),
                     rmse_b     = prmse_theta(fit, "beta", real_theta$beta, "all"),
                     rmse_b0    = prmse_theta(fit, "beta", real_theta$beta, "null"),
                     rmse_bn0   = prmse_theta(fit, "beta", real_theta$beta, "non-null"), 
                     rmse_bpp   = prmse_per_parameter( fit , "beta", real_theta$beta, "all"),
                     rmse_bpp0  = prmse_per_parameter( fit , "beta", real_theta$beta, "null"),
                     rmse_bppn0 = prmse_per_parameter( fit , "beta", real_theta$beta, "non-null"),
                     rmse_sigma = prmse_theta(fit, "sigma", real_theta$sigma ),
                     p0         = sum(real_theta$beta==0),
                     pn0 = sum(real_theta$beta!=0),
                     time = as.numeric( fit$time()$total),
                     ndivergent =   sum(fit$diagnostic_summary( quiet = TRUE )$num_divergent),
                     nmaxtreedepth =  sum(fit$diagnostic_summary( quiet = TRUE)$num_max_treedepth), 
                     ebfmi = mean(fit$diagnostic_summary( quiet = TRUE )$ebfmi),
                     seed = seed)

  
  final_result <- list(seed= seed,
                         model_summary = model_summary,
                         model_performance = model_performance)
  
  return(final_result)
  
}


dataexp_summary <- function(fit_summary_params){
  #fit: stan fit
  #qoi: quantities of interest
  #variables of interest:
  
  fit= fit_summary_params$fit
  standat= fit_summary_params$standat
  
  p= standat$p
  
  voi= fit_summary_params$voi
  moi= fit_summary_params$moi 
  probsoi= fit_summary_params$probsoi
  seed= fit_summary_params$seed
  loocv_flag = fit_summary_params$exp_cond$loocv_flag
  cv_id = fit_summary_params$cv_id
  
  # Define parameter names and real values
  
  params_names <- c("mu_tilde",
                    "log_lik",
                    "mu_tilde_test",
                    "log_lik_test")
  
  # Summarize fit results
  model_summary <-  fit$summary(params_names , 
                                moi,
                                quantiles= ~posterior::quantile2(.,probs=probsoi), 
                                prob_gt_0 = ~ mean(. > 0))  %>% 
    mutate( time = as.numeric(fit$time()$total) ) %>% 
    mutate( seed = seed )
  
  # Compute performance metrics
  
  model_performance <- list(lpd_train =  sum(fit$summary(c("log_lik"), mean)$mean),
                            lpd_test  = sum(fit$summary(c("log_lik_test"), mean)$mean),
                            cv_id = cv_id, 
                            time = as.numeric( fit$time()$total),
                            ndivergent =   sum(fit$diagnostic_summary( quiet = TRUE )$num_divergent),
                            nmaxtreedepth =  sum(fit$diagnostic_summary( quiet = TRUE)$num_max_treedepth), 
                            ebfmi = mean(fit$diagnostic_summary( quiet = TRUE )$ebfmi),
                            seed = seed)
  
  # perf <- as.list(c( sum(fit$summary(c("log_lik"), mean)$mean),
  #                    sum(fit$summary(c("log_lik_test"), mean)$mean),
  #                    as.vector(t(fit$loo()$estimates))[1:4],
  #                    as.vector(table(cut(fit$loo()$diagnostics$pareto_k, 
  #                                        breaks=c(-Inf,0.5, 0.7, 1, Inf)))),
  #                    prmse(standat,fit)))
  # 
  # names(perf) <- c("lpd_train", "lpd_test",
  #                  "elpd_loo","elpd_loo_sd","p_loo", "p_loo_sd",
  #                  "paretok_good", "paretok_ok","paretok_bad", "paretok_verybad",
  #                  "rmse1_train", "rmse2_train", "rmse3_train", "rmse4_train",
  #                  "rmse1_test", "rmse2_test", "rmse3_test", "rmse4_test")
  
  
  
  final_result <- list(seed= seed,
                       model_summary = model_summary,
                       model_performance = model_performance)
  

  
  return(final_result)
}



proj_pred_summary <- function(fit_summary_params){
  
  
  fit <-  fit_summary_params$fit  # stan fit
  standat <-  fit_summary_params$standat #stan data used in fit
  
  voi <-  fit_summary_params$voi #variables of interest
  moi <-  fit_summary_params$moi # metrics of interest
  probsoi <-  fit_summary_params$probsoi #percentiles of interest
  
  seed <-  fit_summary_params$seed #seed
  
  #--- Extract data
  X <-  standat$X
  y <-  standat$y
  p <- standat$p
  Xtest <- standat$Xtest
  ytest <-  standat$ytest
  
  fit_projpred <- list(fit = fit, 
                   data = list(X=X, y=y))
  
  draws_sigma <- as.numeric(posterior::as_draws_matrix(fit$draws(c("sigma"))))
  
  true_non_zero_coefs <-  standat$real_theta$beta!=0 
  real_theta <- standat$real_theta
  
  
  #--- Run projpred
  ref <- projpred::init_refmodel(fit_projpred, 
                                 data= data.frame(X=X, y=y), 
                                 formula = y ~.,
                                 family = gaussian(),
                                 ref_predfun= predfun, 
                                 dis= draws_sigma, 
                                 cvfun = NULL,
                                 cvrefbuilder = NA)
  
  nsel <- min(sum(real_theta$beta!=0), 20) #select the number of nonzero coefficients
  
  # Run projpred
  # With current defaults 
  # L1 search could be beneficial since we have no cv
  fit_vs <- varsel(ref, 
                   nterms_max = nsel, 
                   #ndraws = 2000, #ignored
                   method = "L1", #speed
                   seed = seed) #fit var sel
  
  #TODO: run to check if we have the same results
  # how much speed have we gained? If its fast we could run for all the models
  
  # Running project to get posterior draws of projection
  prj_vs <- project(fit_vs, 
                    nterms =  nsel,
                    #refit_prj = FALSE, # we know how many to select already
                    ndraws = 2000, 
                    seed = seed)
  
  # Get index of selected variables
  # We can also use ranking(fit_vs)$fulldata
  prj_index <- as.numeric(stringr::str_replace_all( prj_vs$predictor_terms , "X.", ""))

  prj_pred_coef_index <-  rep(FALSE, p) # which coefficients were selected
  prj_pred_coef_index[sort(prj_index)] <- TRUE

  prj_pred_correct_selection <- prj_pred_coef_index & true_non_zero_coefs #Has Projpred made a correct selection?
  
  # Extract draws
  prj_vs_draws <- as_draws_df(prj_vs)
  prj_vs_draws <- prj_vs_draws[, 1:(ncol(prj_vs_draws)-3) ]
  
  # Rename 
  names(prj_vs_draws) <- c("b_Intercept", paste0(paste0("beta[", prj_index), "]"),"sigma")
  
  # Complete for nonselected variables. Setting everything to zero.
  prj_vs_draws[ paste0(paste0("beta[", c(1:standat$p)[-prj_index]   ), "]") ] <- 0 
  
  # sort in order b_Intercept, beta, sigma
  prj_vs_draws <- prj_vs_draws %>% 
    relocate( c("b_Intercept", paste0(paste0("beta[",  c(1:standat$p) ), "]"),"sigma")) 
  
  real_values <- c(unlist(real_theta)) #order: b_intercept, b, sigma
  
  # Summarise draws from projpred
  smdf_projpred <- posterior::summarise_draws(prj_vs_draws,
                                              moi, 
                                              quantiles= ~posterior::quantile2(.,probs=probsoi), 
                                              prob_gt_0 = ~ mean(. > 0)) %>% 
    mutate( real_values = real_values, .after = 1) %>% 
    mutate( zscore = (mean-real_values)/sd) %>% 
    mutate( time = NULL) %>% 
    mutate( pp_sel = c(TRUE, prj_pred_coef_index, FALSE), .after = 2   ) %>%  # Was it selected? 
    mutate( pp_csel = c(TRUE, prj_pred_correct_selection, FALSE), .after = 2  ) %>%  # Was it selected correctly? 
    mutate( seed = seed )
   
  
  
  #--- calculate log likehood with projected model

  prj_b_Intercept <- prj_vs_draws %>% dplyr::pull('b_Intercept')
  prj_b <- as.matrix(prj_vs_draws[,2:(ncol(prj_vs_draws)-1)])
  prj_sigma <- prj_vs_draws %>% dplyr::pull('sigma')

  prj_log_lik <- t(dnorm(y, 
                         mean =  t(prj_b_Intercept + prj_b %*% t(X)),
                         sd= prj_sigma, log= TRUE))
  
  prj_log_lik_test <- t(dnorm(ytest, 
                         mean =  t(prj_b_Intercept + prj_b %*% t(Xtest)),
                         sd= prj_sigma, log= TRUE))
  
  #--- RMSE
  
  pbeta <- as_draws_matrix(prj_vs_draws[, 2:(p+1)])
  
  real_beta <- real_theta$beta
  zeroes <- real_beta==0
  nonzeroes <- real_beta!=0
  
  real_beta_zeroes <- real_beta[zeroes]
  real_beta_nonzeroes <- real_beta[nonzeroes]
  
  prmse <- sqrt(mean(rowSums( sweep(pbeta, 2, real_beta) * sweep(pbeta, 2, real_beta) ))) 
  prmse_nonzeroes <- sqrt(mean(rowSums( sweep(pbeta[, nonzeroes], 2, real_beta_nonzeroes) * sweep(pbeta[, nonzeroes], 2, real_beta_nonzeroes) ))) 
  prmse_zeroes <- sqrt(mean(rowSums( sweep(pbeta[, zeroes], 2, real_beta_zeroes) * sweep(pbeta[,  zeroes], 2, real_beta_zeroes) ))) 
  
  #prmse_theta_pp <- 
  
  temp <- c()
  for(i in 1:length(real_beta)){
      #rmse per parameter
      temp[i]= sqrt(mean((pbeta[, i]-real_beta[i])^2)) # how far theta_hat_i is from real_theta
  }
  
  prmse_pp <- mean(temp) 
  
  temp <- c()
  pbetatemp <- pbeta[, zeroes]
  for(i in 1:length(real_beta_zeroes)){
    #rmse per parameter
    temp[i]= sqrt(mean((pbetatemp[, i]-real_beta_zeroes[i])^2)) # how far theta_hat_i is from real_theta
  }
  
  prmse_pp0 <- mean(temp) 
  
  
  temp <- c()
  pbetatemp <- pbeta[, zeroes]
  for(i in 1:length(real_beta_zeroes)){
    #rmse per parameter
    temp[i]= sqrt(mean((pbetatemp[, i]-real_beta_zeroes[i])^2)) # how far theta_hat_i is from real_theta
  }
  
  temp <- c()
  pbetatemp <- pbeta[, nonzeroes]
  for(i in 1:length(real_beta_nonzeroes)){
    #rmse per parameter
    temp[i]= sqrt(mean((pbetatemp[, i]-real_beta_nonzeroes[i])^2)) # how far theta_hat_i is from real_theta
  }
  
  prmse_ppn0 <- mean(temp) 
  
  
  # TODO: add loo, discuss with Paul
  perf_projpred <- as.list(c( sum(colMeans(prj_log_lik)),
                     sum(colMeans(prj_log_lik_test)),
                     prmse,
                     prmse_zeroes,
                     prmse_nonzeroes,
                     prmse_pp,
                     prmse_pp0,
                     prmse_ppn0,
                     sum(real_theta$beta==0),
                     sum(real_theta$beta!=0),
                     seed))
  
  names(perf_projpred) <- c("lpd_train", 
                   "lpd_test",
                   "rmse_b", "rmse_b0", "rmse_bn0",
                   "rmse_bpp", "rmse_bpp0", "rmse_bppn0",
                   "p0", "pn0", 
                   "seed")
  
 
  result_projpred <- list( ppsmdf= smdf_projpred , 
        ppperf = perf_projpred)
  
  return(result_projpred)
   
}


#---- Post simulation summary functions ####

# single simulation 
extract_smdf_id <- function(sim_result,exp_name,ids_list, names_list, sim_cond){
  
  smdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in 1:nsims){
    #Extract row corresponding to nth simulation
    foo_row <- sim_result[i,]
    
    #do it for each column
    for(j in 1:length(ids_list)){
        temp <- foo_row[[j]]$sm %>% #extract sm
          add_column(mcmc_name= names_list[j], .before=1) %>% 
          add_column(id_model= ids_list[j], .before=1) %>% 
          add_column(sim_cond, .before= 1) %>% #add simulation conditions
          add_column(simnr=i, .before=1 ) 
        
        smdf <- rbind(smdf, temp)
    }
  }
  
  smdf 
}

# single simulation 
extract_ppsmdf_id <- function(sim_result,exp_name,ids_list, names_list, sim_cond){
  
  ppsmdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in 1:nsims){
    #Extract row corresponding to nth simulation
    foo_row <- sim_result[i,]
    
    #do it for each column
    for(j in 1:length(ids_list)){
      temp <- foo_row[[j]]$ppsmdf %>% #extract sm projpred
        add_column(mcmc_name= names_list[j], .before=1) %>% 
        add_column(id_model= ids_list[j], .before=1) %>% 
        add_column(sim_cond, .before= 1) %>% #add simulation conditions
        add_column(simnr=i, .before=1 ) 
      
      ppsmdf <- rbind(ppsmdf, temp)
    }
  }
  
  ppsmdf 
}


extract_smdf_sbc_id <- function(sbc_result,ids_list, names_list, sim_cond){
  
  smdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in 1:nsims){
    #Extract row corresponding to nth simulation
    foo_row <- sbc_result[i,]
    
    #do it for each column
    for(j in 1:length(ids_list)){
        
        temp <- foo_row[[j]]$sm %>%
          add_column(id_model= ids_list[j], .before=1) %>% 
          add_column(mcmc_name= names_list[j], .after=1) %>% 
          add_column(sim_cond, .before= 1) %>% #add simulation conditions
          add_column(simnr=i, .before=1 ) %>% 
          unite( "variable_sbc", all_of(c('id_model','variable')),  sep= "_" , remove= FALSE) 
        
        smdf <- rbind(smdf, temp)
    }
    
  }
  
  smdf 
}

extract_preddf_id <- function(sim_result,ids_list, names_list, sim_cond){
  
  preddf <- data.frame() 
  for(i in 1:nsims){
    #Extract row corresponding to nth simulation
    foo_row <- sim_result[i,]
    
    #do it for each column
    for(j in 1:length(ids_list)){
      
      temp <- as_tibble(foo_row[[j]]$model_performance) %>%
        add_column(mcmc_name= names_list[j], .before=1) %>% 
        add_column(id_model= ids_list[j], .before=1) %>% 
        add_column(sim_cond, .before= 1) %>% #add simulation conditions
        add_column(simnr=i, .before=1 ) 
      
      
      preddf <- rbind(preddf, temp)
    }
    
    
  }
  
  preddf 
}


process_pred_result_id  <- function(sim_result, ids_list, names_list, sim_cond){
 
  # Extract the first row from sim_result (assuming sim_result is a matrix or data frame)
  foo_row <- sim_result[1,]
  
  # Use map2 to iterate over ids_list and names_list
  results_list <- map2(ids_list, names_list, function(id, name) {
    # Check if the id exists in foo_row and is not NULL
    if (id %in% names(foo_row) && !is.null(foo_row[[id]])) {
      # Extract 'perf' data frame from the list
      perf_df <- foo_row[[id]]$perf
      
      # Create a tibble with additional columns
      as_tibble(perf_df) %>%
        add_column(mcmc_name = name, .before = 1) %>%
        add_column(id_model = id, .before = 1) %>%
        add_column(sim_cond, .before = 1) %>%
        add_column(simnr = 1, .before = 1)  # Simnr is 1 since nsims is 1
    } else {
      message("Warning: Data missing for id_model ", id)
      NULL  # Return NULL if data is missing
    }
  })
  
  # Combine all results into one data frame, filtering out NULLs
  preddf <- bind_rows(results_list[!sapply(results_list, is.null)])
  
  return(preddf)

}

