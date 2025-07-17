#--- Post processing

modified_performance <- function( summary_name, summaries){
  preddf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in 1:nsims){
    perffit <-  as_tibble(summaries[, summary_name ][[i]]$perf)
    perffit <- perffit %>% add_column(simnr=i, .before=1 )  %>% 
      add_column(name= summary_name, .before=1 )
    preddf <- rbind(preddf,perffit)
  }
  
  preddf
}


modified_summary <- function( summary_name, summaries){
  #returns a modified summary of the fit specified by summary_name including all summaries across
  # different simulations
  
  smdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  
  for(i in 1:nsims){
    smfit <-  summaries[, summary_name ][[i]]$sm  
    smfit <- smfit %>% add_column(simnr=i, .before=1 )  
    smfit <- smfit %>% add_column(rtheta= as.numeric(summaries[, summary_name ][[i]]$rtheta) , .after=2)  
    smdf <- rbind(smdf,smfit)
    
  }
  
  smdf <- smdf  %>% 
    mutate( bias= mean-rtheta ) %>% 
    mutate( sqerror= (mean-rtheta)^2 ) 
  
  smdf %>% add_column(summary_name= summary_name, .before=1)
}

calc_coverage <- function(summary, options_list){
  #options list should include upper and lower quantile that forms the CI
  #we then check if the real value of theta is inside the CI or not
  upper= options_list$upper
  lower= options_list$lower
  
  lower <- as.numeric(as.data.frame(summary[,lower][,1])[,1])
  upper <- as.numeric(as.data.frame(summary[,upper][,1])[,1])
  
  summary <- summary %>% add_column(lower= lower, upper= upper)
  
  summary <- summary  %>% 
    mutate( width= upper-lower ) %>% 
    mutate( coverage = rtheta >= lower & rtheta <= upper ) 
  
  summary
  
}


sim_error <- function(summary, options_list){
  
  # H0: theta=0 
  # Ha: theta!=0
  
  # Type I error: reject H0 given it is true
  # Type II error: do not reject H0 given it is false
  # Power: 1- Type II error
  
  # We are "testing" a hypothesis by using CIs. If 0 is in the CI we proceed
  # to reject H0, if it is not then we do not reject H0.
  
  # Type I error: Given the real value is 0, 0 is not in the CI
  # Type II error: Given the real value is not 0, 0 is in the CI and there is also coverage 
  # We know which params are really 0 and which are not. 
  
  # There is a relationship bw coverage and error. We can have coverage and error
  # at the same time. If we cover but we have error, then we give priority to coverage
  # and consider there is no error. 
  # Why? Because in Bayesian Inference we consider the whole
  # posterior and use it to make inference. 
  
  #There is also the possibility that 0 is outside the CI but there is no coverage of the true theta
  
  
  upper= options_list$upper
  lower= options_list$lower
  if(!"coverage" %in% colnames(summary)){
    summary <- calc_coverage(summary, options_list =  options_list )  
  }
  
  
  tempsummary <- as.data.frame(matrix(0, nrow=0, ncol=  0  ))
  Description <- c("Type I bs zero",
                   "Type II bs non-zero")
  
  stats <- list(sum=sum)
  nl <- nlevels(as.factor(summary$summary_name))
  
  #--- Type I error
  
  
  #bs (zero)
  ebsz <- summary %>% mutate_at( vars(variable) , factor)  %>% 
    filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    filter(rtheta==0) %>% # True value is 0
    group_by(summary_name) %>%  #group by summary 
    mutate( error = !(0 >= lower & 0 <= upper)) %>% # 0 is not in the interval and it should be
    mutate( c_error = error & !coverage ) %>%  #Coverage and error: there is coverage but zero is in the interval (should all be F)
    dplyr::select(summary_name, simnr, coverage,error, c_error)  %>% 
    summarise(across(c(coverage,error, c_error), .fns= stats), n=n())  %>% 
    mutate(error= error_sum/n) %>% #error across sims
    mutate( roc = error) %>% 
    mutate(coverage= coverage_sum/n) %>% #coverage across sims
    mutate(cov_error= c_error_sum/n) %>% 
    mutate( roc_cov = cov_error) %>% 
    mutate( type_of_coef =  0)
  
  #--- Type II error
  #bs (nonzero)
  ebsnz <-  summary %>% mutate_at( vars(variable) , factor)  %>% 
    filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    filter(rtheta!=0) %>%  
    group_by(summary_name) %>% 
    mutate( error = 0 >= lower & 0 <= upper ) %>% # 0 is inside the interval
    mutate( c_error = !coverage & error ) %>% #Coverage and error: there is coverage but zero is in the interval
    dplyr::select(summary_name, simnr,, coverage,error, c_error)  %>% 
    summarise(across(c(coverage,error, c_error), .fns= stats), n=n())  %>% 
    mutate(error= error_sum/n) %>% 
    mutate( roc = 1- error) %>% 
    mutate( coverage= coverage_sum/n) %>% 
    mutate( cov_error= c_error_sum/n) %>% 
    mutate( roc_cov = 1-cov_error) %>% 
    mutate( type_of_coef =  1)
  
  
  errorsummary <- bind_rows(ebsz,ebsnz) %>% 
    add_column( "Description"= rep(Description,each=nl), .before=1) %>% 
    add_column( lower=lower, upper=upper) %>%
    dplyr::select(!c(n,error_sum, c_error_sum,  coverage_sum))
  
  return(errorsummary)
  
}


# Summary for a given data config
summary_config <- function(data_gen_conf){
  
  # List of summaries for a given set of configs
  simL= list()
  mcmcL <- list()
  errorL <- list()
  predL <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  metadataL <- data_gen_conf
  
  for(i in 1:dim(data_gen_conf)[1]){
    temp <- readRDS(paste(path_sim2,'/' ,
                          as.character(data_gen_conf$id[i]),'.RDS',sep=''))
    tempL <- list(simdf= bind_rows(temp[,1]), 
                  predm= bind_rows(temp[,2]),
                  metadata= bind_rows(temp[,3])[1,], 
                  upper= data_gen_conf$upper[1],
                  lower= data_gen_conf$lower[1])
    
    tempL$metadata <- c(tempL$metadata,data_gen_conf[i,])
    
    simL[[i]] <- as_tibble(sim_summary(tempL))
    errorL[[i]] = as_tibble(sim_error(tempL))
    predL= rbind( predL,sim_summary_pred(tempL))
    
  }
  
  predL <- predL %>% add_column(HPid= data_gen_conf$HPid, .before=1) 
  
  summary_list <- list(
    sim= simL,
    error= errorL,
    pred= predL,
    metadata= metadataL)
  
  summary_list
  
}
