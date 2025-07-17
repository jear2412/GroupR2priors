#---- Functions  ----

# Main function to process all directories
process_all_directories <- function(type, result_type = "performance" ) {
  # List all directories that match the pattern
  dirs <- list.dirs("results/group_r2", full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl(paste0("arrayid_[0-9]+_", type), dirs)]
  
  # Process each directory and bind the results into one big dataframe
  
  if( result_type == "performance"){
    final_data <- map_dfr(dirs, process_directory_performance)
    #final_data <- future_map_dfr(dirs, process_directory_performance)  
  }
  
  if( result_type == "summary"){
    final_data <- map_dfr(dirs, process_directory_summary)    
  }
  
  
  
  return(final_data)
}



# Main function to process all directories in parallel
process_all_directories_parallel <- function(path_results , result_type = "performance") {
  # List all directories that match the pattern
  dirs <- list.dirs(path_results, full.names = TRUE, recursive = FALSE)
  #dirs <- dirs[grepl(paste0("arrayid_[0-9]+_", type), dirs)]
  
  # Set up parallel backend
  plan(multisession)  # Choose parallel backend (e.g., multisession, multicore, etc.)
  
  # Process each directory in parallel and bind the results into one big dataframe
  if( result_type == "performance"){
    final_data <- future_map_dfr(dirs, process_directory_performance)  
  }
  
  if( result_type == "summary"){
    final_data <- future_map_dfr(dirs, process_directory_summary)    
  }
  
  return(final_data)
}


process_directory_performance <- function(dir) {
  # List all .RDS files in the directory
  
  files <- list.files(dir, pattern = "^array.*\\.RDS$", full.names = TRUE)
  
  # Apply the function to each file and combine the results
  files_data <- map_dfr(files, ~ process_file_performance(.x))
  
  return(files_data)
}

# Function to process a single directory
process_directory_summary <- function(dir) {
  # List all .RDS files in the directory
  
  files <- list.files(dir, pattern = "^array.*\\.RDS$", full.names = TRUE)
  
  # Apply the function to each file and combine the results
  files_data <- map_dfr(files, ~ process_file_summary(.x))
  
  return(files_data)
}


process_file_summary <- function(filename) {
  tryCatch({
    # Extract simulation number and ID number from the special format we have
    sim_number <- regmatches(filename, regexpr("(?<=arrayid_)([0-9]{1,2}|100)(?=_)", filename, perl = TRUE))
    id_number <- regmatches(filename, regexpr("(?<=_)([1-9][0-9]?|100)(?=\\.RDS$)", filename, perl = TRUE))
    
    # Construct full path to RDS file
    filepath <- filename
    
    # Extract and process prediction data
    foopreddf <- process_summary_result_id(readRDS(filepath), 
                                           ids_list, 
                                           names_list, 
                                           sim_conds[sim_conds$id == id_number,])
    
    # Add simulation number and return
    foopreddf %>% 
      mutate(simnr = sim_number)
    
  }, error = function(e) {
    message("Error processing file ", filename, ": ", e$message)
    NULL  # Return NULL on error to continue processing other files
  })
}



process_file_performance <- function(filename) {
  tryCatch({
    # Extract simulation number and ID number from the special format we have
    sim_number <- regmatches(filename, regexpr("(?<=arrayid_)([0-9]{1,2}|100)(?=_)", filename, perl = TRUE))
    id_number <- regmatches(filename, regexpr("(?<=_)([1-9][0-9]?|100)(?=\\.RDS$)", filename, perl = TRUE))
    
    # Construct full path to RDS file
    filepath <-  filename
    
    # Extract and process prediction data
    foopreddf <- process_performance_result_id(readRDS(filename), 
                                               ids_list, 
                                               names_list, 
                                               sim_conds[sim_conds$id == id_number,])
    
    # Add simulation number and return
    foopreddf %>% 
      mutate(simnr = sim_number)
    
  }, error = function(e) {
    message("Error processing file ", filename, ": ", e$message)
    NULL  # Return NULL on error to continue processing other files
  })
}


process_summary_result_id  <- function(sim_result, ids_list, names_list, sim_cond){
  
  # Extract the first row from sim_result (assuming sim_result is a matrix or data frame)
  foo_row <- sim_result[1,]
  
  # Use map2 to iterate over ids_list and names_list
  results_list <- map2(ids_list, names_list, function(id, name) {
    # Check if the id exists in foo_row and is not NULL
    if (id %in% names(foo_row) && !is.null(foo_row[[id]])) {
      # Extract 'perf' data frame from the list
      summary_df <- foo_row[[id]]$model_summary
      
      # Create a tibble with additional columns
      as_tibble(summary_df) %>%
        add_column(mcmc_name = name, .before = 1) %>%
        add_column(id_model = id, .before = 1) %>%
        add_column(sim_cond, .before = 1) %>%
        add_column(simnr = 1, .before = 1)  # Simnr is 1 since nsims is 1
        #add_column( nchain = 1:nchains )
    } else {
      message("Warning: Data missing for id_model ", id)
      NULL  # Return NULL if data is missing
    }
  })
  
  # Combine all results into one data frame, filtering out NULLs
  summary_df <- bind_rows(results_list[!sapply(results_list, is.null)])
  
  return(summary_df)
  
}

process_performance_result_id  <- function(sim_result, ids_list, names_list, sim_cond){
  
  # Extract the first row from sim_result (assuming sim_result is a matrix or data frame)
  foo_row <- sim_result[1,]
  
  # Use map2 to iterate over ids_list and names_list
  results_list <- map2(ids_list, names_list, function(id, name) {
    # Check if the id exists in foo_row and is not NULL
    if (id %in% names(foo_row) && !is.null(foo_row[[id]])) {
      # Extract 'perf' data frame from the list
      perf_df <- foo_row[[id]]$model_performance
      
      # Create a tibble with additional columns
      as_tibble(perf_df) %>%
        add_column(mcmc_name = name, .before = 1) %>%
        add_column(id_model = id, .before = 1) %>%
        add_column(sim_cond, .before = 1) %>%
        add_column(simnr = 1, .before = 1)   # Simnr is 1 since nsims is 1
        #add_column( nchain = 1:nchains ) 
    } else {
      message("Warning: Data missing for id_model ", id)
      NULL  # Return NULL if data is missing
    }
  })
  
  # Combine all results into one data frame, filtering out NULLs
  preddf <- bind_rows(results_list[!sapply(results_list, is.null)])
  
  return(preddf)
  
}




create_performance_df <- function( sim_conds,  ppdf,  quantities){
  performance_df <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  
  for(i in 1:length(quantities)){
    qoi <- quantities[i]
    
    temp <- ppdf %>% 
      ungroup() %>% 
      dplyr::select(!c(mcmc_name)) %>% 
      group_by(across(all_of( grouping_cols  ))) %>% 
      mutate(rn = row_number(), .before = 1) %>% 
      dplyr::select(  qoi )  %>% 
      pivot_wider(names_from = id_model, 
                  values_from = qoi) %>% 
      ungroup()
    
    # Caution: clumsy code that is hard coded!
    # columns are subtracted even if they have na, check this someday
    
    temp <- temp %>% 
      pivot_longer(starts_with("dif" , ignore.case = FALSE) ,
                   names_to = "quantity",
                   values_to = "value") %>%  
      mutate( var = qoi) %>% 
      mutate( quantity = as.factor(quantity))
    
    
    difperformance_df <- rbind(difperformance_df, temp)
    
  }
  
  difperformance_df 
  
}


create_preddf <- function( sim_conds,  ppdf,  quantities){
  preddf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  grouping_cols <- c(names(sim_conds), "simnr",  "id_model" )
  
  for(i in 1:length(quantities)){
    qoi <- quantities[i]
    
    temp <- ppdf %>% 
      ungroup() %>% 
      dplyr::select(!c(mcmc_name)) %>% 
      group_by(across(all_of( grouping_cols  ))) %>% 
      mutate(rn = row_number(), .before = 1) %>% 
      dplyr::select(  qoi )  %>% 
      pivot_wider(names_from = id_model, 
                  values_from = qoi) %>% 
      ungroup()
    
    # Caution: clumsy code that is hard coded!
    # columns are subtracted even if they have na, check this
    
    temp <- temp %>% 
      pivot_longer(starts_with("dif" , ignore.case = FALSE) ,
                   names_to = "quantity",
                   values_to = "value") %>%  
      mutate( var = qoi) %>% 
      mutate( quantity = as.factor(quantity))
    
    
    difpreddf <- rbind(difpreddf, temp)
    
  }
  
  difpreddf 
  
}

calculate_differences <- function(df, diffs) {
  # Helper function to handle differences dynamically
  for (name in names(diffs)) {
    cols <- diffs[[name]]
    if (all(cols %in% names(df))) { 
      df <- df %>%
        mutate( !!sym(name) := !!sym(cols[1]) - !!sym(cols[2]))
    }
  }
  return(df)
}

create_difpreddf <- function(sim_conds, ppdf, quantities) {
  
  # Define the column pairs for differences in a dynamic list
  # HARD coded
  # New model, Old model
  diffs <- list(
    difD2d1 = c("D2def1g", "D2default"), 
    difD2d2 = c("D2def2g", "D2default"), 
    difD2d3 = c("D2def3g", "D2default"),
    difD2u  = c("D2unifg", "D2unif"),
    difD2n1 = c("D2new1g", "D2new"),
    difD2n2 = c("D2new2g", "D2new"),
    difD2n3 = c("D2new3g", "D2new3"),
    difD2n4 = c("D2new4g", "D2new4"),
    difD2n5 = c("D2new5g", "D2new5"),
    difD2n6 = c("D2new6g", "D2new6"),
    difGIG = c("gigg", "BP"),
    difGIGnn = c("giggnn", "BP"),
    difGIGh = c("gigghier", "BPhier"),
    difGIGnh = c("giggnh", "BPhier")
  )
  
  grouping_cols <- c(names(sim_conds), "simnr", "id_model")

  # Calculate differences
  # New model - old model (according to diffs list)
  difpreddf <- map_dfr(quantities, function(qoi) {
    
    reshaped_data <- ppdf %>%
      dplyr::select(-mcmc_name) %>%
      mutate(rn = row_number(), .before = 1) %>%
      dplyr::select(qoi, id_model, all_of(grouping_cols)) %>%
      pivot_wider(names_from = id_model, values_from = qoi)

    
    # Apply the dynamic difference calculations using the helper function
    reshaped_data <- calculate_differences(reshaped_data, diffs)
    
    # Reshape the data into long format
    reshaped_data %>%
      pivot_longer(starts_with("dif", ignore.case = FALSE), 
                   names_to = "quantity", 
                   values_to = "value") %>%
      mutate(var = qoi, 
             quantity = as.factor(quantity))
  })
  
  return(difpreddf)
}

qoi_plot <- function(dat, qoi, form, 
                     xlabel = "",
                     ylabel = "", title = "", 
                     plot_prefix = "",
                     save_flag = TRUE, 
                     width = 17, 
                     height= 10,
                     delta_flag = FALSE,
                     bound = NULL ){
  
  # Apply the filter only if `bound` is not NULL
  if (!is.null(bound)) {
    dat <- dat %>% 
      filter( quantity == qoi) %>% 
      filter(abs(value) < bound)
  }
  
  if(delta_flag){
    plot_qoi <- dat %>% 
      filter( quantity == qoi) %>% 
      ggplot(aes(x= id_model , 
                 y= value,
                 fill= id_model))+
      geom_violin() +
      geom_boxplot(width = 0.1, outlier.size = 0.3) +
      stat_summary(fun=median, 
                   geom = "point",
                   color = "black", 
                   size = 1)+
      scale_fill_paletteer_d(my_palette)+
      ggh4x::facet_nested( form  ,
                           labeller = labeller(R2 = label_R2, 
                                               p = label_p, 
                                               n = label_n,
                                               rho = label_rho, 
                                               type = label_gcf
                           ),
                           independent = "y",
                           scales = "free"
      ) +
      #geom_hline( yintercept = 0,  
      #            linetype = 2, 
      #            color= my_hline_color ) +
      labs(
        y = bquote(Delta ~ .(ylabel)),
        x = xlabel
      )+
      ggtitle(title)+
      paper_theme +
      theme(legend.position="none")
  }else{
    
    plot_qoi <- dat %>% 
      filter( quantity == qoi) %>% 
      ggplot(aes(x= id_model , 
                 y= value,
                 fill= id_model))+
      geom_violin() +
      geom_boxplot(width = 0.2, outlier.size = 0.3) +
      stat_summary(fun=median, geom = "point", color = "black", size = 1)+
      scale_fill_paletteer_d(my_palette)+
      ggh4x::facet_nested( form  ,
                           labeller = labeller(R2 = label_R2, 
                                               p = label_p, 
                                               n = label_n,
                                               rho = label_rho, 
                                               type = label_gcf
                           ),
                           independent = "y",
                           scales = "free"
      ) +
      labs(
        y = ylabel,
        x = xlabel
      )+
      ggtitle(title)+
      paper_theme +
      theme(legend.position="none")
  }
  
  
  if(save_flag){
    # Construct filename
    filename <- file.path(results_folder,
                          paste0( plot_prefix, "_",
                                  qoi, "-", 
                                  experiment_name, ".pdf"))
    
    # Save the plot
    save_plot(plot_qoi, filename, width = width, height = height)
  }
  
  plot_qoi
  
}



# Plot difference between grouped and nongrouped versions
diff_qoi_plot <- function(dat, qoi, form, 
                          xlabel = "",
                          ylabel = "", 
                          title = "", 
                          plot_prefix = "",
                          save_flag = TRUE, 
                          width = 17, 
                          height= 10.5,
                          delta_flag = FALSE,
                          bound = NULL ){
  
  # Apply the filter only if `bound` is not NULL
  if (!is.null(bound)) {
    dat <- dat %>% 
      filter( var == qoi) %>% 
      filter(abs(value) < bound)
  }
  
  
  plot_qoi <-  dat %>% 
    filter( var == qoi) %>% 
    ggplot(aes(x = quantity , 
               y = value,
               fill = quantity))+
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    stat_summary(fun=median, 
                 geom = "point", 
                 color = "black",
                 size = 1)+
    scale_fill_paletteer_d(my_palette)+
    ggh4x::facet_nested( form  ,
                         labeller = labeller(R2 = label_R2, 
                                             p = label_p, 
                                             rho = label_rho, 
                                             n = label_n,
                                             type = label_gcf ),
                         independent = "y", 
                         scales = "free") +
    geom_hline( yintercept = 0,  
                linetype = 2, 
                color= my_hline_color ) +
    paper_theme +
    labs(
      y = bquote(Delta ~ .(ylabel)),
      x = xlabel, 
    )+
    ggtitle(title)+
    #theme(strip.background=element_rect(color="grey30", 
    #                                    fill="grey90"))+
    theme(legend.position="none", 
          axis.text.x = element_text(angle= 45, hjust = 0.5)
          #plot.margin = margin(t = 10, r = 10, b = 30, l = 10) 
    )
  
  
  
  if(save_flag){
    # Construct filename
    filename <- file.path(results_folder,
                          paste0( plot_prefix, "_delta_",
                                  qoi, "-", 
                                  experiment_name, ".pdf"))
    
    # Save the plot
    save_plot(plot_qoi, filename, width = width, height = height)
  }
  
  plot_qoi
  
}


save_plot <- function(plot, filename, width = 17, height = 11) {
  # Helper function to save a plot
  ggsave(filename = filename, plot = plot, width = width, height = height)
}

process_and_save_plots <- function(dat, quantities, form, results_folder, 
                                   plot_prefix, plot_func, 
                                   delta_flag = FALSE) {
  
  # Function to create and save plots for each quantity
  
  # Ensure results_folder exists
  if (!dir.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
  }
  
  # Iterate over each quantity and create/save the plot
  
  if( !delta_flag){
    walk(quantities, function(quantity) {
      
      # Create the plot
      plot <- plot_func(dat, quantity, form) +
        ggtitle(quantity) +
        labs(
          y = bquote(Delta ~ .(quantity)),
          x = "Model"
        )
      
      # Construct filename
      filename <- file.path(results_folder,
                            paste0( plot_prefix,
                                    quantity, ".pdf"))
      
      # Save the plot
      save_plot(plot, filename)
    })
  } else {
    walk(quantities, function(quantity) {
      
      # Create the plot
      plot <- plot_func(dat, quantity, form) +
        ggtitle(quantity) +
        labs(
          y = quantity,
          x = "Model"
        )
      
      # Construct filename
      filename <- file.path(results_folder,
                            paste0( plot_prefix,
                                    quantity, ".pdf"))
      
      # Save the plot
      save_plot(plot, filename)
    })
  }
  
  
}

delta_metrics_long_format <- function(data, id_model_vector) {
  # Filter dataset based on the provided id_model vector
  
  grouping_cols <- c(names(sim_conds), "simnr")
  
  processed_data <- data %>%
    mutate(id_model = case_when(
      id_model == "D2dd" ~ "D2",
      id_model == "D2du" ~ "D2du",
      id_model == "D2ud" ~ "D2ud",
      id_model == "D2uu" ~ "D2uu",
      id_model == "L2dd" ~ "LNF",
      id_model == "L2dd_s" ~ "LNS",
      id_model == "L2du" ~ "LNdu",
      id_model == "L2du_s" ~ "LNdus",
      id_model == "L2ud" ~ "LNud",
      id_model == "L2ud_s" ~ "Luds",
      id_model == "L2uu" ~ "LNuu",
      id_model == "L2uu_s" ~ "LNuus",
      TRUE ~ id_model  # Keep original value if no match
    )) %>% 
    filter(id_model %in% id_model_vector) %>%
    mutate(id_model = factor(id_model, levels = id_model_vector)) %>%
    # Convert simulation condition columns to factors
    mutate(across(all_of(grouping_cols), as.factor)) %>%
    # Group by the specified columns
    group_by(across(all_of(grouping_cols))) %>%
    # Calculate delta quantities
    mutate(delta_rmse_b = rmse_b - min(rmse_b, na.rm = TRUE),
           delta_rmse_b0 = rmse_b0 - min(rmse_b0, na.rm = TRUE),
           delta_rmse_bn0 = rmse_bn0 - min(rmse_bn0, na.rm = TRUE),
           delta_rmse_bpp = rmse_bpp - min(rmse_bpp, na.rm = TRUE),
           delta_rmse_bpp0 = rmse_bpp0 - min(rmse_bpp0, na.rm = TRUE),
           delta_rmse_bppn0 = rmse_bppn0 - min(rmse_bppn0, na.rm = TRUE),
           delta_lpd_train = lpd_train - max(lpd_train, na.rm = TRUE),
           delta_lpd_test = lpd_test - max(lpd_test, na.rm = TRUE)) %>%
    ungroup() %>%
    select(-c(mcmc_name)) %>% 
    # Pivot to long format
    pivot_longer(
      cols = -all_of( c( grouping_cols, "id_model" )), # Keep grouping columns unchanged
      names_to = "var",           # Column for original column names
      values_to = "value"            # Column for their values
    )
  
  return(processed_data)
}



# Define a function to compute TPR and FPR at different thresholds
compute_roc_data <- function(data, lower_q, upper_q) {
  data %>%
    group_by(p, n,  R2, id_model, type) %>%  
    mutate(predicted_positive = as.numeric(get(lower_q) > 0 | get(upper_q) < 0),  # CI excludes 0
           actual_positive = as.numeric(real_values != 0)) %>%  # True effect present
    summarise(
      TP = sum(predicted_positive == 1 & actual_positive == 1),
      FP = sum(predicted_positive == 1 & actual_positive == 0),
      TN = sum(predicted_positive == 0 & actual_positive == 0),
      FN = sum(predicted_positive == 0 & actual_positive == 1)
    ) %>%
    mutate(TPR = TP / (TP + FN),  # Sensitivity
           FPR = FP / (FP + TN))  # 1 - Specificity
}

