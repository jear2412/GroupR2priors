#---- Setup ----

library(tidyverse)
library(tidybayes)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(corrplot)
library(ggcorrplot)
library(paletteer)
library(patchwork)
library(ggh4x)
library(purrr)
library(future)
library(furrr)
library(viridis)
library(patchwork)


# scripts
source("scripts/aux_functions/all_auxfunctions.R") #script containing auxiliary functions
source("scripts/aux_functions/post_process_functions.R") # functions used in post processing data
source("scripts/aux_functions/summary_functions.R") #script containing summary auxiliary functions
source("scripts/sim_group/plot_utils.R") #script containing summary auxiliary functions

#--- Read conditions ----

# Indexed by n and p
# Options
# (n, p) = (500, 50),  
# (n, p) = (50, 200) 
# (n, p) = (200, 100)
# (n, p) = (200, 500)

#--- Fixed conditions df
# The following assumes that all conditionsdf were fit with the same set models

# Fill the conditionsdf to proceed 
conditionsdf <- data.frame( n = c(200), 
                            p = c(100), 
                            rho_out = 0.2, 
                            spars = 0.75)

directory_names <- with(conditionsdf, paste0("n_", n, "_p_", p))

nchains <- 4 # how many chains were used

experiment_name <- "exp1"

# Initialize lists
sim_conds_list <- list()
ppdf_list <- list()
summarydf_list <- list()

# Flags
performance_flag <- TRUE 
only_D2_flag <- TRUE #should we only plot D2 models?

for( i in seq_along(directory_names)){
  
  directory_name <- directory_names[i]
  
  path_results <- file.path("results/group_r2", directory_name)  # Hardcoded DO NOT CHANGE
  
  # Extract simulation conditions
  sim_setup <- readRDS(file.path(path_results, paste0(directory_name, "_sim_setup.rds")))
  nsims <- sim_setup$nsims
  
  sim_conds_list[[i]] <- read.table(file= paste0(path_results,"/",
                                       paste0(directory_name,
                                              "_sim_conds.txt")), header= TRUE)
  
  nconds <- nrow(sim_conds_list[[i]])
  
  sim_conds <- sim_conds_list[[i]]
  fit_params_list <- readRDS(file.path(path_results, paste0(directory_name, "_mcmc_lists.rds")))
  
  ids_list <- sapply(fit_params_list, function(x) x$id, USE.NAMES = FALSE)
  names_list <- sapply(fit_params_list, function(x) x$mcmc_name, USE.NAMES = FALSE)
  
  # Create performance dataframe
  
  if(performance_flag){
    ppdf_list[[i]] <- process_all_directories_parallel( path_results = path_results, 
                                              result_type = "performance" ) %>%
                      mutate(nsims = nsims, nconds = nconds)
  }
    
}

ppdf <- do.call(rbind, ppdf_list )
sim_conds <- do.call(rbind, sim_conds_list) %>% 
              mutate(id = row_number())
# Process results 

if(performance_flag){
  
  
  ppdf <- ppdf %>%
    mutate(across(all_of( names(sim_conds) ), as.factor)) %>% 
    group_by(across(all_of( c( names(sim_conds), "simnr", "id")))) %>% 
    mutate(delta_rmse_b = rmse_b - min(rmse_b, na.rm = TRUE),
           delta_rmse_b0 = rmse_b0 - min(rmse_b0, na.rm = TRUE),
           delta_rmse_bn0 = rmse_bn0 - min(rmse_bn0, na.rm = TRUE),
           delta_rmse_bpp = rmse_bpp - min(rmse_bpp, na.rm = TRUE),
           delta_rmse_bpp0 = rmse_bpp0 - min(rmse_bpp0, na.rm = TRUE),
           delta_rmse_bppn0 = rmse_bppn0 - min(rmse_bppn0, na.rm = TRUE),
           delta_lpd_train= lpd_train - max(lpd_train, na.rm = TRUE),
           delta_lpd_test= lpd_test - max(lpd_test, na.rm = TRUE)) %>% 
    ungroup()
  
  
  grouping_cols <- c(names(sim_conds), 
                     "simnr", "id",  "id_model", "mcmc_name")
  
  preddf <- ppdf %>% 
    pivot_longer(!all_of(grouping_cols) , 
                 names_to = "quantity",
                 values_to = "value") %>% 
    mutate( quantity = as.factor(quantity))
  
  # Create preddf differences 
  quantities= c("lpd_train", 
                "lpd_test", 
                "rmse_b", 
                "rmse_b0", 
                "rmse_bn0",
                "rmse_bpp", 
                "rmse_bpp0",
                "rmse_bppn0")
  
  # fix
  difpreddf <- create_difpreddf(sim_conds, ppdf, quantities)
  
}

#---- ggplot2 settings ####

my_palette <- "LaCroixColoR::paired"

paper_theme <- theme_bw()+ 
  theme(plot.title = element_text(size=16),
        #axis.line.y = element_blank(),
        #panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 17),
        strip.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(angle= 90, vjust= 0.5, hjust= 1, size = 16),
        axis.text.y = element_text(size = 16))

my_hline_color <- "#9a031e" #nice color for vertical or horizontal lines

results_folder <- file.path("results/group_r2", "plots") 

if(!dir.exists(results_folder)){
  dir.create(file.path(results_folder))
}

# Recode your df 
# Named vectors for recoding

label_n <- c("500" = "n = 500", 
             "200" = "n = 200")

label_p <- c("50" = "p = 50",
             "100" = "p = 100",
             "200" = "p = 200", 
             "500" = "p = 500")

label_rho <- c("0.5" = "rho = 0.5", 
               "0.9" = "rho = 0.9")

label_R2 <- c("0.25" = "R2 = 0.25", 
              "0.8" = "R2 = 0.8")

label_gcf <- c("dist" = "Dist", 
               "con" = "Con", 
               "rdist" = "Rdist",
               "rcon" = "Rcon",
               "rcoef" = "Rcoef" )

#---- Filter data ----


#---- Plots ----

form_left <- "type"
form_right <- "R2"
form <- as.formula( paste(form_left, "~", form_right))

bound = 10000

levels(as.factor(preddf$id_model))


# Recode model names
plotdf <- preddf %>% 
  mutate(id_model = recode_factor(id_model, 
                                  "D2def1g" = "D2gd1", 
                                  "D2def2g" = "D2gd2", 
                                  "D2default" = "D2d", 
                                  "D2new" = "D2n",
                                  "D2new1g"  = "D2gn1",
                                  "D2new2g"  = "D2gn2", 
                                  "D2unif" = "D2u", 
                                  "D2unifg" = "D2gu", 
                                  "gigg"= "GIGG", 
                                  "gigghier"= "GIGGh", 
                                  "giggnh" = "GIGGnh", 
                                  "giggnn" = "GIGGnn"
                                  )) 

difpreddf$quantity <- as.factor(difpreddf$quantity)

difplotdf <- difpreddf 

if( only_D2_flag  ){
  
  plotdf <- plotdf %>% 
    filter(str_starts(id_model, "D2"))
  
  difplotdf <- difpreddf  %>% 
    filter(str_starts(quantity, "difD2"))
  
}



difplotdf <- difplotdf %>% 
  filter( !quantity %in% c("difD2d1", 
                           "difD2d2", 
                           #"difD2n1", 
                           "difGIGnn",
                           "difGIGnh",
                           "difD2n6"
                           )) %>% 
  mutate(quantity = recode_factor(quantity, 
                                  "difD2d1" = "def1", 
                                  "difD2d2" = "def2", 
                                  "difD2n1" = "R2-c", 
                                  "difD2n2" = "R2-d",
                                  "difD2n3" = "R2-0.1",
                                  "difD2n4" = "R2-0.5",
                                  "difD2n5" = "R2-1.0",
                                  "difD2n6" = "R2-5.0",
                                  "difD2u"  = "R2-u",
                                  "difGIG"  = "GIGG", 
                                  "difGIGh" = "GIGGh", 
                                  "difGIGnn" = "GIGGnn", 
                                  "difGIGnh" = "GIGGnh"
                                  ))



plotdf <- preddf %>% 
  mutate(id_model = recode_factor(id_model, 
                                  "D2def1g" = "D2gd1", 
                                  "D2def2g" = "D2gd2", 
                                  "D2default" = "R2-w", 
                                  "D2new" = "R2-c",
                                  "D2new2" = "R2-d",
                                  "D2new3" = "R2-0.1",
                                  "D2new4" = "R2-0.5",
                                  "D2new5" = "R2-1.0",
                                  "D2new6" = "R2-5.0",
                                  "D2newg" = "GR2-c",
                                  "D2new1g" = "GR2-c",
                                  "D2new2g" = "GR2-d",
                                  "D2new3g" = "GR2-0.1",
                                  "D2new4g" = "GR2-0.5",
                                  "D2new5g" = "GR2-1.0",
                                  "D2new6g" = "GR2-5.0",
                                  "D2unif" = "R2-u", 
                                  "D2unifg" = "GR2-u", 
                                  "gigg"= "GIGG", 
                                  "gigghier"= "GIGGh", 
                                  "giggnh" = "GIGGnh", 
                                  "giggnn" = "GIGGnn",
                                  "BPhier" = "BPh"
  )) %>% 
  filter( !id_model %in% c( "R2-w", "R2-5.0", "GR2-5.0"
  )) 



asinh_df <- difplotdf %>%
  filter(var == "lpd_test") %>%
  mutate(
    value = asinh(value),
    var = "lpd_test_asinh"
  )

# Bind to original
difplotdf <- bind_rows(difplotdf, asinh_df)

# lpd ----

# lpd train

qoi_plot(dat = plotdf, 
         xlabel = "Prior",
         ylabel="LPD train" , 
         qoi = "lpd_train",  
         delta_flag = FALSE, 
         save_flag = TRUE, 
         width = 17, 
         height = 10,
         form = form, 
         bound = 1500 )

# lpd test
qoi_plot(dat = plotdf, 
         xlabel = "Prior",
         ylabel="LPD test" , 
         qoi = "lpd_test",  
         delta_flag = FALSE, 
         save_flag = TRUE,
         form = form, 
         bound = 5000 )

# lpd test delta
diff_qoi_plot(dat = difplotdf, 
              xlabel = "Prior",
              ylabel="ELPD" , 
              qoi = "lpd_test",  
              delta_flag = TRUE, 
              save_flag = TRUE, 
              plot_prefix = "dif",
              form = form, 
              bound = 500 )

diff_qoi_plot(dat = difplotdf, 
              xlabel = "Prior",
              ylabel="ELPD" , 
              qoi = "lpd_test_asinh",  
              delta_flag = TRUE, 
              save_flag = TRUE, 
              plot_prefix = "dif",
              form = form, 
              bound = 15 )


# RMSE ----
# Different types of RMSEs considered
# RMSE bpp nonzero

qoi_plot(dat = preddf, 
         xlabel = "Prior",
         ylabel="RMSE" , 
         qoi = "rmse_bppn0",  
         title = "Nonzero coefficients",
         delta_flag = FALSE, 
         form = form, bound = 10 )


diff_qoi_plot(dat = difplotdf, 
              xlabel = "Prior",
              ylabel="RMSE" , 
              qoi = "rmse_bppn0",  
              title = "Nonzero coefficients",
              delta_flag = TRUE, 
              form = form, bound = 5000 )

# RMSE bpp zeros

qoi_plot(dat = preddf, 
         xlabel = "Prior",
         ylabel="RMSE" , 
         qoi = "rmse_bpp0",  
         title = "Zero coefficients",
         delta_flag = FALSE, 
         form = form, bound = 5000 )


diff_qoi_plot(dat = difplotdf,
              xlabel = "Prior",
              ylabel="RMSE" , 
              qoi = "rmse_bpp0",  
              title = "Zero coefficients",
              delta_flag = TRUE, 
              form = form, bound = 3 )

# RMSE bpp

qoi_plot(dat = preddf, ylabel="RMSE" , 
         xlabel = "Prior",
         qoi = "rmse_bpp",  
         title = "All coefficients",
         delta_flag = FALSE, 
         form = form, bound = 3 )

diff_qoi_plot(dat = difplotdf, 
              xlabel = "Prior",
              ylabel="RMSE" , 
              qoi = "rmse_bpp",  
              title = "All coefficients",
              delta_flag = TRUE, 
              form = form, bound = 3 )


# delta rmse b partitioned and by coef type

coef_types <- levels(as.factor(sim_conds$type))

for( ctype in coef_types){
  
  difplotdf  %>%
    filter(str_detect(var, "^rmse_b$|^rmse_b0$|^rmse_bn0$")) %>% 
    filter( type == ctype ) %>% 
    mutate( var = recode(var, 
                         rmse_b = "all",
                         rmse_b0 = "zeroes",
                         rmse_bn0 = "nonzeroes")) %>% 
    ggplot(aes(x= quantity, y= value, fill= quantity))+
    geom_violin() +
    geom_boxplot(width = 0.25, outlier.size = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 1)+
    scale_fill_paletteer_d(my_palette)+
    geom_hline( yintercept = 0,  
                linetype = 2, 
                color= my_hline_color ) +
    ggh4x::facet_nested( var ~ p + R2,
                         labeller = labeller(R2 = label_R2, 
                                             p = label_p, 
                                             n = label_n
                                             #gencoef_function = label_gcf 
                         ),,
                         remove_labels = "all",
                         independent = "y",
                         scales = "free"
    ) +
    labs(
      y = bquote(Delta ~ "RMSE"),
      x = "Prior")+
    ggtitle( label_gcf[as.character(ctype)])+
    paper_theme  +
    theme(legend.position="none", 
          panel.spacing.y = unit(0.5, "lines"),
          axis.text.x = element_text(angle= 45, hjust = 0.5)
          #plot.margin = margin(t = 10, r = 10, b = 30, l = 10) 
    )
  
  
  fname= paste0(results_folder,"/",
                paste0("delta_rmse_postb", 
                       experiment_name, 
                       "_", 
                       ctype,
                       ".pdf",sep=""))  
  
  ggsave(filename=fname, 
         plot = last_plot(),
         width = 17, 
         height= 11)
  
  
  
}

# Coverage Tables ----

# Create summarydf

summarydf <- process_all_directories_parallel( path_results =  experiment_name ,
                                               result_type = "summary" ) 

coverage_stats <- summarydf %>%
  filter(grepl("^b", variable)) %>% 
  mutate(coverage = as.numeric(real_values < q97.5 & real_values > q2.5)) %>% 
  mutate( ci_width = q97.5-q2.5)


# Coverage: proportion, average width
coverage_summary <- summarydf %>%
  filter(grepl("^beta", variable)) %>% 
  #filter(id_model %in% id_model_vector) %>% 
  mutate(coverage = as.numeric(real_values < q97.5 & real_values > q2.5)) %>% 
  mutate( ci_width = q97.5-q2.5) %>% 
  group_by(id, p, type, R2, id_model, gencoef_function) %>%
  summarize(
    coverage_rate = mean(coverage, na.rm = TRUE),
    coverage_sd = sd(coverage, na.rm = TRUE),
    average_ci_width = mean(ci_width, na.rm = TRUE),
    average_ci_width_sd = sd(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

coverage_non_zero <- summarydf %>%
  filter(grepl("^beta", variable)) %>% 
  filter(real_values != 0) %>%  # Condition on theta ≠ 0
  #filter(id_model %in% id_model_vector) %>% 
  mutate(coverage = as.numeric(real_values < q97.5 & real_values > q2.5)) %>% 
  mutate( ci_width = q97.5-q2.5) %>% 
  group_by(id, p, type, R2, id_model, gencoef_function) %>%
  summarize(
    coverage_nonzero_rate = mean(coverage, na.rm = TRUE),
    coverage_nonzero_sd = sd(coverage, na.rm = TRUE),
    average_nonzero_ci_width = mean(ci_width, na.rm = TRUE),
    average_nonzero_ci_width_sd = sd(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

coverage_zero <- summarydf %>%
  filter(grepl("^beta", variable)) %>% 
  filter(real_values == 0) %>%  # Condition on theta ≠ 0
  #filter(id_model %in% id_model_vector) %>% 
  mutate(coverage = as.numeric(real_values < q97.5 & real_values > q2.5)) %>% 
  mutate( ci_width = q97.5-q2.5) %>% 
  group_by(id, p, type, R2, id_model, gencoef_function) %>%
  summarize(
    coverage_zero_rate = mean(coverage, na.rm = TRUE),
    coverage_zero_sd = sd(coverage, na.rm = TRUE),
    average_zero_ci_width = mean(ci_width, na.rm = TRUE),
    average_zero_ci_width_sd = sd(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

# Type I: Given that it is zero and is declared not as zero
# Specificity = 1- Type I
coverage_specificity <- summarydf %>%
  filter(grepl("^beta", variable)) %>% 
  filter(real_values == 0) %>%  # Condition on theta = 0
  #filter(id_model %in% id_model_vector) %>% 
  mutate(specificity = as.numeric(q2.5 <= 0 & q97.5 >= 0)) %>%  # CI contains 0 (fail to reject)
  mutate(ci_width = q97.5 - q2.5) %>% 
  group_by(id, p, type, R2, id_model, gencoef_function) %>%
  summarize(
    specificity_rate = mean(specificity, na.rm = TRUE),  # Compute specificity
    specificity_sd = sd(specificity, na.rm = TRUE),
    average_spec_ci_width = mean(ci_width, na.rm = TRUE),
    average_spec_ci_width_sd = sd(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

# Type II error: Given that it is not zero and is declared as such
# sensitivity = power = 1-type II error 
coverage_sensitivity <- summarydf %>%
  filter(grepl("^beta", variable)) %>% 
  filter(real_values != 0) %>%  # Condition on theta ≠ 0
  #filter(id_model %in% id_model_vector) %>% 
  mutate(sensitivity = as.numeric(!(q2.5 <= 0 & q97.5 >= 0))) %>%  # Check if 0 is inside and negate
  mutate(ci_width = q97.5 - q2.5) %>% 
  group_by(id, p, type, R2, id_model, gencoef_function) %>%
  summarize(
    sensitivity_rate = mean(sensitivity, na.rm = TRUE),  # Compute sensitivity (power)
    sensitivity_sd = sd(sensitivity, na.rm = TRUE),
    average_sens_ci_width = mean(ci_width, na.rm = TRUE),
    averages_sens_ci_width_sd = sd(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

# Combine summary statistics into a single dataframe
summary_table <- coverage_summary %>%
  left_join(coverage_specificity, by = c("id",  "p", "type", "R2", "id_model", "gencoef_function")) %>%
  left_join(coverage_sensitivity, by = c("id", "p", "type", "R2", "id_model", "gencoef_function")) %>%
  left_join(coverage_non_zero, by = c("id", "p", "type", "R2", "id_model", "gencoef_function")) %>%
  left_join(coverage_zero, by = c("id", "p", "type", "R2", "id_model", "gencoef_function")) %>%
  select(id, p, type, id_model, gencoef_function, 
         coverage_rate, specificity_rate, sensitivity_rate, 
         coverage_nonzero_rate, coverage_zero_rate,
         average_ci_width, average_zero_ci_width, average_nonzero_ci_width) %>%
  rename("Coverage" = coverage_rate,
         "Specificity" = specificity_rate,
         "Sensitivity" = sensitivity_rate,
         "Coverage Zero" = coverage_zero_rate,
         "Coverage Nonzero" = coverage_nonzero_rate,
         "Avg. CI Width" = average_ci_width,
         "Avg. Zero CI Width"= average_zero_ci_width,
         "Avg. Nonzero CI Width"= average_nonzero_ci_width,) %>%
  mutate(across(c(Coverage, Specificity, `Sensitivity`, `Coverage Nonzero` , `Coverage Zero`,
                  `Avg. CI Width`, `Avg. Zero CI Width`, `Avg. Zero CI Width` ), round, 3))


# Split summary_table by `id` and apply function
sm_tab_fixed <- summary_table %>%
  filter( gencoef_function == "gen_coef_fixed3" ) %>% 
  group_split(id) %>%
  map(create_latex_table)


sm_tab_fixed_filename <- file.path(results_folder,
                                   paste0( "tab_fixed_", experiment_name, ".RDS"))


saveRDS(sm_tab_fixed, file = sm_tab_fixed_filename)

sm_tab_random <- summary_table %>%
  filter( gencoef_function == "gen_coef_random" ) %>% 
  group_split(id) %>%
  map(create_latex_table)


sm_tab_random_filename <- file.path(results_folder,
                                    paste0( "tab_random_", experiment_name, ".RDS"))

saveRDS(sm_tab_random, file = sm_tab_random_filename)



# Generate LaTeX table 
# condition 

# ROC curves ----

my_palette <- "LaCroixColoR::paired"

quantile_pairs <- list(
  c("q1", "q99"),
  c("q2.5", "q97.5"),
  c("q5", "q95"),
  c("q10", "q90"),
  c("q15", "q85"),
  c("q20", "q80"),
  c("q25", "q75"),
  c("q33", "q67"),
  c("q35", "q65"),
  c("q40", "q60"),
  c("q50", "q50")  # Median-based decision
)

# Compute ROC data for each threshold
roc_data <- bind_rows(lapply(quantile_pairs, function(q) {
  compute_roc_data(summarydf, q[1], q[2]) %>%
    mutate(threshold = paste(q[1], "-", q[2]))
}))

# Plot ROC

roc_data <- roc_data %>% 
  mutate(id_model = recode_factor(id_model, 
                                            "D2def1g" = "D2gd1", 
                                            "D2def2g" = "D2gd2", 
                                            "D2default" = "R2-w", 
                                            "D2new" = "R2-c",
                                            "D2new2" = "R2-d",
                                            "D2new3" = "R2-0.1",
                                            "D2new4" = "R2-0.5",
                                            "D2new5" = "R2-1.0",
                                            "D2new6" = "R2-5.0",
                                            "D2newg" = "GR2-c",
                                            "D2new1g" = "GR2-c",
                                            "D2new2g" = "GR2-d",
                                            "D2new3g" = "GR2-0.1",
                                            "D2new4g" = "GR2-0.5",
                                            "D2new5g" = "GR2-1.0",
                                            "D2new6g" = "GR2-5.0",
                                            "D2unif" = "R2-u", 
                                            "D2unifg" = "GR2-u", 
                                            "gigg"= "GIGG", 
                                            "gigghier"= "GIGGh", 
                                            "giggnh" = "GIGGnh", 
                                            "giggnn" = "GIGGnn",
                                            "BPhier" = "BPh"
)) %>% 
  filter( !id_model %in% c( "R2-w", "R2-5.0", "BP", "GIGG", 
                            "GR2-5.0", "BPh", "GIGGh",  "RHS"
  )) 


proc <- roc_data %>%
  mutate(
    grouping = case_when(
      grepl("^GR2-", id_model) ~ "grouped",
      grepl("^R2-", id_model) ~ "nongrouped",
      TRUE ~ NA_character_
    ),
    base_model = paste0("R2-", sub("^(GR2-|R2-)", "", id_model))
  ) %>%
  filter(!is.na(grouping)) %>%  # optional: exclude models like BP, RHS, etc.
  group_by(p, type, R2, id_model) %>%
  ggplot(aes(
    x = FPR, y = TPR,
    color = base_model,
    linetype = grouping
  )) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_color_paletteer_d(my_palette) +
  scale_linetype_manual(values = c("grouped" = "dashed", "nongrouped" = "solid")) +
  ggh4x::facet_nested(
    R2 ~ type,
    labeller = labeller(
      R2 = label_R2,
      type = label_gcf
    ),
    independent = "y",
    scales = "free"
  ) +
  labs(
    title = "ROC curves",
    x = "False Positive Rate (1 - specificity)",
    y = "True Positive Rate (sensitivity)",
    color = "Prior", 
    linetype = "Version"
  ) +
  paper_theme +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18)
  )


proc

roc_filename <- file.path( paste0( "roc", experiment_name, ".pdf"))

# Save the plot
save_plot(proc, roc_filename, width = 20, height = 7)


# MCMC ----
 
 qoi_plot(dat = plotdf, 
          xlabel = "prior",
          ylabel="time" , 
          qoi = "time",  
          title = "",
          delta_flag = FALSE, 
          plot_prefix = "HMC",
          form = form, bound = 500 )










