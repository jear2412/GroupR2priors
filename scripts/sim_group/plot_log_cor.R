# Load necessary packages
library(ggplot2)
library(dplyr)
library(MCMCpack)
library(paletteer)

# r = corr(log( phi_g varphi_gj), log(phi_g varphi_gk))

# Montecarlo ----

# one observation of r
simulate_r <- function(a_g, c_g, p_g = 5, n_sim = 1000, j = 1, k = 2, G = 10, log_flag = TRUE) {
  
  log_j <- numeric(n_sim)
  log_k <- numeric(n_sim)
  
  for (i in seq_len(n_sim)) {
    
    phi <- as.numeric(MCMCpack::rdirichlet(1, rep(a_g, G))) # Sample phi_g ~ Dir(a_g, ..., a_g)
    g <- sample(seq_len(G), 1) # Choose one group g at random
    
    varphi <- as.numeric(MCMCpack::rdirichlet(1, rep(c_g, p_g))) # Sample varphi_g ~ Dir(c_g, ..., c_g), length p_g
    
    # Store logs
    log_j[i] <- log(phi[1] * varphi[j])
    log_k[i] <- log(phi[1] * varphi[k])
  }
  
  if(log_flag){
    return( r = cor(log_j, log_k))  
  }else{
    return( r = cor(exp(log_j), exp(log_k) ))
  }
  
}


# many observations of r
simulate_r_distribution <- function(a_g, c_g, p_g = 5, j = 1, k = 2, G = 5, 
                                    n_sim_inner = 1000, n_rep = 100, log_flag = TRUE) {
  r <- replicate(n_rep, {
    simulate_r(a_g = a_g, 
               c_g= c_g, 
               p_g = p_g, 
               n_sim = n_sim_inner, 
               G = G, 
               log_flag = log_flag)
  })
  
  r <- data.frame(r = r, a_g = a_g, c_g = c_g, G = G, log_flag = log_flag)
  
}

# Define hyperparameter combinations
params <- expand.grid(a_g = c(0.01, 0.1, 0.5 , 1, 5, 10), 
                      c_g = c(0.01, 0.1, 0.5, 1, 5, 10), 
                      log_flag = c(TRUE, FALSE))


r_data <- purrr::map_dfr(1:nrow(params), function(i) {
  
  with(params[i, ], simulate_r_distribution(a_g, c_g,  p_g = 5, log_flag = log_flag))
  
})

table( r_data$log_flag ) #should be the same

# Custom labellers 
ag_labeller <- function(x) paste0("a[g] == ", x)
cg_labeller <- function(x) paste0("c[g] == ", x)

# plot that considers the log
p_log <- r_data %>% 
  dplyr::filter( log_flag == TRUE) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(bins = 25, color = "black", fill = "white") +
  ggh4x::facet_nested( a_g ~ c_g , 
                      labeller = labeller(
                        a_g = ag_labeller,
                        c_g = cg_labeller,
                        .default = label_parsed),
                      remove_labels = "all",
                      independent = "y",
                      scales = "free") +
  theme_minimal()+ 
  theme( 
    axis.text.y = element_blank(), 
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 20)
    #axis.text.y = element_text(size = 15)
  ) + 
  labs(title = "Prior distribution on correlations",
       x = expression(r == cor(log(phi[g] * varphi[gj]), log(phi[g] * varphi[gk]))),
       y = NULL)

p_log
ggsave("scripts/sim_group/plots/grouped_r2_log_cor.pdf", plot = p_log, width = 20, height = 9)

# plot that considers raw correlations
p_non_log <- r_data %>% 
  dplyr::filter( log_flag == FALSE) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(bins = 25, color = "black", fill = "white") +
  ggh4x::facet_nested( a_g ~ c_g , 
                       labeller = labeller(
                         a_g = ag_labeller,
                         c_g = cg_labeller,
                         .default = label_parsed),
                       remove_labels = "all",
                       independent = "y",
                       scales = "free") +
  theme_minimal()+ 
  theme( 
    axis.text.y = element_blank(), 
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 20)
    #axis.text.y = element_text(size = 15)
  ) + 
  labs(title = "Prior distribution on correlations",
       x = expression(r == cor( phi[g] * varphi[gj], phi[g] * varphi[gk])),
       y = NULL)

p_non_log

ggsave("scripts/sim_group/plots/grouped_r2_nonlog_cor.pdf", plot = p_non_log, width = 20, height = 9)



# Analytical ----


corr_logprod <- function(a_g, c_g, G = 10, p_g = 5) {
  num <- trigamma(a_g) - trigamma(G * a_g) - trigamma(p_g * c_g)
  den <- trigamma(a_g) - trigamma(G * a_g) + trigamma(c_g) - trigamma(p_g * c_g)
  return(num / den)
}

corr_prod <- function(a_g, c_g, G = 10, p_g = 5) {
  term1 <- ((a_g + 1) * c_g) / (G *(G * a_g + 1) * p_g * (c_g * p_g + 1))
  term2 <- 1 / (G^2 * p_g^2)

  num <- term1-term2
  
  term4 <- ( (a_g + 1)) / (G * (G * a_g + 1))
  term5 <- ( (c_g + 1)) / (p_g * (p_g * c_g + 1))
  term6 <- 1 / (G^2 * p_g^2)
  
  den <- term4 * term5 - term6
  return(num / den)
}

my_palette_d <- "LaCroixColoR::paired"

a_g_vals <- c(0.1, 0.25, 0.5, 1, 3, 5)
c_g_seq <- seq(0.01, 5, length.out = 500)

# line plots
plot_correlation_curves <- function(G = 10, p_g = 5) {
  
  df <- expand.grid(a_g = a_g_vals,
                    c_g = c_g_seq) %>%
    mutate(corr = mapply(corr_logprod, a_g, c_g,
                         MoreArgs = list(G = G, p_g = p_g)))
  
  ggplot(df, aes(x = c_g, y = corr, color = factor(a_g))) +
    geom_line(size = 1) +
    paletteer::scale_color_paletteer_d(my_palette_d)+
    labs(title = "Correlation of marginal log variances",
         x = expression(c[g]),
         y = "Correlation",
         color = expression(a[g])) +
    theme_minimal()+ 
    theme( 
      #axis.text.y = element_blank(), 
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18), 
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18), 
      plot.title = element_text(size = 20)
    )

}

plot_corr_product_curves <- function(G = 10, p_g = 5) {
  
  df <- expand.grid(a_g = a_g_vals, c_g = c_g_seq) %>%
    mutate(corr = mapply(corr_prod, a_g, c_g, MoreArgs = list(G = G, p_g = p_g)))
  
  ggplot(df, aes(x = c_g, y = corr, color = factor(a_g))) +
    geom_line(size = 1.2) +
    paletteer::scale_color_paletteer_d(my_palette_d)+
    labs(title = "Correlation of marginal variances",
         x = expression(c[g]),
         y = "Correlation",
         color = expression(a[g])) +
    theme_minimal()+ 
    theme( 
      #axis.text.y = element_blank(), 
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18), 
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18), 
      plot.title = element_text(size = 20)
    )

  }


# Contour plots
plot_correlation_contour <- function(G = 10, p_g = 5) {
  
  a_g_seq <- seq(0.05, 5, length.out = 100)
  c_g_seq <- seq(0.05, 5, length.out = 100)
  
  df <- expand.grid(a_g = a_g_seq, 
                    c_g = c_g_seq) %>%
    mutate(corr = mapply(corr_logprod,
                         a_g, c_g, 
                         MoreArgs = list(G = G, p_g = p_g)))
  
  ggplot(df, aes(x = c_g, y = a_g, z = corr)) +
    geom_contour_filled(bins = 15) + 
    scale_fill_viridis_d(option = "D") +
    labs(title = "Contour plot of log correlation",
         x = expression(c[g]),
         y = expression(a[g]),
         fill = "Correlation") +
    theme_minimal()+ 
    theme( 
      axis.text.y = element_blank(), 
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 15),
      #axis.text.y = element_text(size = 15)
    )

}

p1 <- plot_correlation_curves(G = 10, p_g = 5)
p1

ggsave("scripts/sim_group/plots/log_cor_analytical.pdf", plot = p1, width = 14, height = 4)



p2 <- plot_corr_product_curves(G = 10, p_g = 5)
p2

p3 <- plot_correlation_contour( G = 10, p_g = 10)
p3
ggsave("scripts/sim_group/plots/cor_analytical.pdf", plot = p1, width = 14, height = 4)


# Montecarlo all pairwise

# One simulation of all pairwise correlations
simulate_r_all_pairs <- function(a_g, c_g, p_g = 5, G = 10, n_sim = 100, log_flag = TRUE) {
  mat_log <- matrix(NA, nrow = n_sim, ncol = p_g)
  
  for (i in seq_len(n_sim)) {
    phi <- as.numeric(MCMCpack::rdirichlet(1, rep(a_g, G)))
    g <- sample(seq_len(G), 1)
    varphi <- as.numeric(MCMCpack::rdirichlet(1, rep(c_g, p_g)))
    x <- phi[g] * varphi
    mat_log[i, ] <- if (log_flag) log(x) else x
  }
  
  cor_matrix <- cor(mat_log)
  
  # Extract upper triangle without diagonal
  cor_values <- cor_matrix[upper.tri(cor_matrix)]
  return(cor_values)
}

# Repeat simulation multiple times and store all pairwise correlations
simulate_r_all_distribution <- function(a_g, c_g, p_g = 5, G = 10,
                                        n_sim_inner = 1000, n_rep = 100,
                                        log_flag = TRUE) {
  all_corrs <- replicate(n_rep, {
    simulate_r_all_pairs(a_g, c_g, p_g, G, n_sim_inner, log_flag)
  })
  
  # Reshape into tidy format
  df <- data.frame(
    r = as.vector(all_corrs),
    rep = rep(seq_len(n_rep), each = choose(p_g, 2)),
    a_g = a_g,
    c_g = c_g,
    G = G,
    log_flag = log_flag
  )
  
  return(df)
}


# Define hyperparameter combinations
params <- expand.grid(a_g = c(0.01, 0.1, 0.5 , 1, 5, 10), 
                      c_g = c(0.01, 0.1, 0.5, 1, 5, 10), 
                      log_flag = c(TRUE, FALSE))


r_all_data <- purrr::map_dfr(1:nrow(params), function(i) {
  
  with(params[i, ], simulate_r_all_distribution(a_g, c_g,  p_g = 5, G = 10,
                                                log_flag = log_flag))
  
})

table( r_all_data$log_flag ) #should be the same

# Custom labellers 
ag_labeller <- function(x) paste0("a[g] == ", x)
cg_labeller <- function(x) paste0("c[g] == ", x)

# plot that considers the log
p_log_all <- r_all_data %>% 
  dplyr::filter( log_flag == TRUE) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(bins = 25, color = "black", fill = "white") +
  ggh4x::facet_nested( a_g ~ c_g , 
                       labeller = labeller(
                         a_g = ag_labeller,
                         c_g = cg_labeller,
                         .default = label_parsed),
                       remove_labels = "all",
                       independent = "y",
                       scales = "free") +
  theme_minimal()+ 
  theme( 
    axis.text.y = element_blank(), 
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 20)
    #axis.text.y = element_text(size = 15)
  ) + 
  labs(title = "Prior distribution on correlations",
       x = expression(r == cor(log(phi[g] * varphi[gj]), log(phi[g] * varphi[gk]))),
       y = NULL)

p_log_all

ggsave("scripts/sim_group/plots/grouped_r2_log_cor.pdf", plot = p_log, width = 20, height = 9)

# plot that considers raw correlations
p_non_log_all <- r_all_data %>% 
  dplyr::filter( log_flag == FALSE) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(bins = 25, color = "black", fill = "white") +
  ggh4x::facet_nested( a_g ~ c_g , 
                       labeller = labeller(
                         a_g = ag_labeller,
                         c_g = cg_labeller,
                         .default = label_parsed),
                       remove_labels = "all",
                       independent = "y",
                       scales = "free") +
  theme_minimal()+ 
  theme( 
    axis.text.y = element_blank(), 
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 20)
    #axis.text.y = element_text(size = 15)
  ) + 
  labs(title = "Prior distribution on correlations",
       x = expression(r == cor( phi[g] * varphi[gj], phi[g] * varphi[gk])),
       y = NULL)

p_non_log_all

ggsave("scripts/sim_group/plots/grouped_r2_nonlog_cor.pdf", plot = p_non_log, width = 20, height = 9)




