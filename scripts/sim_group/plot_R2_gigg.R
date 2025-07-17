# Load necessary packages
library(ggplot2)
library(dplyr)


# simulate r2 for one group from GIGG
simulate_r2_per_group <- function(a_g, c_g, p = 50, n_sim = 10000, tau2 = 1, sigma2 = 1) {
  r2_values <- replicate(n_sim, {
    gamma2 <- rgamma(1, shape = a_g, rate = 1) # gamma 
    lambda2 <- 1 / rgamma(p, shape = c_g, rate = 1) #inverse gamma
    sd_b <- sqrt(tau2 * gamma2 * lambda2) # scale of b
    b <- rnorm(p, mean = 0, sd = sd_b)
    r2 <- sum(b^2) / (sum(b^2) + sigma2) # r2 in normal means setting
    return(r2)
  })
  
  data.frame(R2 = r2_values, a_g = a_g, c_g = c_g)
}

# simulate r2 for all from GIGG

simulate_r2 <- function(a_g, c_g, p = 50, G = 10, n_sim = 10000, tau2 = 1, sigma2 = 1) {
  stopifnot(p %% G == 0)  # ensure equal group sizes
  p_g <- p / G            # number of coefficients per group
  
  r2_values <- replicate(n_sim, {
    b <- numeric(p)
    
    for (g in 1:G) {
      gamma2 <- rgamma(1, shape = a_g, rate = 1)         #  gamma2 ~ Gamma(a_g, 1)
      lambda2 <- 1 / rgamma(p_g, shape = c_g, rate = 1)  #  lambda2 ~ InvGamma(c_g, 1)
      sd_b <- sqrt(tau2 * gamma2 * lambda2)              # SD for b_gl
      b[((g - 1) * p_g + 1):(g * p_g)] <- rnorm(p_g, mean = 0, sd = sd_b)
    }
    
    r2 <- sum(b^2) / (sum(b^2) + sigma2)
    return(r2)
  })
  
  data.frame(R2 = r2_values, a_g = a_g, c_g = c_g, G = G)
}

# Define hyperparameter combinations
params <- expand.grid(a_g = c(0.1, 0.5 ), 
                      c_g = c(0.1, 0.5, 1))

r2_data <- purrr::map_dfr(1:nrow(params), function(i) {
  
  with(params[i, ], simulate_r2_per_group(a_g, c_g,  p = 50,  n_sim = 10000, tau2= 1, sigma2= 1))
  
})

# Custom labellers 
ag_labeller <- function(x) paste0("a[g] == ", x)
cg_labeller <- function(x) paste0("c[g] == ", x)


p1 <- ggplot(r2_data, aes(x = R2)) +
  geom_histogram(bins = 20, color = "black", fill = "white") +
  ggh4x::facet_nested(a_g ~ c_g, 
                      labeller = labeller(
                        a_g = ag_labeller,
                        c_g = cg_labeller,
                        .default = label_parsed
                      ),
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
  labs(title = "Implied Prior on R2 with the GIGG Prior",
       x = expression(R^2), 
       y = NULL)


p1  

ggsave("scripts/sim_group/plots/gigg_r2.pdf", plot = p1, width = 14, height = 6)



