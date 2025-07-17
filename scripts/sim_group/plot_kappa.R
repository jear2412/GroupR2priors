# Load libraries
library(stats)
library(ggplot2)
library(gsl)
library(future)
library(furrr)

# Parameters
aG <- 0.5
pg <- 2
b <- 0.1

alpha <- aG - pg / 2
beta <- aG + b
phis <- c(0.5, 0.5 )


# Analytical integral function using confluent hypergeometric function U
compute_integral_analytic <- function(c_val, alpha, beta, log_flag = TRUE) {
  lambda <- beta-alpha
  gamma_lambda <- gamma(lambda)
  z <- c_val / 2
  U_val <- hyperg_U(lambda, 1 -alpha, z)
  if(log_flag){
    result <- log(gamma_lambda)+log(U_val)
  }else{
    result <- gamma_lambda*U_val
  }
  return(result)
}

compute_integral_numerical <- function(c_val, alpha, beta) {
  integrand <- function(t) {
    t^(alpha - 1) * (1 + t)^(-beta) * exp(-c_val / (2 * t))
  }
  result <- tryCatch({
    integrate(integrand, lower = 0, upper = Inf)$value
  }, error = function(e) NA)
  return(result)
}



# Define grid
b_seq <- seq(-1, 1, length.out = 200)
grid <- expand.grid(b1 = b_seq, b2 = b_seq)

# Compute c for each grid point
grid$c <- grid$b1^2 / phis[1] + grid$b2^2 / phis[2]

# Compute analytical integral for each c
grid$integral <- (sapply(grid$c, function(c_val) {
  compute_integral_analytic(c_val, alpha = alpha, beta = beta, log_flag = FALSE)
}))

lapis_lazuli <- "#26619c"

p1 <- ggplot(grid, aes(x = b1, y = b2, z = integral)) +
  geom_contour(bins = 20, color = lapis_lazuli) +
  #geom_contour(color = lapis_lazuli) +
  labs(
    title = "Marginal joint priors log scale",
    x = expression(b[1]),  # Set the x-axis label to b1
    y = expression(b[2])   # Set the y-axis label to b2
  ) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.title = element_text(size = 18),  # Axis label size
        axis.text = element_text(size = 13) ,   # Axis number size
        strip.text = element_text(size = 18),
        plot.title = element_text(size = 20), 
        panel.spacing = unit(0.75, "lines")
  ) 

p1



# Montecarlo approximation to p(b_g)

# Number of Monte Carlo samples
n_samples <- 1000

# Sample phis from Dirichlet
sample_phis <- function(a, size) {
  g1 <- rgamma(size, shape = 1)
  g2 <- rgamma(size, shape = 0.5)
  sum_g <- g1 + g2
  phis <- cbind(g1 / sum_g, g2 / sum_g)
  return(phis)
}

# Get Dirichlet samples
dirichlet_samples <- sample_phis(aG/2, n_samples)

plan(multisession, workers = 12) 

# Parallel computation of marginal p(b1, b2) using Monte Carlo with future
grid$marginal_mc <- future_map_dbl(1:nrow(grid), function(i) {
  b1 <- grid$b1[i]
  b2 <- grid$b2[i]
  mean(sapply(1:n_samples, function(j) {
    phi <- dirichlet_samples[j, ]
    c_val <- b1^2 / phi[1] + b2^2 / phi[2]
    compute_integral_analytic(c_val, alpha = alpha, beta = beta, log_flag = FALSE)
  }))
}, .progress = TRUE)  # .progress = TRUE for progress bar!


p2 <- ggplot(grid, aes(x = b1, y = b2, z = log(marginal_mc))) +
  geom_contour(bins = 20, color = lapis_lazuli) +
  labs(
    title = "Marginal joint priors log scale",
    x = expression(b[1]),  # Set the x-axis label to b1
    y = expression(b[2])   # Set the y-axis label to b2
  ) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.title = element_text(size = 18),  # Axis label size
        axis.text = element_text(size = 13) ,   # Axis number size
        strip.text = element_text(size = 18),
        plot.title = element_text(size = 20), 
        panel.spacing = unit(0.75, "lines")
  ) 

p2



