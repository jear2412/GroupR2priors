library(gtools)
library(dplyr)
library(tidyr)
library(ggplot2)


# ggplot2 settings

my_palette <- "LaCroixColoR::paired"

paper_theme <- theme_bw()+ 
  theme(plot.title = element_text(size=23),
        #axis.line.y = element_blank(),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 19),
        strip.text.y = element_text(size = 19),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        #axis.text.x = element_text(angle= 90, vjust= 0.5, hjust= 1, size = 17),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

my_hline_color <- "#9a031e" #nice color for vertical or horizontal lines


# Settings
sigma2 <- 1
group_size <- 10
G <- 20 # number of groups
a_G_vals <- c(0.1, 0.5, 1, 5)
c_g_vals <- c(0.01, 0.1, .5, 1, 5 )
n_samples <- 5000

r2_means <- G*a_G_vals/(G*a_G_vals+0.5)
r2_scales <- G*a_G_vals  + 0.5 
r2_group_mean <- a_G_vals/(a_G_vals+0.5)*1/G

r2_df <- data.frame( R2_mean =  r2_means,
                     R2_scale = r2_scales, 
                     R2group_mean = r2_group_mean)

print(round(r2_df, 2))


# Simulate meff for one group
simulate_meff_one_group <- function(a_G, c_g) {
  meff <- numeric(n_samples)
  for (i in 1:n_samples) {
    # tau^2 (G a_g, 0.5)
    # tau_g^2 bp(a_G)
    tau2 <- rgamma(1, a_G, 1) / rgamma(1, 0.5, 1)  # BetaPrime(a_G, 0.5)
    varphi <- as.numeric(rdirichlet(1, rep(c_g, group_size)))
    lambda2 <- varphi * tau2 *sigma2
    kappa <- 1 / (1 + lambda2)
    meff[i] <- sum(1 - kappa)
  }
  tibble(a_G = a_G, c_g = c_g, meff = meff)
}

# Simulate across (a_G, c_g) combinations
grid <- expand.grid(a_G = a_G_vals, c_g = c_g_vals)

results <- purrr::pmap_dfr(grid, simulate_meff_one_group) %>%
  mutate(
    a_G = factor(a_G),
    c_g = factor(c_g)
  )

# Custom labellers for math expressions
ag_labeller <- function(x) paste0("a[G] == ", x)
cg_labeller <- function(x) paste0("c[g] == ", x)


# Plot
p1 <- ggplot(results, aes(x = meff)) +
  geom_histogram(bins = 20, color = "black", fill = "white") +
  ggh4x::facet_nested(a_G ~ c_g,
                      labeller = labeller(
                        a_G = ag_labeller,
                        c_g = cg_labeller,
                        .default = label_parsed
                      ),,
                      remove_labels = "all",
                      #independent = "y",
                      independent = "y",
                      scales = "free"
  ) +
  labs(
    x = expression(m[eff * "," * g]~"(effective number of coefficients per group)"),
    y = ""
  ) + 
  theme_minimal()+ 
  #paper_theme +
  theme( 
  axis.text.y = element_blank(), 
  strip.text.x = element_text(size = 19),
  strip.text.y = element_text(size = 19),
  axis.title.x = element_text(size = 20),
  #axis.title.y = element_text(size = 18),
  axis.text.x = element_text(size = 16.5),
  #axis.text.y = element_text(size = 15)
  )

p1  
#
ggsave("scripts/sim_group/plots/meff_group.pdf", plot = p1, width = 14, height = 6)

