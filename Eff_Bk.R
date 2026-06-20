############################################################
# R-optimal design algorithm on B_k
# Efficiency calculation and plotting
############################################################

rm(list = ls())

library(nleqslv)
library(ggplot2)
library(dplyr)
library(grid)

############################################################
# Basic functions
############################################################

delta_fun <- function(c, d, k) {
  c * (1 - k * c) + c^2 - d
}

############################################################
# Product of diagonal elements of M(c,d,t)
############################################################

diag_product_M_cd <- function(c, d, t, k) {
  
  delta <- delta_fun(c, d, k)
  
  # validity checks
  if(delta <= 0) return(NA)
  if(c <= 0 || d <= 0) return(NA)
  if(c - k * d <= 0) return(NA)
  if(1 - k * c <= 0) return(NA)
  if(t >= 1) return(NA)
  
  ##########################################################
  # Diagonal elements
  ##########################################################
  
  # first block
  d1 <- (1 - k * c * t) /
    ((1 - t) * (1 - k * c))
  
  # second block
  d2 <- (1 / (c - k * d)) *
    (1 + (c^2 - d) /
       (c * (1 - k * c)))
  
  # third block
  d3 <- 1 / c
  
  # fourth block
  d4 <- 1 / d
  
  ##########################################################
  # Product
  ##########################################################
  
  prod_diag <- d1 *
    (d2^k) *
    (d3^k) *
    (d4^(k * (k - 1) / 2))
  
  return(prod_diag)
}

############################################################
# g-functions
############################################################

g0_fun <- function(c, d, t, k) {
  
  delta <- delta_fun(c, d, k)
  
  k * c * (c - k * d) /
    (delta * (1 - k * c)) +
    
    (1 - k * c * t) /
    ((1 - t) * (1 - k * c))
}

g1_fun <- function(c, d, t, k) {
  
  delta <- delta_fun(c, d, k)
  
  -(delta + c - k * d) /
    (delta * (1 - k * c))
}

g2_fun <- function(c, d, t, k) {
  
  delta <- delta_fun(c, d, k)
  
  c * (1 - k * c) /
    (delta * (c - k * d))
}

g3_fun <- function(c, d, t, k) {
  
  delta <- delta_fun(c, d, k)
  
  2 * (c^2 - d) /
    (delta * (c - k * d)) +
    
    k * (c^2 - d)^2 /
    (delta * c * (1 - k * c) * (c - k * d)) +
    
    (1 - t) /
    ((1 - k * c) * (1 - k * c * t))
}

############################################################
# Derived functions
############################################################

g1_1_fun <- function(c, d, t, k) {
  
  g1 <- g1_fun(c, d, t, k)
  g2 <- g2_fun(c, d, t, k)
  g3 <- g3_fun(c, d, t, k)
  
  2 * g1 +
    1 / c -
    2 * t * (g1 + c * g2 + k * c * g3)
}

g0_1_fun <- function(c, d, t, k) {
  
  g1_1 <- g1_1_fun(c, d, t, k)
  
  -0.5 * c * k * (g1_1 + (k + 2) / c)
}

g2_1_fun <- function(c, d, t, k) {
  
  g2 <- g2_fun(c, d, t, k)
  g3 <- g3_fun(c, d, t, k)
  g1_1 <- g1_1_fun(c, d, t, k)
  
  g2 + g3 + g1_1
}

g3_1_fun <- function(c, d, t, k) {
  
  g3 <- g3_fun(c, d, t, k)
  g1_1 <- g1_1_fun(c, d, t, k)
  
  2 * g3 + 1 / d + 2 * g1_1
}

############################################################
# Feasibility conditions
############################################################

check_conditions <- function(c, d, k) {
  
  cond1 <- (k * d < c)
  cond2 <- (c < 1 / k)
  
  all(cond1, cond2)
}

############################################################
# Equation system
############################################################

system_eq <- function(x, t, k) {
  
  c_val <- x[1]
  d_val <- x[2]
  
  c(
    g0_1_fun(c_val, d_val, t, k),
    g2_1_fun(c_val, d_val, t, k)
  )
}

############################################################
# Multi-start search
############################################################

solve_cd_search <- function(t,
                            k,
                            n_start = 50) {
  
  best <- NULL
  best_g3 <- Inf
  
  for(i in 1:n_start) {
    
    c0 <- runif(1, 0.02, 1 / k)
    d0 <- runif(1, 1e-6, c0 / k)
    
    sol <- try(
      nleqslv(
        c(c0, d0),
        system_eq,
        t = t,
        k = k,
        method = "Broyden",
        control = list(btol = 1e-10)
      ),
      silent = TRUE
    )
    
    if(inherits(sol, "try-error"))
      next
    
    if(sol$termcd != 1)
      next
    
    c_val <- sol$x[1]
    d_val <- sol$x[2]
    
    if(!check_conditions(c_val, d_val, k))
      next
    
    g3v <- abs(g3_1_fun(c_val, d_val, t, k))
    
    if(g3v < best_g3) {
      
      best_g3 <- g3v
      
      best <- list(
        c = c_val,
        d = d_val,
        g3_abs = g3v
      )
    }
  }
  
  return(best)
}

############################################################
# Main computation
############################################################

t_grid <- seq(0, 0.99, length = 100)

result_df <- data.frame()

for(k in 2:6) {
  
  cat("Processing k =", k, "\n")
  
  ##########################################################
  # Baseline solution at t = 0
  ##########################################################
  
  base_sol <- solve_cd_search(t = 0, k = k)
  
  if(is.null(base_sol)) {
    cat("No baseline solution for k =", k, "\n")
    next
  }
  
  c0 <- base_sol$c
  d0 <- base_sol$d
  
  ##########################################################
  # Loop over t
  ##########################################################
  
  for(t in t_grid) {
    
    sol <- solve_cd_search(t = t, k = k)
    
    if(is.null(sol))
      next
    
    c_val <- sol$c
    d_val <- sol$d
    
    prod_now <- diag_product_M_cd(c_val, d_val, t, k)
    
    prod_base <- diag_product_M_cd(c0, d0, t, k)
    
    Eff <- prod_now / prod_base
    
    result_df <- rbind(
      result_df,
      data.frame(
        k = factor(k),
        t = t,
        c = c_val,
        d = d_val,
        Eff = Eff
      )
    )
  }
}

############################################################
# Plot
############################################################

ggplot(result_df,
       aes(x = t,
           y = Eff,
           color = k,
           linetype = k)) +
  
  geom_line(size = 1) +
  
  scale_linetype_manual(
    values = c("solid",
               "dashed",
               "dotted",
               "dotdash",
               "longdash"),
    
    labels = c(
      expression(k == 2),
      expression(k == 3),
      expression(k == 4),
      expression(k == 5),
      expression(k == 6)
    )
  ) +
  
  ##########################################################
# Academic color palette
##########################################################

scale_color_manual(
  values = c(
    "#1F77B4",  # blue
    "#FF7F0E",  # orange
    "#2CA02C",  # green
    "#D62728",  # red
    "#9467BD"   # purple
  ),
  
  labels = c(
    expression(k == 2),
    expression(k == 3),
    expression(k == 4),
    expression(k == 5),
    expression(k == 6)
  )
) +
  
  labs(
    x = "t",
    y = "Efficiency",
    color = NULL,
    linetype = NULL
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    
    ########################################################
    # Legend inside lower-left corner
    ########################################################
    
    legend.position = c(0.18, 0.20),
    
    ########################################################
    # Clean legend appearance
    ########################################################
    
    legend.background = element_blank(),
    legend.key = element_blank(),
    
    ########################################################
    # Longer legend lines
    ########################################################
    
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.7, "cm"),
    
    ########################################################
    # Typography
    ########################################################
    
    legend.text = element_text(size = 12),
    
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  )

############################################################
# Optional: view results
############################################################

head(result_df)