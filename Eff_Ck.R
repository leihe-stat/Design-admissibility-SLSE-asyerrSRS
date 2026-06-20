############################################################
# Batch computation of Eff for k = 2,...,6
# and t in [0,0.99]
############################################################

rm(list = ls())

library(nleqslv)
library(ggplot2)
library(dplyr)

############################################################
# Basic functions
############################################################

tau_fun <- function(a, b, k) {
  a + (k - 1) * b - k * a^2
}

lambda_fun <- function(a, b, k) {
  a + (k - 2) * b - (k - 1) * a^2
}

############################################################
# Product of diagonal elements of M(a,b,t)
############################################################

diag_product_M <- function(a, b, t, k) {
  
  tau <- tau_fun(a, b, k)
  
  if(tau <= 0) return(NA)
  if(a <= 0 || b <= 0) return(NA)
  if(a <= b) return(NA)
  if(t >= 1) return(NA)
  
  # First block
  d1 <- (a - b + k * (b - a^2 * t)) /
    ((1 - t) * tau)
  
  # Second block
  d2 <- (1 / (a - b)) *
    (1 + (a^2 - b) / tau)
  
  # Third block
  d3 <- 1 / a
  
  # Fourth block
  d4 <- 1 / b
  
  # Product
  prod_diag <- d1 *
    (d2^k) *
    (d3^k) *
    (d4^(k * (k - 1) / 2))
  
  return(prod_diag)
}

############################################################
# h-functions
############################################################

h0_fun <- function(a, b, t, k) {
  
  lam <- lambda_fun(a, b, k)
  tau <- tau_fun(a, b, k)
  
  (a - b + b * k) / tau +
    a^2 * (a - b) * k / (lam * tau) +
    t / (1 - t)
}

h1_fun <- function(a, b, t, k) {
  
  tau <- tau_fun(a, b, k)
  lam <- lambda_fun(a, b, k)
  
  -a / tau -
    a / lam -
    a * (a^2 - b) * k / (lam * tau)
}

h2_fun <- function(a, b, t, k) {
  
  tau <- tau_fun(a, b, k)
  lam <- lambda_fun(a, b, k)
  
  tau / (lam * (a - b))
}

h3_fun <- function(a, b, t, k) {
  
  tau <- tau_fun(a, b, k)
  lam <- lambda_fun(a, b, k)
  
  a^2 * (1 - t) /
    ((a - b + b * k - a^2 * k * t) * tau) +
    
    2 * (a^2 - b) /
    (lam * (a - b)) +
    
    (a^2 - b)^2 * k /
    (lam * tau * (a - b))
}

h1_1_fun <- function(a, b, t, k) {
  
  h1 <- h1_fun(a, b, t, k)
  h2 <- h2_fun(a, b, t, k)
  h3 <- h3_fun(a, b, t, k)
  
  2 * h1 +
    1 / a -
    2 * t * (h1 + a * k * h3 + a * h2)
}

h0_1_fun <- function(a, b, t, k) {
  
  h1_1 <- h1_1_fun(a, b, t, k)
  
  -0.5 * a * k * (h1_1 + (k + 2) / a)
}

h2_1_fun <- function(a, b, t, k) {
  
  h2 <- h2_fun(a, b, t, k)
  h3 <- h3_fun(a, b, t, k)
  h1_1 <- h1_1_fun(a, b, t, k)
  
  h2 + h3 + h1_1
}

h3_1_fun <- function(a, b, t, k) {
  
  h3 <- h3_fun(a, b, t, k)
  
  2 * h3 + 1 / b
}

############################################################
# Feasibility conditions
############################################################

check_conditions <- function(a, b, k) {
  
  cond1 <- (0 < b)
  cond2 <- (b < a)
  cond3 <- (a < 1)
  cond4 <- (k * a^2 < a + b * (k - 1))
  
  all(cond1, cond2, cond3, cond4)
}

############################################################
# Equation system
############################################################

system_eq <- function(x, t, k) {
  
  a <- x[1]
  b <- x[2]
  
  c(
    h0_1_fun(a, b, t, k),
    h2_1_fun(a, b, t, k)
  )
}

############################################################
# Multi-start search
############################################################

solve_ab_search <- function(t,
                            k,
                            n_start = 50) {
  
  best <- NULL
  best_h3 <- Inf
  
  for(i in 1:n_start) {
    
    a0 <- runif(1, 0.1, 0.95)
    b0 <- runif(1, 0.01, a0 - 0.01)
    
    sol <- try(
      nleqslv(
        c(a0, b0),
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
    
    a <- sol$x[1]
    b <- sol$x[2]
    
    if(!check_conditions(a, b, k))
      next
    
    h3v <- abs(h3_1_fun(a, b, t, k))
    
    if(h3v < best_h3) {
      
      best_h3 <- h3v
      
      best <- list(
        a = a,
        b = b,
        h3_abs = h3v
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
  
  base_sol <- solve_ab_search(t = 0, k = k)
  
  if(is.null(base_sol)) {
    cat("No baseline solution for k =", k, "\n")
    next
  }
  
  a0 <- base_sol$a
  b0 <- base_sol$b
  
  base_prod <- diag_product_M(a0, b0, 0, k)
  
  ##########################################################
  # Loop over t
  ##########################################################
  
  for(t in t_grid) {
    
    sol <- solve_ab_search(t = t, k = k)
    
    if(is.null(sol))
      next
    
    a <- sol$a
    b <- sol$b
    
    prod_now <- diag_product_M(a, b, t, k)
    
    prod_base <- diag_product_M(a0, b0, t, k)
    
    Eff <- prod_now / prod_base
    
    result_df <- rbind(
      result_df,
      data.frame(
        k = factor(k),
        t = t,
        a = a,
        b = b,
        Eff = Eff
      )
    )
  }
}

############################################################
# Publication-style color palette specified by user
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
# User-specified academic color palette
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
# Optional: view data
############################################################

head(result_df)