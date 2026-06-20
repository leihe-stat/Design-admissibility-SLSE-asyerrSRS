############################################################
# R-optimal design algorithm on B_k (Corrected Version)
############################################################

rm(list = ls())
library(nleqslv)

# ==========================================================
# Basic functions
# ==========================================================

delta_fun <- function(c, d, k) {
  c * (1 - k * c) + c^2 - d
}

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

# ==========================================================
# Derived functions
# ==========================================================

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
  
  # 修正了原论文 Equation (A.3) 中遗漏的交叉项补偿 2 * g1_1
  2 * g3 + 1 / d + 2 * g1_1
}

# ==========================================================
# Feasibility conditions
# ==========================================================

check_conditions <- function(c, d, k) {
  
  cond1 <- (k * d < c)
  cond2 <- (c < 1 / k)
  all(cond1, cond2)
}

# ==========================================================
# Equation system
# ==========================================================

system_eq <- function(x, t, k) {
  
  c <- x[1]
  d <- x[2]
  
  c(
    g0_1_fun(c, d, t, k),
    g2_1_fun(c, d, t, k)
  )
}

# ==========================================================
# Multi-start search
# ==========================================================

solve_cd_search <- function(t, k, n_start = 300) {
  
  best <- NULL
  best_g3 <- Inf
  
  for(i in 1:n_start) {
    
    c0 <- runif(1, 0.02, 1/k)
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

# ==========================================================
# Compute weights directly from equation (14)
# ==========================================================

compute_weights <- function(c, d, k) {
  
  xi_Sk <- k^2 * d
  
  xi_S1 <- k * c - k^2 * d
  
  xi_S0 <- 1 - k * c
  
  c(
    xi_S0 = xi_S0,
    xi_S1 = xi_S1,
    xi_Sk = xi_Sk
  )
}

# ==========================================================
# Main batch computation
# ==========================================================

k_values <- 2:11
t_values <- c(0, 0.5, 0.9)

all_results <- list()

counter <- 1

for(k in k_values) {
  
  for(t in t_values) {
    
    cat("\n=================================================\n")
    cat("k =", k, ", t =", t, "\n")
    cat("=================================================\n")
    
    cd_res <- solve_cd_search(t, k)
    
    if(is.null(cd_res)) {
      
      cat("No feasible solution found.\n")
      next
    }
    
    weights <- compute_weights(
      c = cd_res$c,
      d = cd_res$d,
      k = k
    )
    
    # Check weights
    if(any(weights <= 0) || any(weights >= 1)) {
      
      cat("Weights not feasible.\n")
      next
    }
    
    res <- data.frame(
      
      k = k,
      t = t,
      
      c = round(cd_res$c, 6),
      d = round(cd_res$d, 6),
      
      xi_S0 = round(weights[1], 6),
      xi_S1 = round(weights[2], 6),
      xi_Sk = round(weights[3], 6),
      
      g3_abs = signif(cd_res$g3_abs, 4)
    )
    
    print(res)
    
    all_results[[counter]] <- res
    
    counter <- counter + 1
  }
}

# ==========================================================
# Final summary
# ==========================================================

final_results <- do.call(rbind, all_results)

cat("\n=====================================\n")
cat("FINAL RESULTS\n")
cat("=====================================\n")
final_results$xi_S0 <- round(final_results$xi_S0, 4)
final_results$xi_S1 <- round(final_results$xi_S1, 4)
final_results$xi_Sk <- round(final_results$xi_Sk, 4)
print(final_results)