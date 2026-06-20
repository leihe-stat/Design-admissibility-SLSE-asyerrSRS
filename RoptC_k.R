############################################################
# Batch computation for multiple k and t
############################################################

rm(list = ls())
library(nleqslv)

# ==========================================================
# Basic functions
# ==========================================================

tau_fun <- function(a, b, k) {
  a + (k - 1) * b - k * a^2
}

lambda_fun <- function(a, b, k) {
  a + (k - 2) * b - (k - 1) * a^2
}

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

# ==========================================================
# Feasibility conditions
# ==========================================================

check_conditions <- function(a, b, k) {
  
  cond1 <- (0 < b)
  cond2 <- (b < a)
  cond3 <- (a < 1)
  cond4 <- (k * a^2 < a + b * (k - 1))
  
  all(cond1, cond2, cond3, cond4)
}

# ==========================================================
# Equation system
# ==========================================================

system_eq <- function(x, t, k) {
  
  a <- x[1]
  b <- x[2]
  
  c(
    h0_1_fun(a, b, t, k),
    h2_1_fun(a, b, t, k)
  )
}

# ==========================================================
# Multi-start search
# ==========================================================

solve_ab_search <- function(t,
                            k,
                            n_start = 300) {
  
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

# ==========================================================
# Solve weights
# ==========================================================

solve_weights <- function(r, a, b, k) {
  
  A <- matrix(
    c(
      1, 1, 1,
      r[1], r[2], r[3],
      r[1]^2, r[2]^2, r[3]^2
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  rhs <- c(
    1,
    k * a,
    k * (k - 1) * b + k * a
  )
  
  solve(A, rhs)
}

# ==========================================================
# Find minimum support design
# ==========================================================

find_best_design <- function(a, b, k) {
  
  combin <- combn(k + 1, 3) - 1
  
  candidate_list <- list()
  
  for(j in 1:ncol(combin)) {
    
    r <- combin[, j]
    
    w <- try(
      solve_weights(r, a, b, k),
      silent = TRUE
    )
    
    if(inherits(w, "try-error"))
      next
    
    if(any(w <= 0) || any(w >= 1))
      next
    
    support_num <- sum(
      choose(k, r) * 2^r
    )
    
    candidate_list[[length(candidate_list)+1]] <- list(
      r = r,
      w = w,
      support_num = support_num
    )
  }
  
  if(length(candidate_list) == 0)
    return(NULL)
  
  support_sizes <- sapply(candidate_list,
                          function(x) x$support_num)
  
  best_id <- max(which(support_sizes == min(support_sizes)))
  #best_id <- which.min(support_sizes)
  
  candidate_list[[best_id]]
}

# ==========================================================
# Main batch computation
# ==========================================================

k_values <- 2:6
t_values <- c(0, 0.5, 0.9)

all_results <- list()

counter <- 1

for(k in k_values) {
  
  for(t in t_values) {
    
    cat("\n=================================================\n")
    cat("k =", k, ", t =", t, "\n")
    cat("=================================================\n")
    
    ab_res <- solve_ab_search(t, k)
    
    if(is.null(ab_res)) {
      
      cat("No feasible (a,b) found.\n")
      next
    }
    
    best_design <- find_best_design(
      a = ab_res$a,
      b = ab_res$b,
      k = k
    )
    
    if(is.null(best_design)) {
      
      cat("No feasible support design.\n")
      next
    }
    
    res <- data.frame(
      k = k,
      t = t,
      a = round(ab_res$a, 6),
      b = round(ab_res$b, 6),
      
      r1 = best_design$r[1],
      r2 = best_design$r[2],
      r3 = best_design$r[3],
      
      w1 = round(best_design$w[1], 6),
      w2 = round(best_design$w[2], 6),
      w3 = round(best_design$w[3], 6),
      
      support_num = best_design$support_num,
      
      h3_abs = signif(ab_res$h3_abs, 4)
    )
    
    print(res)
    
    all_results[[counter]] <- res
    
    counter <- counter + 1
  }
}

# ==========================================================
# Final summary table
# ==========================================================

final_results <- do.call(rbind, all_results)

cat("\n\n=====================================\n")
cat("FINAL RESULTS\n")
cat("=====================================\n")
final_results$w1 <- round(final_results$w1, 4)
final_results$w2 <- round(final_results$w2, 4)
final_results$w3 <- round(final_results$w3, 4)
print(final_results)