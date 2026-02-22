library(MASS)
library(dplyr)
library(ggplot2)
library(latex2exp)

penalize_diag = FALSE
#define the KL loss
evaluate_kl <- function(Omega_est, Sigma_target) {
  
  p <- ncol(Sigma_target)
  # Log-determinant of Estimated Precision
  det_est <- determinant(Omega_est, logarithm = TRUE)
  if (det_est$sign != 1) return(Inf) 
  log_det_C <- det_est$modulus[1]
  
  # Trace Term 
  trace_term <- sum(diag(Omega_est %*% Sigma_target))
  

  # Entropy Term: 
  det_S <- determinant(Sigma_target, logarithm = TRUE)
  log_det_S <- det_S$modulus[1]
  entropy_term <- -(log_det_S + p)
  
  
  
  # Final KL
  kl_val <- -log_det_C + trace_term + entropy_term
  return(as.numeric(kl_val))
  
}



# Model Generation

generate_precision <- function(model_type, p) {
  Omega <- matrix(0, p, p)
  if (model_type == "Heterogeneous") { diag(Omega) <- 1:p }
  else if (model_type == "AR1") { for (i in 1:p) for (j in 1:p) { if (i == j) Omega[i, j] <- 1; if (abs(i - j) == 1) Omega[i, j] <- 0.5 } }
  else if (model_type == "AR2") { for (i in 1:p) for (j in 1:p) { if (i == j) Omega[i, j] <- 1; if (abs(i - j) == 1) Omega[i, j] <- 0.5; if (abs(i - j) == 2) Omega[i, j] <- 0.25 } }
  else if (model_type == "AR3") { for (i in 1:p) for (j in 1:p) { if (i == j) Omega[i, j] <- 1; if (abs(i - j) == 1) Omega[i, j] <- 0.4; if (abs(i - j) == 2) Omega[i, j] <- 0.2; if (abs(i - j) == 3) Omega[i, j] <- 0.2 } }
  else if (model_type == "AR4") { for (i in 1:p) for (j in 1:p) { if (i == j) Omega[i, j] <- 1; if (abs(i - j) == 1) Omega[i, j] <- 0.4; if (abs(i - j) == 2) Omega[i, j] <- 0.2; if (abs(i - j) == 3) Omega[i, j] <- 0.2; if (abs(i - j) == 4) Omega[i, j] <- 0.1 } }
  else if (model_type == "Full") { Omega <- matrix(1, p, p); diag(Omega) <- 2 }
  else if (model_type == "Star") { diag(Omega) <- 1; w <- min(0.2, 0.99 / sqrt(p - 1)); for (i in 2:p) { Omega[1, i] <- w; Omega[i, 1] <- w } }
  else if (model_type == "Circle") { for (i in 1:p) for (j in 1:p) { if (i == j) Omega[i, j] <- 1; if (abs(i - j) == 1) Omega[i, j] <- 0.5 }; Omega[1, p] <- 0.4; Omega[p, 1] <- 0.4 }
  return(Omega)
  
}

Abar_second_moment <- function(X) {
  n <- nrow(X)
  A <- crossprod(X) / n            
}



#gradient-sign attack
attack_linf <- function(X, eps, ridge = 1e-4) {
  p <- ncol(X)
  S <- cov(X)
  Omega_hat <- solve(S + ridge * diag(p))
  X + eps * sign(X %*% Omega_hat)
}


Lambda_from_delta <- function(X, delta, penalize_diag = FALSE) {
  # Column-wise scale proxy (mean absolute value per feature)
  mus <- colMeans(abs(X))
  
  # Element-wise penalty: delta * (mus_i * mus_j)
  Lam <- delta * outer(mus, mus) + delta *delta
  
  if (!penalize_diag) diag(Lam) <- 0
  
  Lam
}

# (C) Build a standard scalar l1 penalty matrix Rho from lambda (the usual graphical lasso)
Rho_from_lambda <- function(lambda, p, penalize_diag = FALSE) {
  # Off-diagonals all share the same penalty level "lambda"
  Rho <- matrix(lambda, p, p)
  
  # Typically do not penalize diagonal
  if (!penalize_diag) diag(Rho) <- 0
  
  Rho
}
fit_quic <- function(A, rho) {
  p <- ncol(A)
  
  out <- tryCatch(
    QUIC(A, rho = rho, msg = 0),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  # QUIC typically returns the precision matrix in out$X (sometimes out$Theta)
  if (!is.null(out$X))     return(as.matrix(out$X))
  if (!is.null(out$Theta)) return(as.matrix(out$Theta))
  
  NULL
}
# ======================
# Experiment Settings
# ======================

set.seed(2025)
p <- 20
n_trials <- 100
sample_sizes <- c(20, 40, 60)

test_mult <- 5
val_mult <- 1

models <- c("Heterogeneous","AR1","AR2","AR3","AR4","Full","Star","Circle")

delta_grid   <- seq(0.01, 1, length.out = 30)
lambda_grid  <- seq(0.01, 1, length.out = 30)
epsilon_grid <- seq(0, 1, length.out = 10)

penalize_diag <- FALSE



tuning_results <- data.frame()
robust_results <- data.frame()



cat("Starting experiment...\n")


# ==========
# Main Loop
# ==========

for (model_type in models) {
  cat("Model:", model_type, "\n")
  Omega_true <- generate_precision(model_type, p)
  Sigma_true <- solve(Omega_true) 
  
  
  for (n in sample_sizes) {
    
    n_train <- n
    n_val   <- val_mult * n
    n_test  <- test_mult * n
    
    cat("  n_train =", n_train, "\n")
    for (trial in 1:n_trials) {
      #generate training, validation, and testing data
      X_train <- mvrnorm(n_train, mu = rep(0, p), Sigma = Sigma_true)
      X_val   <- mvrnorm(n_val,   mu = rep(0, p), Sigma = Sigma_true)
      X_test  <- mvrnorm(n_test,  mu = rep(0, p), Sigma = Sigma_true)
      
      A_train <- Abar_second_moment(X_train) 
      
      
      
      # ---- Tune Method 1 (Delta) ----
      
      val_m1 <- rep(NA_real_, length(delta_grid))
      for (i in seq_along(delta_grid)) {
        delta <- delta_grid[i]
        Lam <- Lambda_from_delta(X_train, delta, penalize_diag = penalize_diag)
        C <- fit_quic(A_train, rho = Lam)
        
        #Use Sigma_true for clean validation
        
        if (!is.null(C)) val_m1[i] <- evaluate_kl(C, Sigma_true)
        
      }
      
      if (all(is.na(val_m1))) next
      
      delta_star <- delta_grid[which.min(val_m1)]
      Lam_star <- Lambda_from_delta(X_train, delta_star, penalize_diag = penalize_diag)
      
      #print(delta_star )
      
      C1_star <- fit_quic(A_train, rho = Lam_star)
      
      if (is.null(C1_star)) next
      
      
      
      # ---- Tune Method 2 (Lambda) ----
      
      val_m2 <- rep(NA_real_, length(lambda_grid))
      
      for (j in seq_along(lambda_grid)) {
        lambda <- lambda_grid[j]
        Rho <- Rho_from_lambda(lambda, p, penalize_diag = penalize_diag)
        C <- fit_quic(A_train, rho = Rho)
        # [CHANGE] Use Sigma_true for clean validation
        if (!is.null(C)) val_m2[j] <- evaluate_kl(C, Sigma_true)
        
      }
      
      if (all(is.na(val_m2))) next
      lambda_star <- lambda_grid[which.min(val_m2)]
      Rho_star <- Rho_from_lambda(lambda_star, p, penalize_diag = penalize_diag)
      C2_star <- fit_quic(A_train, rho = Rho_star)
      #print(lambda_star)
      
      if (is.null(C2_star)) next
      tuning_results <- rbind(tuning_results,
                              data.frame(Model=model_type, SampleSize=n, Trial=trial,
                                         delta_star=delta_star, lambda_star=lambda_star))
      
      
      
      # ---- Robustness Evaluation ----

      # 1. Clean Risk: Compare Est against Sigma_true (Theoretical KL)
      r1_clean <- evaluate_kl(C1_star, Sigma_true)
      r2_clean <- evaluate_kl(C2_star, Sigma_true)
      
      
      
      for (eps in epsilon_grid) {
        # 2. Attack Data
        X_test_adv <- attack_linf(X_test, eps = eps, ridge = 1e-4)
        # 3. Adversarial Risk
        # Note: The distribution has physically shifted. The "True" Sigma of this 
        # new data is unknown, so we MUST approximate it with cov(X_test_adv).
        
        Sigma_adv_est <- cov(X_test_adv)
        r1_adv <- evaluate_kl(C1_star, Sigma_adv_est)
        r2_adv <- evaluate_kl(C2_star, Sigma_adv_est)
        robust_results <- rbind(robust_results,
                                data.frame(Model=model_type, SampleSize=n, Trial=trial,
                                           Method="Perturbed", eps=eps,
                                           RiskClean=r1_clean, RiskAdv=r1_adv, riskgap=r1_adv - r1_clean),
                                data.frame(Model=model_type, SampleSize=n, Trial=trial,
                                           Method="Standard", eps=eps,
                                           RiskClean=r2_clean, RiskAdv=r2_adv, riskgap=r2_adv - r2_clean)
        )
      }
    }
  }
}
cat("Done.\n")









# ==============================================================================
# 8) Summaries + Final Visualization
# ==============================================================================
# 1. Summarize Data (Mean & Standard Error)

summ_riskadv <- robust_results %>%
  group_by(Model, SampleSize, eps, Method) %>%
  summarise(
    riskadvMean = mean(riskgap, na.rm=TRUE),
    riskadvSE   = sd(riskgap, na.rm=TRUE) / sqrt(sum(!is.na(riskadv))),
    .groups = "drop"
  ) %>%
  
  # 2. Rename Factors to "Perturbed" and "Standard"
  
  mutate(
    Method_Label = Method,
    
    # Create a nice label for sample sizes (e.g., "n = 50")
    SampleSize_Label = paste0("n = ", SampleSize)
)


# 3. Create a PDF to save all plots
pdf("Robustness_Analysis_Plots.pdf", width = 8, height = 6)



# 4. Loop through models and plot
unique_models <- unique(summ_riskadv$Model)



for (m in unique_models) {
  df_m <- summ_riskadv %>% filter(Model == m)
  p <- ggplot(df_m, aes(x = eps, y = riskadvMean, color = Method_Label, fill = Method_Label)) +
    
    
    
    # --- Geometries ---
    # Confidence Interval Ribbon (Lighter transparency)
    geom_ribbon(aes(ymin = riskadvMean - riskadvSE, ymax = riskadvMean + riskadvSE),
                alpha = 0.2, color = NA) +

    # Main Line
    geom_line(aes(linetype = Method_Label), linewidth = 0.5) +

    # Points (optional, helps see the grid steps)
    geom_point(size = 0.8) +
    
    
    
    # --- Faceting ---
    # Facet by Sample Size (Columns) so we can compare N side-by-side
    facet_wrap(~SampleSize_Label, nrow = 1, scales = "free_y") +
    
    
    # --- Scales & Colors ---
    # Standard (Red, Dashed), Perturbed (Blue, Solid)
    
    scale_color_manual(
      values = c("Standard"  = "#E41A1C", "Perturbed" = "#377EB8"),
      labels = c("Standard" = TeX("$l_1$"), "Perturbed" = "Perturbed")
    ) +
    
    scale_fill_manual(
      values = c("Standard"  = "#E41A1C", "Perturbed" = "#377EB8"),
      labels = c("Standard" = TeX("$l_1$"), "Perturbed" = "Perturbed")
    ) +
    
    scale_linetype_manual(
      values = c("Standard"  = "dashed", "Perturbed" = "solid"),
      labels = c("Standard" = TeX("$l_1$"), "Perturbed" = "Perturbed")
    ) +
    
    
    
    # --- Labels & Theme ---
    labs(
      title = paste0(m, " Model "),
      #subtitle = expression(paste("Performance degradation as perturbation ", epsilon, " increases")),
      x = TeX("Values of $\\epsilon$"),
      y = "KL",
      color = "Method",
      fill = "Method",
      linetype = "Method"
    ) +
    
    
    
    theme_bw(base_size = 11) + # Cleaner white background
    theme(
      plot.title = element_text(hjust = 0.5,size = 14),
      legend.position = "top",
      legend.title = element_text(size = 9),
      strip.background = element_rect(fill = "gray95"), # Light gray header background
      panel.grid.minor = element_blank() # Remove minor grid lines for clarity
    )
  
  
  
  # Print to PDF
  
  print(p)
  file_name <- paste0("Robustness_", m, ".pdf")
  ggsave(filename = file_name, plot = p, width = 6, height = 5)
  
}
