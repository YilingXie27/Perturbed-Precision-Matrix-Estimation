
# Packages
library("QUIC")
library("MASS")
library(dplyr)
library(ggplot2)
library(ggh4x)
library(latex2exp)
library(foreach)
library(doParallel)

rm(list = ls())


n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterSetRNGStream(cl, 2025)  # reproducible parallel RNG
cat(sprintf("Parallel backend: %d cores\n", n_cores))

# Function to generate 6 Models
generate_matrix <- function(model_type, p) {
  
  # Initialize a zero matrix
  Omega <- matrix(0, p, p)
  
  if (model_type == "AR1") {
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) Omega[i, j] <- 1
        if (abs(i - j) == 1) Omega[i, j] <- 0.5
      }
    }
  }
  
  else if (model_type == "AR2") {
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) Omega[i, j] <- 1
        if (abs(i - j) == 1) Omega[i, j] <- 0.5
        if (abs(i - j) == 2) Omega[i, j] <- 0.25
      }
    }
  }
  
  else if (model_type == "AR3") {
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) Omega[i, j] <- 1
        if (abs(i - j) == 1) Omega[i, j] <- 0.4
        if (abs(i - j) == 2) Omega[i, j] <- 0.2
        if (abs(i - j) == 3) Omega[i, j] <- 0.2
      }
    }
  }
  
  else if (model_type == "AR4") {
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) Omega[i, j] <- 1
        if (abs(i - j) == 1) Omega[i, j] <- 0.4
        if (abs(i - j) == 2) Omega[i, j] <- 0.2
        if (abs(i - j) == 3) Omega[i, j] <- 0.2
        if (abs(i - j) == 4) Omega[i, j] <- 0.1
      }
    }
  }
  
  else if (model_type == "Star") {
    diag(Omega) <- 1
    w <- 0.2
    w <- min(w, 0.99 / sqrt(p - 1)) # promise positive definiteness
    for (i in 2:p) {
      Omega[1, i] <- w
      Omega[i, 1] <- w
    }
  }
  
  else if (model_type == "Circle") {
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) Omega[i, j] <- 1
        if (abs(i - j) == 1) Omega[i, j] <- 0.5
      }
    }
    Omega[1, p] <- 0.4
    Omega[p, 1] <- 0.4
  }
  Omega <- (Omega + t(Omega)) / 2
  eigmin <- min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values)
  if (eigmin <= 1e-8) {
    Omega <- Omega + (abs(eigmin) + 0.1) * diag(p)
  }
  
  return(Omega)
}

# Function to calculate performance metrics
evaluate_performance <- function(true_mat, est_mat, threshold = 1e-5) {
  
  #Frobenius Norm ---
  frob_error <- norm(true_mat - est_mat, type = "F")
  
  #Selection Accuracy ---
  true_edges <- abs(true_mat[upper.tri(true_mat)]) > threshold
  est_edges  <- abs(est_mat[upper.tri(est_mat)]) > threshold
  
  TP <- sum(true_edges == 1 & est_edges == 1)
  FP <- sum(true_edges == 0 & est_edges == 1)
  TN <- sum(true_edges == 0 & est_edges == 0)
  FN <- sum(true_edges == 1 & est_edges == 0)
  
  ACC <- (TP + TN) / (TP + FP + TN + FN)
  
  TNR <- ifelse((TN + FP) == 0, 0, TN / (TN + FP))
  TPR <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- ifelse(denominator == 0, 0, numerator / denominator)
  
  p <- nrow(true_mat)
  true_cov <- tryCatch(solve(true_mat), error = function(e) NULL)
  
  if (is.null(true_cov)) {
    KL_loss <- NA
  } else {
    log_det_est  <- determinant(est_mat, logarithm = TRUE)$modulus[1]
    log_det_true <- determinant(true_mat, logarithm = TRUE)$modulus[1]
    trace_term <- sum(diag(est_mat %*% true_cov))
    KL_loss <- -log_det_est + trace_term - (-log_det_true + p)
  }
  
  return(c(ACC = ACC, TNR = TNR, KL = KL_loss, MCC = MCC, TPR = TPR, Frob = frob_error))
}

# Function to calculate BIC (Data-driven selection metric)
calculate_ebic <- function(S, Omega_Est, n, p, gamma = 0) {
  log_det <- determinant(Omega_Est, logarithm = TRUE)$modulus[1]
  trace_term <- sum(diag(S %*% Omega_Est))
  log_lik <- log_det - trace_term
  
  n_edges <- sum(abs(Omega_Est[upper.tri(Omega_Est)]) > 1e-5)
  
  ebic <- -n * log_lik + n_edges * log(n) + 4 * n_edges * gamma * log(p)
  
  return(ebic)
}

# ==============================================================================
# 2. Simulation Parameters
# ==============================================================================

p <- 30
n_trials <- 100 # repeat 100 times
sample_sizes <- c(10, 20, 30, 40, 50, 60)

# Gird for the parameter values
delta_values <- seq(0.01, 1, length.out = 20)
lambda_values <- seq(0.01, 1, length.out = 20)

models <- c("AR1", "AR2", "AR3", "AR4", "Star", "Circle")

results <- data.frame()

# ==============================================================================
# 3. Main Experiment Loop
# ==============================================================================

set.seed(2025)
cat("Starting Simulation...\n")

results_list <- list()
counter <- 1

# Heterogeneity Factors
scale_factors <- c(rep(10, 5), rep(1, p - 5))
Scaling_Matrix <- diag(scale_factors)
Inverse_Scaling <- diag(1/scale_factors)

for (mod in models) {
  cat(paste("Processing Model:", mod, "\n"))
  
  Base_Omega <- generate_matrix(mod, p)
  Base_Sigma <- solve(Base_Omega)
  
  # Inject Heterogeneity
  True_Sigma <- Scaling_Matrix %*% Base_Sigma %*% Scaling_Matrix
  
  # True Omega = D^-1 * Base_Omega * D^-1
  True_Omega <- Inverse_Scaling %*% Base_Omega %*% Inverse_Scaling
  
  for (n in sample_sizes) {
    print(n)
    datasets <- list()
    for (t in 1:n_trials) {
      raw_X <- mvrnorm(n, mu = rep(0, p), Sigma = True_Sigma)
      datasets[[t]] <- scale(raw_X, center = TRUE, scale = FALSE)
    }
    
    # --- Algorithm 1: Perturbed --- (parallel over trials)
    res_perturbed <- foreach(
      t = 1:n_trials,
      .combine = rbind,
      .packages = c("QUIC"),
      .export = c("evaluate_performance", "calculate_ebic")
    ) %dopar% {
      
      X <- datasets[[t]]
      S <- (t(X) %*% X) / n
      mean_abs_x <- colMeans(abs(X))
      
      best_ebic <- Inf
      best_fit_metrics <- NULL
      best_param_val <- NA
      
      for (d_val in delta_values) {
        Lambda_Mat <- outer(mean_abs_x, mean_abs_x, function(x, y) d_val*x + d_val*y + d_val^2)
        diag(Lambda_Mat) <- 0
        
        fit <- tryCatch({ QUIC(S, rho = Lambda_Mat, msg = 0) }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          current_ebic <- calculate_ebic(S, fit$X, n, p)
          if (current_ebic < best_ebic) {
            best_ebic <- current_ebic
            best_param_val <- d_val
            best_fit_metrics <- evaluate_performance(True_Omega, fit$X)
          }
        }
      }
      
      if (is.null(best_fit_metrics)) return(NULL)
      
      data.frame(
        Model = mod, Algorithm = "Perturbed", SampleSize = as.factor(n),
        Trial = t, Selected_Param = best_param_val,
        Metric = names(best_fit_metrics), Value = as.numeric(best_fit_metrics)
      )
    }
    
    if (!is.null(res_perturbed) && nrow(res_perturbed) > 0) {
      results_list[[counter]] <- res_perturbed
      counter <- counter + 1
    }
    
    #GLasso (raw) --- (parallel over trials)
    res_l1 <- foreach(
      t = 1:n_trials,
      .combine = rbind,
      .packages = c("QUIC"),
      .export = c("evaluate_performance", "calculate_ebic")
    ) %dopar% {
      
      X <- datasets[[t]]
      S <- (t(X) %*% X) / n
      
      best_ebic <- Inf
      best_fit_metrics <- NULL
      best_param_val <- NA
      
      for (l_val in lambda_values) {
        Rho_Standard <- matrix(l_val, p, p)
        diag(Rho_Standard) <- 0
        
        fit <- tryCatch({ QUIC(S, rho = Rho_Standard, msg = 0) }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          current_ebic <- calculate_ebic(S, fit$X, n, p)
          if (current_ebic < best_ebic) {
            best_ebic <- current_ebic
            best_param_val <- l_val
            best_fit_metrics <- evaluate_performance(True_Omega, fit$X)
          }
        }
      }
      
      if (is.null(best_fit_metrics)) return(NULL)
      
      data.frame(
        Model = mod, Algorithm = "ell_1", SampleSize = as.factor(n),
        Trial = t, Selected_Param = best_param_val,
        Metric = names(best_fit_metrics), Value = as.numeric(best_fit_metrics)
      )
    }
    
    if (!is.null(res_l1) && nrow(res_l1) > 0) {
      results_list[[counter]] <- res_l1
      counter <- counter + 1
    }
    
    # --- Algorithm 3 (Baseline): GLasso-STD --- (parallel over trials)
    res_l1_std <- foreach(
      t = 1:n_trials,
      .combine = rbind,
      .packages = c("QUIC"),
      .export = c("evaluate_performance", "calculate_ebic")
    ) %dopar% {
      
      X <- datasets[[t]]          # already centered
      S <- (t(X) %*% X) / n
      
      sdv <- apply(X, 2, sd)
      sdv[sdv == 0 | !is.finite(sdv)] <- 1
      Dinv <- diag(1 / sdv)
      
      X_std <- scale(X, center = FALSE, scale = sdv)
      S_std <- (t(X_std) %*% X_std) / n
      
      best_ebic <- Inf
      best_fit_metrics <- NULL
      best_param_val <- NA
      
      for (l_val in lambda_values) {
        Rho_Standard <- matrix(l_val, p, p)
        diag(Rho_Standard) <- 0
        
        fit <- tryCatch({ QUIC(S_std, rho = Rho_Standard, msg = 0) }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          Omega_std <- fit$X
          Omega_back <- Dinv %*% Omega_std %*% Dinv
          
          current_ebic <- calculate_ebic(S, Omega_back, n, p)
          
          if (current_ebic < best_ebic) {
            best_ebic <- current_ebic
            best_param_val <- l_val
            best_fit_metrics <- evaluate_performance(True_Omega, Omega_back)
          }
        }
      }
      
      if (is.null(best_fit_metrics)) return(NULL)
      
      data.frame(
        Model = mod, Algorithm = "ell_1_std", SampleSize = as.factor(n),
        Trial = t, Selected_Param = best_param_val,
        Metric = names(best_fit_metrics), Value = as.numeric(best_fit_metrics)
      )
    }
    
    if (!is.null(res_l1_std) && nrow(res_l1_std) > 0) {
      results_list[[counter]] <- res_l1_std
      counter <- counter + 1
    }
    
    # --- Algorithm 4 (Baseline): Perturbed-STD --- (parallel over trials)
    # Standardize columns by sample SD, build Lambda using mean_abs on standardized X,
    # fit QUIC on S_std, then map precision back to original scale.
    
    res_perturbed_std <- foreach(
      t = 1:n_trials,
      .combine = rbind,
      .packages = c("QUIC"),
      .export = c("evaluate_performance", "calculate_ebic")
    ) %dopar% {
      
      X <- datasets[[t]]          # already centered
      S <- (t(X) %*% X) / n
      
      # 1) Standardize by sample SD
      sdv <- apply(X, 2, sd)
      sdv[sdv == 0 | !is.finite(sdv)] <- 1
      Dinv <- diag(1 / sdv)
      
      X_std <- scale(X, center = FALSE, scale = sdv)
      S_std <- (t(X_std) %*% X_std) / n
      
      mean_abs_x_std <- colMeans(abs(X_std))
      
      best_ebic <- Inf
      best_fit_metrics <- NULL
      best_param_val <- NA
      
      # 2) Search delta in standardized space, BIC on original scale after mapping back
      for (d_val in delta_values) {
        Lambda_Mat_std <- outer(
          mean_abs_x_std, mean_abs_x_std,
          function(x, y) d_val * x + d_val * y + d_val^2
        )
        diag(Lambda_Mat_std) <- 0
        
        fit <- tryCatch({ QUIC(S_std, rho = Lambda_Mat_std, msg = 0) },
                        error = function(e) NULL)
        
        if (!is.null(fit)) {
          Omega_std  <- fit$X
          Omega_back <- Dinv %*% Omega_std %*% Dinv  # map back to original scale
          
          current_ebic <- calculate_ebic(S, Omega_back, n, p)
          
          if (current_ebic < best_ebic) {
            best_ebic <- current_ebic
            best_param_val <- d_val
            best_fit_metrics <- evaluate_performance(True_Omega, Omega_back)
          }
        }
      }
      
      if (is.null(best_fit_metrics)) return(NULL)
      
      data.frame(
        Model = mod, Algorithm = "Perturbed_std", SampleSize = as.factor(n),
        Trial = t, Selected_Param = best_param_val,
        Metric = names(best_fit_metrics), Value = as.numeric(best_fit_metrics)
      )
    }
    
    if (!is.null(res_perturbed_std) && nrow(res_perturbed_std) > 0) {
      results_list[[counter]] <- res_perturbed_std
      counter <- counter + 1
    }
    
  }
}

# Combine Results
raw_results <- do.call(rbind, results_list)

#stop parallel backend
parallel::stopCluster(cl)


# Summarize (Average over 100 trials)
final_summary <- raw_results %>%
  group_by(Model, Algorithm, SampleSize, Metric) %>%
  summarise(
    Avg_Value = mean(Value),
    SE_Value = sd(Value) / sqrt(n()),
    Avg_Selected_Param = mean(Selected_Param)
  )

print(head(final_summary))


table_mean_se <- final_summary %>%
  filter(Metric %in% c("ACC","MCC","TPR","TNR")) %>%
  filter(Model %in% c("AR3")) %>%
  filter(Algorithm %in% c("Perturbed","ell_1","ell_1_std")) %>%
  filter(SampleSize %in% c(20,30,40)) %>%
  mutate(MeanSE = sprintf("%.3f (%.3f)", Avg_Value, SE_Value)) %>%
  select(Model, Algorithm, SampleSize, Metric, MeanSE) %>%
  pivot_wider(
    names_from = Metric,
    values_from = MeanSE
  ) %>%
  arrange(Model, SampleSize, Algorithm)

print(table_mean_se)