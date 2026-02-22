suppressPackageStartupMessages({
  library(QUIC)
  library(MASS)
  library(plsgenomics)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
})

set.seed(2025)

# --- 1. Settings ---
n_cores <- parallel::detectCores() - 1
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
registerDoParallel(cores = n_cores)

# 100 Repeats of 10-Fold CV
N_REPEATS <- 100
cat(sprintf("Running %d Replications on %d cores...\n", N_REPEATS, n_cores))



get_Perturbed_rho <- function(X, delta) {
  mus <- colMeans(abs(X))
  Rho <- (delta * outer(mus, mus, "+")) + (delta^2)
  diag(Rho) <- 0
  return(Rho)
}

score_lda_prob <- function(Omega, X_tr, y_tr, X_te, class_levels) {
  K <- length(class_levels)
  
  if (is.null(Omega) || all(Omega == 0) || any(!is.finite(Omega))) {
    return(list(
      class = factor(rep(NA_character_, nrow(X_te)), levels = class_levels),
      prob  = matrix(NA_real_, nrow(X_te), K,
                     dimnames = list(NULL, class_levels))
    ))
  }
  
  mus <- lapply(class_levels, function(k) {
    idx <- which(y_tr == k)
    if (length(idx) == 0) return(rep(0, ncol(X_tr)))
    colMeans(X_tr[idx, , drop = FALSE])
  })
  
  class_counts <- table(factor(y_tr, levels = class_levels))
  priors <- log(pmax(class_counts, 1e-10) / length(y_tr))
  
  quads <- sapply(1:K, function(k) {
    -0.5 * as.numeric(t(mus[[k]]) %*% Omega %*% mus[[k]])
  })
  
  XOm <- X_te %*% Omega
  log_scores <- matrix(0, nrow(X_te), K)
  for (k in 1:K) log_scores[, k] <- (XOm %*% mus[[k]]) + quads[k] + priors[k]
  
  probs <- t(apply(log_scores, 1, function(x) {
    if (any(!is.finite(x))) return(rep(NA_real_, K))
    x_safe <- x - max(x)
    ex <- exp(x_safe)
    ex / sum(ex)
  }))
  
  pred_class <- factor(class_levels[max.col(probs, ties.method = "first")],
                       levels = class_levels)
  
  bad <- apply(probs, 1, function(r) all(is.na(r)))
  pred_class[bad] <- NA
  
  list(class = pred_class, prob = probs)
}

make_stratified_folds <- function(y, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- factor(y)
  idx_by_class <- split(seq_along(y), y)
  
  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- integer(0)
  
  for (cls in names(idx_by_class)) {
    idx <- sample(idx_by_class[[cls]])
    parts <- split(idx, rep(1:K, length.out = length(idx)))
    for (k in 1:K) folds[[k]] <- c(folds[[k]], parts[[k]])
  }
  
  folds
}

inner_cv_mcc <- function(pred, truth, class_levels) {
  ok <- !is.na(pred) & !is.na(truth)
  if (sum(ok) < 2) return(0)
  
  pred_ok  <- factor(pred[ok],  levels = class_levels)
  truth_ok <- factor(truth[ok], levels = class_levels)
  
  tbl <- table(pred_ok, truth_ok)
  if (!all(dim(tbl) == c(2, 2))) return(0)
  
  tp <- tbl[1, 1]; tn <- tbl[2, 2]; fp <- tbl[1, 2]; fn <- tbl[2, 1]
  numer <- (tp * tn) - (fp * fn)
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(0)
  return(numer / denom)
}

tune_parameter <- function(X, y, type, class_levels, n_folds = 5, seed = NULL) {
  grid_vals <- seq(0.01, 1.0, length.out = 15)
  if (!is.null(seed)) set.seed(seed)
  
  folds <- make_stratified_folds(y, K = n_folds, seed = seed)
  scores <- numeric(length(grid_vals))
  
  for (i in seq_along(grid_vals)) {
    fold_scores <- rep(NA_real_, n_folds)
    
    for (f in 1:n_folds) {
      idx_val <- folds[[f]]
      idx_tr <- setdiff(seq_len(nrow(X)), idx_val)
      
      if (length(unique(y[idx_tr])) < 2) next
      
      X_inner <- X[idx_tr, , drop = FALSE]
      inner_center <- colMeans(X_inner)
      X_inner_centered <- sweep(X_inner, 2, inner_center, "-")
      S <- crossprod(X_inner_centered) / nrow(X_inner_centered) + diag(1e-6, ncol(X))
      
      X_val_centered <- sweep(X[idx_val, , drop = FALSE], 2, inner_center, "-")
      
      Rho <- if (type == "std") grid_vals[i] else get_Perturbed_rho(X_inner_centered, grid_vals[i])
      
      fit <- try(QUIC(S, rho = Rho, msg = 0), silent = TRUE)
      if (!inherits(fit, "try-error") && !is.null(fit$X)) {
        res <- score_lda_prob(fit$X, X_inner_centered, y[idx_tr], X_val_centered, class_levels)
        # Use MCC for tuning 
        fold_scores[f] <- inner_cv_mcc(res$class, y[idx_val], class_levels)
      }
    }
    
    scores[i] <- mean(fold_scores, na.rm = TRUE)
    if (is.nan(scores[i])) scores[i] <- 0
  }
  
  best_idx <- which.max(scores)
  if (length(best_idx) == 0 || scores[best_idx] == 0) {
    warning("Tuning failed, using default parameter 0.3")
    return(0.3)
  }
  return(grid_vals[best_idx])
}

calc_metrics <- function(pred, prob, truth, class_levels) {
  pred2 <- as.character(pred)
  pred2[is.na(pred2)] <- class_levels[2]  # ALL 当默认（负类）
  pred2 <- factor(pred2, levels = class_levels)
  
  truth2 <- factor(truth, levels = class_levels)
  
  acc <- mean(pred2 == truth2)
  
  tbl <- table(pred2, truth2)
  tp <- tbl[1,1]; tn <- tbl[2,2]; fp <- tbl[1,2]; fn <- tbl[2,1]
  numer <- (tp * tn) - (fp * fn)
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denom == 0) 0 else numer / denom
  tnr <- if ((tn + fp) == 0) 0 else tn / (tn + fp)
  tpr <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  
  c(ACC = acc, MCC = mcc, TNR = tnr, TPR = tpr)
}


# --- 3. Empty result template for robust foreach combining (Fix #3) ---
empty_result <- data.frame(
  Replicate = integer(0), p = integer(0), Fold = integer(0),
  Lambda_Std = numeric(0), Delta_Perturbed = numeric(0),
  Std_ACC = numeric(0), Perturbed_ACC = numeric(0),
  Std_MCC = numeric(0), Perturbed_MCC = numeric(0),
  Std_TNR = numeric(0), Perturbed_TNR = numeric(0),
  Std_TPR = numeric(0), Perturbed_TPR = numeric(0)
)

# --- 4. Main Experiment ---

run_benchmark <- function(X_full, y_full, p_list) {
  # Fix #4: Explicitly set AML as positive class (minority), ALL as negative
  class_levels <- c("AML", "ALL")
  y_full <- factor(y_full, levels = class_levels)
  
  final_res <- foreach(
    rep = 1:N_REPEATS,
    .combine = rbind,
    .packages = c("QUIC", "dplyr"),
    .export = c(
      "get_Perturbed_rho", "score_lda_prob",
      "make_stratified_folds",
      "tune_parameter", "calc_metrics",
      "inner_cv_mcc", "empty_result"
    ),
    # Fix #3: Use "remove" to drop failed workers cleanly
    .errorhandling = "remove"
  ) %dopar% {
    
    out <- tryCatch({
      set.seed(2025 + rep)
      rep_res <- list()
      
      for (p in p_list) {
        folds <- make_stratified_folds(y_full, K = 10, seed = 2025 + rep)
        
        for (f in 1:10) {
          idx_te <- folds[[f]]
          idx_tr <- setdiff(seq_len(nrow(X_full)), idx_te)
          
          train_vars <- apply(X_full[idx_tr, ], 2, var)
          top_idx <- order(train_vars, decreasing = TRUE)[1:p]
          
          # Fix #1: Compute train center AND train sd, apply both to train and test
          X_tr_raw <- X_full[idx_tr, top_idx]
          train_center <- colMeans(X_tr_raw)
          train_sd     <- apply(X_tr_raw, 2, sd)
          train_sd[train_sd == 0] <- 1  # prevent division by zero
          
          X_tr <- scale(X_full[idx_tr, top_idx], center = train_center, scale = train_sd)
          X_te <- scale(X_full[idx_te, top_idx], center = train_center, scale = train_sd)
          
          # Strip scale attributes to avoid downstream surprises
          X_tr <- X_tr[, , drop = FALSE]
          X_te <- X_te[, , drop = FALSE]
          attr(X_tr, "scaled:center") <- NULL
          attr(X_tr, "scaled:scale")  <- NULL
          attr(X_te, "scaled:center") <- NULL
          attr(X_te, "scaled:scale")  <- NULL
          
          y_tr <- y_full[idx_tr]
          y_te <- y_full[idx_te]
          
          if (length(unique(y_tr)) < 2) {
            warning(sprintf("Rep %d, p=%d, Fold %d: Missing class in training", rep, p, f))
            next
          }
          
          S_tr <- crossprod(X_tr) / nrow(X_tr) + diag(1e-6, p)
          
          # 1 GLASSO
          lam <- tune_parameter(
            X_tr, y_tr, "std", class_levels,
            seed = 2025 + rep * 1000 + f
          )
          fit_s <- try(QUIC(S_tr, rho = lam, msg = 0), silent = TRUE)
          Om_s <- if (!inherits(fit_s, "try-error") && !is.null(fit_s$X)) fit_s$X else NULL
          res_s <- score_lda_prob(Om_s, X_tr, y_tr, X_te, class_levels)
          
          # 2. Perturbed
          del <- tune_parameter(
            X_tr, y_tr, "Perturbed", class_levels,
            seed = 2025 + rep * 1000 + f + 500
          )
          Rho_pert <- get_Perturbed_rho(X_tr, del)
          fit_a <- try(QUIC(S_tr, rho = Rho_pert, msg = 0), silent = TRUE)
          Om_a <- if (!inherits(fit_a, "try-error") && !is.null(fit_a$X)) fit_a$X else NULL
          res_a <- score_lda_prob(Om_a, X_tr, y_tr, X_te, class_levels)
          
          m_s <- calc_metrics(res_s$class, res_s$prob, y_te, class_levels)
          m_a <- calc_metrics(res_a$class, res_a$prob, y_te, class_levels)
          
          rep_res[[length(rep_res) + 1]] <- data.frame(
            Replicate = rep,
            p = p,
            Fold = f,
            Lambda_Std = lam,
            Delta_Perturbed = del,
            Std_ACC = m_s["ACC"],
            Perturbed_ACC = m_a["ACC"],
            Std_MCC = m_s["MCC"],
            Perturbed_MCC = m_a["MCC"],
            Std_TNR = m_s["TNR"],
            Perturbed_TNR = m_a["TNR"],
            Std_TPR = m_s["TPR"],
            Perturbed_TPR = m_a["TPR"],
            row.names = NULL
          )
        }
      }
      
      if (length(rep_res) > 0) do.call(rbind, rep_res) else empty_result
    }, error = function(e) {
      warning(sprintf("Replicate %d failed: %s", rep, conditionMessage(e)))
      empty_result
    })
    
    out
  }
  
  final_res
}



library(golubEsets)
data(Golub_Merge)

X_in <- t(exprs(Golub_Merge))
y_in <- factor(pData(Golub_Merge)$ALL.AML)

print(dim(X_in))      # 72 x 7129
print(length(y_in))   # 72
print(table(y_in))

cat("Starting...\n")
cat(sprintf("Data dimensions: n=%d, p_full=%d\n", nrow(X_in), ncol(X_in)))
cat(sprintf("Class distribution: %s\n",
            paste(names(table(y_in)), table(y_in), sep = "=", collapse = ", ")))
cat("Positive class (first level): AML\n")

df_res <- run_benchmark(X_in, y_in, p_list = c(20, 30, 40, 50, 60))

if (is.null(df_res) || nrow(df_res) == 0) {
  stop("No results returned from benchmark")
}
cat(sprintf("\nCompleted: %d total observations\n", nrow(df_res)))

# --- Final Summary with Proper Standard Error ---
cat("\n=== Final Statistical Summary ===\n")

summary_by_rep <- df_res %>%
  group_by(p, Replicate) %>%
  summarise(
    Rep_Std_MCC = mean(Std_MCC, na.rm = TRUE),
    Rep_Perturbed_MCC = mean(Perturbed_MCC, na.rm = TRUE),
    Rep_Std_ACC = mean(Std_ACC, na.rm = TRUE),
    Rep_Perturbed_ACC = mean(Perturbed_ACC, na.rm = TRUE),
    Rep_Std_TNR = mean(Std_TNR, na.rm = TRUE),
    Rep_Perturbed_TNR = mean(Perturbed_TNR, na.rm = TRUE),
    Rep_Std_TPR = mean(Std_TPR, na.rm = TRUE),
    Rep_Perturbed_TPR = mean(Perturbed_TPR, na.rm = TRUE),
    .groups = "drop"
  )

summary_proper <- summary_by_rep %>%
  group_by(p) %>%
  summarise(
    N_Replicates = n(),
    
    Mean_Perturbed_ACC = mean(Rep_Perturbed_ACC, na.rm = TRUE),
    SE_Perturbed_ACC   = sd(Rep_Perturbed_ACC, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Perturbed_ACC))),
    
    Mean_Perturbed_MCC = mean(Rep_Perturbed_MCC, na.rm = TRUE),
    SE_Perturbed_MCC   = sd(Rep_Perturbed_MCC, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Perturbed_MCC))),
    
    Mean_Perturbed_TNR = mean(Rep_Perturbed_TNR, na.rm = TRUE),
    SE_Perturbed_TNR   = sd(Rep_Perturbed_TNR, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Perturbed_TNR))),
    
    Mean_Perturbed_TPR = mean(Rep_Perturbed_TPR, na.rm = TRUE),
    SE_Perturbed_TPR   = sd(Rep_Perturbed_TPR, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Perturbed_TPR))),
    
    Mean_Std_ACC = mean(Rep_Std_ACC, na.rm = TRUE),
    SE_Std_ACC   = sd(Rep_Std_ACC, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Std_ACC))),
    
    Mean_Std_MCC = mean(Rep_Std_MCC, na.rm = TRUE),
    SE_Std_MCC   = sd(Rep_Std_MCC, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Std_MCC))),
    
    Mean_Std_TNR = mean(Rep_Std_TNR, na.rm = TRUE),
    SE_Std_TNR   = sd(Rep_Std_TNR, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Std_TNR))),
    
    Mean_Std_TPR = mean(Rep_Std_TPR, na.rm = TRUE),
    SE_Std_TPR   = sd(Rep_Std_TPR, na.rm = TRUE) / sqrt(sum(!is.na(Rep_Std_TPR))),
    
    .groups = "drop"
  )

cat("\n--- Summary with Proper SE (based on replicate means) ---\n")
print(as.data.frame(summary_proper), digits = 4)

cat("\n=== Benchmark Complete ===\n")

# --- Cleanup parallel backend ---
try(stopImplicitCluster(), silent = TRUE)