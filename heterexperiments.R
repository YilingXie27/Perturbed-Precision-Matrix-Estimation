# ============================================================
# Simulation: heterogeneous marginal scales advantage demo
# Methods: (1) glasso-raw (scalar rho), (2) glasso-std (scalar rho on standardized data),
#          (3) weighted (your adaptive rho matrix)
# Selection: EBIC over a grid (same EBIC form for all)
# Targeted diagnostic: FP_H (false edges touching any high-variance node)
# ============================================================

suppressPackageStartupMessages({
  library(QUIC)
  library(MASS)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------
# 0. Utilities
# ----------------------------

ar1_corr <- function(p, rho = 0.4) {
  idx <- 1:p
  outer(idx, idx, function(i, j) rho^abs(i - j))
}

make_block_sigma <- function(dL = 60, dH = 20, M = 4, rho = 0.4) {
  Rl <- ar1_corr(dL, rho)
  Dl <- diag(rep(1, dL))
  Dh <- diag(rep(M, dH))
  SigL <- Dl %*% Rl %*% Dl
  SigH <- Dh %*% diag(dH) %*% Dh
  Sigma <- Matrix::bdiag(SigL, SigH)
  as.matrix(Sigma)
}

sample_cov_1n <- function(X) {
  # 1/n covariance (matches your bar{A} notation)
  n <- nrow(X)
  crossprod(X) / n
}

logdet_spd <- function(Omega, tol = 1e-10) {
  # robust logdet for SPD matrix using Cholesky
  # return -Inf / error as NA; caller handles
  Omega <- (Omega + t(Omega)) / 2
  # quick eigenvalue guard (optional but helps)
  # use chol directly (fast, fails if not SPD)
  R <- chol(Omega, pivot = FALSE)
  2 * sum(log(diag(R) + tol))
}

ebic_score <- function(Omega, S, n, gamma = 0.5, tol_edge = 1e-6) {
  # EBIC: n * (tr(S Omega) - logdet(Omega)) + df log n + 4 gamma df log d
  # S is 1/n covariance; Omega should be SPD
  d <- ncol(S)
  Omega <- (Omega + t(Omega)) / 2
  
  ld <- tryCatch(logdet_spd(Omega), error = function(e) NA_real_)
  if (!is.finite(ld)) return(Inf)
  
  nll <- sum(S * Omega) - ld
  
  off <- abs(Omega); diag(off) <- 0
  k_edges <- sum(off[upper.tri(off)] > tol_edge)
  df <- d + k_edges
  
  n * nll + df * log(n) + 4 * gamma * df * log(d)
}

safe_quic <- function(S, rho) {
  # rho can be scalar or matrix
  out <- tryCatch(
    QUIC(S, rho = rho, msg = 0),
    error = function(e) NULL
  )
  out
}

fit_quic_scalar <- function(S, n, lambda_grid, gamma = 0.5, tol_edge = 1e-6) {
  best_score <- Inf
  best_Omega <- NULL
  best_lambda <- NA_real_
  
  for (lam in lambda_grid) {
    fit <- safe_quic(S, rho = lam)
    if (is.null(fit)) next
    Omega <- fit$X
    sc <- ebic_score(Omega, S, n = n, gamma = gamma, tol_edge = tol_edge)
    if (is.finite(sc) && sc < best_score) {
      best_score <- sc
      best_Omega <- Omega
      best_lambda <- lam
    }
  }
  
  list(score = best_score, Omega = best_Omega, lambda = best_lambda)
}

fit_quic_weighted <- function(S, X, delta_grid, gamma = 0.5, tol_edge = 1e-6) {
  n <- nrow(X)
  d <- ncol(X)
  
  omega_hat <- colMeans(abs(X))  # \hat{\omega}_i
  
  best_score <- Inf
  best_Omega <- NULL
  best_delta <- NA_real_
  
  for (delta in delta_grid) {
    Lam <- matrix(0, d, d)
    # vectorized construction
    # Lam_ij = 2 delta (omega_i + omega_j) + 2 delta^2, i != j
    Om <- matrix(rep(omega_hat, each = d), d, d)
    Lam <- 2 * delta * (Om + t(Om)) + 2 * delta^2
    diag(Lam) <- 0
    Lam <- (Lam + t(Lam)) / 2
    
    fit <- safe_quic(S, rho = Lam)
    if (is.null(fit)) next
    Omega <- fit$X
    sc <- ebic_score(Omega, S, n = n, gamma = gamma, tol_edge = tol_edge)
    if (is.finite(sc) && sc < best_score) {
      best_score <- sc
      best_Omega <- Omega
      best_delta <- delta
    }
  }
  
  list(score = best_score, Omega = best_Omega, delta = best_delta)
}

support_metrics <- function(Omega_hat, Omega_true, idxH, tol_edge = 1e-6) {
  # returns TPR/FPR/MCC + FP_H (edges touching high-variance nodes)
  d <- ncol(Omega_true)
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2
  Omega_true <- (Omega_true + t(Omega_true)) / 2
  
  Ahat <- (abs(Omega_hat) > tol_edge); diag(Ahat) <- FALSE
  Atru <- (abs(Omega_true) > tol_edge); diag(Atru) <- FALSE
  
  ut <- upper.tri(Atru)
  ah <- Ahat[ut]; at <- Atru[ut]
  
  TP <- sum(ah & at)
  FP <- sum(ah & !at)
  FN <- sum(!ah & at)
  TN <- sum(!ah & !at)
  
  # FP involving any high-variance node
  ij <- which(ut, arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]
  touchH <- (i %in% idxH) | (j %in% idxH)
  FP_H <- sum(ah & !at & touchH)
  
  # Optional decomposition (nice for diagnostics)
  HH <- (i %in% idxH) & (j %in% idxH)
  HL <- touchH & !HH
  FP_HH <- sum(ah & !at & HH)
  FP_HL <- sum(ah & !at & HL)
  
  # Rates
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), NA_real_)
  
  # MCC (avoid integer overflow by forcing numeric)
  TP <- as.numeric(TP); FP <- as.numeric(FP); FN <- as.numeric(FN); TN <- as.numeric(TN)
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- ifelse(is.finite(den) & den > 0, (TP * TN - FP * FN) / den, NA_real_)
  
  list(TP = TP, FP = FP, FN = FN, TN = TN,
       TPR = TPR, FPR = FPR, MCC = MCC,
       FP_H = FP_H, FP_HH = FP_HH, FP_HL = FP_HL)
}

# ----------------------------
# 1. One replicate for a given M
# ----------------------------
run_one <- function(n = 200, dL = 60, dH = 20, M = 4, rho = 0.4,
                    lambda_grid, delta_grid,
                    gamma = 0.5, tol_edge = 1e-6, diag_jitter = 1e-8) {
  
  d <- dL + dH
  idxH <- (dL + 1):d
  
  Sigma <- make_block_sigma(dL, dH, M, rho)
  Omega_true <- solve(Sigma)
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = Sigma)
  S <- sample_cov_1n(X)
  # jitter to avoid numerical issues (usually unnecessary when n>d, but safe)
  S <- S + diag(diag_jitter, d)
  
  # (1) glasso-raw
  fit_raw <- fit_quic_scalar(S, n = n, lambda_grid = lambda_grid, gamma = gamma, tol_edge = tol_edge)
  if (is.null(fit_raw$Omega)) stop("glasso-raw failed on all lambdas; enlarge grid or reduce gamma.")
  met_raw <- support_metrics(fit_raw$Omega, Omega_true, idxH, tol_edge = tol_edge)
  
  # (2) glasso-std: standardize columns by sample sd, fit on standardized S, then map back
  sds <- apply(X, 2, sd)
  Xs <- scale(X, center = FALSE, scale = sds)
  Ss <- sample_cov_1n(Xs) + diag(diag_jitter, d)
  fit_std <- fit_quic_scalar(Ss, n = n, lambda_grid = lambda_grid, gamma = gamma, tol_edge = tol_edge)
  if (is.null(fit_std$Omega)) stop("glasso-std failed on all lambdas; enlarge grid or reduce gamma.")
  Dhat <- diag(sds)
  Omega_std_back <- solve(Dhat) %*% fit_std$Omega %*% solve(Dhat)
  met_std <- support_metrics(Omega_std_back, Omega_true, idxH, tol_edge = tol_edge)
  
  # (3) weighted (your method) on raw scale
  fit_w <- fit_quic_weighted(S, X, delta_grid = delta_grid, gamma = gamma, tol_edge = tol_edge)
  if (is.null(fit_w$Omega)) stop("weighted failed on all deltas; enlarge grid or reduce gamma.")
  met_w <- support_metrics(fit_w$Omega, Omega_true, idxH, tol_edge = tol_edge)
  
  data.frame(
    M = M,
    method = c("glasso-raw", "glasso-std", "weighted"),
    FP = c(met_raw$FP, met_std$FP, met_w$FP),
    FP_H = c(met_raw$FP_H, met_std$FP_H, met_w$FP_H),
    FP_HH = c(met_raw$FP_HH, met_std$FP_HH, met_w$FP_HH),
    FP_HL = c(met_raw$FP_HL, met_std$FP_HL, met_w$FP_HL),
    TPR = c(met_raw$TPR, met_std$TPR, met_w$TPR),
    FPR = c(met_raw$FPR, met_std$FPR, met_w$FPR),
    MCC = c(met_raw$MCC, met_std$MCC, met_w$MCC),
    lambda_raw = c(fit_raw$lambda, NA, NA),
    lambda_std = c(NA, fit_std$lambda, NA),
    delta = c(NA, NA, fit_w$delta)
  )
}

# ----------------------------
# 2. Full simulation over M values
# ----------------------------
run_simulation <- function(Ms = c(1, 2, 4, 8),
                           reps = 30,
                           n = 200, dL = 60, dH = 20, rho = 0.4,
                           gamma = 0.5, tol_edge = 1e-6,
                           lambda_grid = 10^seq(-5, -0.2, length.out = 40),
                           delta_grid  = 10^seq(-5, -0.2, length.out = 40),
                           seed = 123) {
  
  set.seed(seed)
  out <- vector("list", length(Ms) * reps)
  k <- 1L
  
  for (M in Ms) {
    for (r in 1:reps) {
      df <- run_one(n = n, dL = dL, dH = dH, M = M, rho = rho,
                    lambda_grid = lambda_grid, delta_grid = delta_grid,
                    gamma = gamma, tol_edge = tol_edge)
      df$rep <- r
      out[[k]] <- df
      k <- k + 1L
      if (r %% 5 == 0) cat(sprintf("M=%s: finished %d/%d reps\n", M, r, reps))
    }
  }
  
  all_res <- bind_rows(out)
  all_res
}

# ----------------------------
# 3. Run + summarize + plot
# ----------------------------

Ms <- c(1, 2, 4, 8)
reps <- 30
n <- 200
dL <- 60
dH <- 20
rho <- 0.4
gamma <- 0.5      # EBIC gamma; if too sparse, try 0.25
tol_edge <- 1e-6

lambda_grid <- 10^seq(-5, -0.2, length.out = 40)
delta_grid  <- 10^seq(-5, -0.2, length.out = 40)

all_res <- run_simulation(Ms = Ms, reps = reps, n = n, dL = dL, dH = dH, rho = rho,
                          gamma = gamma, tol_edge = tol_edge,
                          lambda_grid = lambda_grid, delta_grid = delta_grid,
                          seed = 2026)

# Targeted normalization for FP_H
NH <- dH * dL + choose(dH, 2)
all_res <- all_res %>% mutate(FPR_H = FP_H / NH)

# --- Quick sanity checks ---
cat("\nCounts:\n")
print(with(all_res, table(method, M)))

# --- Aggregate table (mean + SE + 95% CI) ---
summ <- all_res %>%
  group_by(method, M) %>%
  summarise(
    FP_mean = mean(FP), FP_se = sd(FP)/sqrt(n()),
    FP_H_mean = mean(FP_H), FP_H_se = sd(FP_H)/sqrt(n()),
    FPR_H_mean = mean(FPR_H), FPR_H_se = sd(FPR_H)/sqrt(n()),
    TPR_mean = mean(TPR), TPR_se = sd(TPR)/sqrt(n()),
    FPR_mean = mean(FPR), FPR_se = sd(FPR)/sqrt(n()),
    MCC_mean = mean(MCC, na.rm = TRUE), MCC_se = sd(MCC, na.rm = TRUE)/sqrt(sum(!is.na(MCC))),
    .groups = "drop"
  ) %>%
  mutate(
    FP_H_lo = FP_H_mean - 1.96 * FP_H_se,
    FP_H_hi = FP_H_mean + 1.96 * FP_H_se,
    FPR_H_lo = FPR_H_mean - 1.96 * FPR_H_se,
    FPR_H_hi = FPR_H_mean + 1.96 * FPR_H_se
  )

cat("\nSummary (means):\n")
print(summ %>% arrange(M, method))

# --- Plot 1: FP_H vs M ---
p1 <- ggplot(summ, aes(x = M, y = FP_H_mean, group = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = FP_H_lo, ymax = FP_H_hi), width = 0.15, linewidth = 0.3) +
  scale_x_continuous(breaks = sort(unique(summ$M))) +
  labs(x = "Scale ratio M (high/low variance)", y = "FP_H (false edges touching high-variance nodes)")

# --- Plot 2: FPR_H vs M (normalized) ---
p2 <- ggplot(summ, aes(x = M, y = FPR_H_mean, group = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = FPR_H_lo, ymax = FPR_H_hi), width = 0.15, linewidth = 0.3) +
  scale_x_continuous(breaks = sort(unique(summ$M))) +
  labs(x = "Scale ratio M (high/low variance)", y = "FPR_H = FP_H / (dH*dL + choose(dH,2))")

print(p1)
print(p2)

# Optional: save results
# write.csv(all_res, "all_res_hetero_scales.csv", row.names = FALSE)
# ggsave("FPH_vs_M.pdf", p1, width = 6, height = 4)
# ggsave("FPRH_vs_M.pdf", p2, width = 6, height = 4)