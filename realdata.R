#============================================================
# 0. Setup
#============================================================
pkgs <- c("QUIC","igraph","ggplot2","gridExtra","bnlearn","dplyr","tidyr")
inst <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(inst)) install.packages(inst)

library(QUIC)
library(igraph)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

set.seed(123)

#============================================================
# 1. Load Sachs
#============================================================
# 0. Setup
#============================================================
pkgs <- c("QUIC","igraph","ggplot2","gridExtra","dplyr","tidyr")
inst <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(inst)) install.packages(inst)

library(QUIC)
library(igraph)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

set.seed(123)

#============================================================
# 1. Load Leukemia (Golub)
#   Priority: golubEsets::Golub_Merge (72 x 7129)
#   Fallback: multtest::golub (38 x 3051)
#============================================================
cat(">>> Loading Leukemia (Golub)...\n")

X_raw <- NULL

# --- Try Bioconductor golubEsets first
if (!requireNamespace("golubEsets", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  # If Bioconductor unavailable on your machine, this may fail; then we fallback.
  try(BiocManager::install("golubEsets", update=FALSE, ask=FALSE), silent=TRUE)
}

if (requireNamespace("golubEsets", quietly=TRUE) &&
    requireNamespace("Biobase", quietly=TRUE)) {
  suppressPackageStartupMessages(library(golubEsets))
  suppressPackageStartupMessages(library(Biobase))
  data(Golub_Merge, package="golubEsets")
  X_raw <- exprs(Golub_Merge)   # genes x samples
  X_raw <- t(X_raw)             # samples x genes
  cat(">>> Loaded Golub_Merge: ", nrow(X_raw), "samples x", ncol(X_raw), "genes\n")
}

# --- Fallback to multtest::golub
if (is.null(X_raw)) {
  if (!requireNamespace("multtest", quietly=TRUE)) {
    install.packages("multtest")
  }
  suppressPackageStartupMessages(library(multtest))
  data(golub, package="multtest")
  X_raw <- t(golub)  # multtest::golub is genes x samples; transpose to samples x genes
  cat(">>> Loaded multtest::golub: ", nrow(X_raw), "samples x", ncol(X_raw), "genes\n")
}

n_full <- nrow(X_raw)
p_full <- ncol(X_raw)

#============================================================
# 2. Choose p_use genes (top variance) to make QUIC feasible
#============================================================
p_use <- 200     # change to 100/300/500 as you like
if (p_use >= p_full) p_use <- p_full

vars <- apply(X_raw, 2, var)
keep <- order(vars, decreasing=TRUE)[1:p_use]
X_raw <- X_raw[, keep, drop=FALSE]

# Standardize (important for gene expression)
X_raw <- scale(X_raw)

n_full <- nrow(X_raw)
p_dim  <- ncol(X_raw)
cat(">>> Using: ", n_full, "samples x", p_dim, "genes (top-variance)\n")

#============================================================
# 3. Utilities: QUIC fit, EBIC selection, TNR/MCC, topology plot
#============================================================
safe_logdet <- function(M) {
  ch <- tryCatch(chol(M), error=function(e) NULL)
  if (is.null(ch)) return(NA_real_)
  2*sum(log(diag(ch)))
}

ggm_loglik <- function(Theta, S, n) {
  ld <- safe_logdet(Theta)
  if (is.na(ld)) return(NA_real_)
  (n/2) * (ld - sum(S * Theta))
}

edge_count <- function(adj) sum(adj[upper.tri(adj)])

ebic_score <- function(Theta, adj, S, n, gamma=0.5) {
  p <- ncol(S)
  ll <- ggm_loglik(Theta, S, n)
  if (is.na(ll)) return(NA_real_)
  k <- edge_count(adj)
  -2*ll + k*log(n) + 4*gamma*k*log(p)
}

fit_quic_adj <- function(S_train, rho_mat, thr=1e-5) {
  fit <- try(QUIC(S_train, rho=rho_mat, msg=0), silent=TRUE)
  if (inherits(fit, "try-error") || is.null(fit$X)) return(NULL)
  Theta <- fit$X
  adj <- (abs(Theta) > thr) * 1
  diag(adj) <- 0
  list(Theta=Theta, adj=adj)
}

calc_tnr_mcc <- function(est, true) {
  upper <- upper.tri(est)
  e <- est[upper]; t <- true[upper]
  TP <- sum(e==1 & t==1)
  FP <- sum(e==1 & t==0)
  FN <- sum(e==0 & t==1)
  TN <- sum(e==0 & t==0)
  
  TNR <- TN / (TN + FP + 1e-12)
  denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) + 1e-24)
  MCC <- (TP*TN - FP*FN) / denom
  c(TP=TP, FP=FP, FN=FN, TN=TN, TNR=TNR, MCC=MCC)
}

get_W <- function(delta, svec) {
  P <- (outer(svec, svec, "+") * delta) + delta^2
  diag(P) <- 0
  P
}

select_by_ebic <- function(method=c("Standard","Weighted"),
                           X_sub, tuning_seq, thr=1e-5, gamma=0.5) {
  method <- match.arg(method)
  n <- nrow(X_sub); p <- ncol(X_sub)
  S <- cov(X_sub)
  
  best <- list(score=Inf, tuning=NA_real_, fit=NULL)
  path <- data.frame()
  
  for (tuning in tuning_seq) {
    if (method=="Standard") {
      P <- matrix(tuning, p, p); diag(P) <- 0
    } else {
      svec <- apply(X_sub, 2, sd)
      P <- get_W(tuning, svec)
    }
    
    fit <- fit_quic_adj(S, P, thr=thr)
    if (is.null(fit)) next
    
    sc <- ebic_score(fit$Theta, fit$adj, S, n, gamma=gamma)
    if (is.na(sc) || !is.finite(sc)) next
    
    path <- rbind(path, data.frame(Method=method, Tuning=tuning, EBIC=sc, Edges=edge_count(fit$adj)))
    if (sc < best$score) best <- list(score=sc, tuning=tuning, fit=fit)
  }
  
  if (nrow(path)==0) stop(paste0("No EBIC results for method=", method))
  list(best=best, path=path)
}

plot_net <- function(est, true, main_title) {
  g_est  <- graph_from_adjacency_matrix(est,  mode="undirected", diag=FALSE)
  g_true <- graph_from_adjacency_matrix(true, mode="undirected", diag=FALSE)
  g_u <- graph.union(g_true, g_est)
  lay <- layout_with_kk(g_u)
  
  plot(g_u, layout=lay, vertex.label=NA, vertex.size=5,
       vertex.color="white", vertex.frame.color="black",
       edge.color=NA, main=main_title)
  
  g_fn <- graph.difference(g_true, g_est)
  if(ecount(g_fn)>0) plot(g_fn, layout=lay, add=TRUE, edge.color="gray", edge.lty=2, edge.width=1)
  g_fp <- graph.difference(g_est, g_true)
  if(ecount(g_fp)>0) plot(g_fp, layout=lay, add=TRUE, edge.color="#D55E00", edge.width=2)
  g_tp <- graph.intersection(g_est, g_true, keep.all.vertices=TRUE)
  if(ecount(g_tp)>0) plot(g_tp, layout=lay, add=TRUE, edge.color="#009E73", edge.width=3)
}

#============================================================
# 4. Build a reference graph (silver standard) from FULL data
#   - This is ONLY for evaluation as "true_adj"
#   - sample_sizes still EXCLUDES full data
#============================================================
thr_adj <- 1e-5
gamma_ebic <- 0.5     # EBIC default; set 0 for BIC

grid_len_ref <- 40
lambda_seq_ref <- seq(0, 1, length.out=grid_len_ref)

cat(">>> Building reference graph using FULL data (Standard+EBIC)...\n")
sel_ref <- select_by_ebic("Standard", X_raw, lambda_seq_ref, thr=thr_adj, gamma=gamma_ebic)
ref_adj <- sel_ref$best$fit$adj
cat(">>> Reference edges = ", edge_count(ref_adj), "  (lambda*=", sel_ref$best$tuning, ")\n")

# If you prefer a less biased reference, you can replace ref_adj by:
#   ref_adj <- (ref_adj + sel_wgt_full$best$fit$adj >= 1) * 1   # union
# or intersection, or stability graph, etc.

#============================================================
# 5. Experiment controls (vary sample size, EXCLUDING full data)
#============================================================
n_trials <- 50
grid_len <- 25

sample_sizes <- c(15,20,30,40,50,60,70)
sample_sizes <- sort(unique(sample_sizes))
sample_sizes <- sample_sizes[sample_sizes < n_full]  # exclude full
sample_sizes <- sample_sizes[sample_sizes >= 10]
if (length(sample_sizes)==0) stop("All sample_sizes were removed. Please choose sizes < n_full.")
cat(">>> Using sample sizes (excluding full): ", paste(sample_sizes, collapse=", "), "\n")

#============================================================
# 6. Run varying-n experiment (EBIC tuning)
#============================================================
results <- data.frame()

cat(">>> Running varying-n experiment...\n")
for (n_sub in sample_sizes) {
  cat("  - n =", n_sub, "\n")
  
  for (tt in 1:n_trials) {
    idx <- sample(1:n_full, n_sub, replace=FALSE)
    X_sub <- X_raw[idx, , drop=FALSE]
    S_sub <- cov(X_sub)
    
    lambda_seq <- seq(0, 1, length.out=grid_len)
    delta_seq  <- seq(0, 1, length.out=grid_len)
    
    # --- Standard
    sel_std <- select_by_ebic("Standard", X_sub, lambda_seq, thr=thr_adj, gamma=gamma_ebic)
    fit_std <- sel_std$best$fit
    met_std <- calc_tnr_mcc(fit_std$adj, ref_adj)
    
    results <- rbind(results, data.frame(
      n=n_sub, Trial=tt, Method="Standard",
      Tuning=sel_std$best$tuning, SelScore=sel_std$best$score, Edges=edge_count(fit_std$adj),
      TNR=met_std["TNR"], MCC=met_std["MCC"]
    ))
    
    # --- Weighted
    sel_wgt <- select_by_ebic("Weighted", X_sub, delta_seq, thr=thr_adj, gamma=gamma_ebic)
    fit_wgt <- sel_wgt$best$fit
    met_wgt <- calc_tnr_mcc(fit_wgt$adj, ref_adj)
    
    results <- rbind(results, data.frame(
      n=n_sub, Trial=tt, Method="Weighted",
      Tuning=sel_wgt$best$tuning, SelScore=sel_wgt$best$score, Edges=edge_count(fit_wgt$adj),
      TNR=met_wgt["TNR"], MCC=met_wgt["MCC"]
    ))
  }
}

#============================================================
# 7. Summaries + Plots: TNR/MCC vs n
#============================================================
summ <- results %>%
  group_by(Method, n) %>%
  summarise(
    TNR_mean=mean(TNR, na.rm=TRUE),
    TNR_se=sd(TNR, na.rm=TRUE)/sqrt(sum(!is.na(TNR))),
    MCC_mean=mean(MCC, na.rm=TRUE),
    MCC_se=sd(MCC, na.rm=TRUE)/sqrt(sum(!is.na(MCC))),
    Edges_mean=mean(Edges, na.rm=TRUE),
    Edges_se=sd(Edges, na.rm=TRUE)/sqrt(sum(!is.na(Edges))),
    .groups="drop"
  )

print(summ)

p_tnr <- ggplot(summ, aes(x=n, y=TNR_mean, color=Method, fill=Method)) +
  geom_ribbon(aes(ymin=TNR_mean-TNR_se, ymax=TNR_mean+TNR_se), alpha=0.18, color=NA) +
  geom_line(linewidth=1.0) + geom_point(size=2) +
  theme_bw() +
  labs(title=paste0("TNR vs n (vs Reference Graph, EBIC gamma=", gamma_ebic, ")"),
       x="Sample size n", y="TNR") +
  theme(legend.position="bottom")

p_mcc <- ggplot(summ, aes(x=n, y=MCC_mean, color=Method, fill=Method)) +
  geom_ribbon(aes(ymin=MCC_mean-MCC_se, ymax=MCC_mean+MCC_se), alpha=0.18, color=NA) +
  geom_line(linewidth=1.0) + geom_point(size=2) +
  theme_bw() +
  labs(title=paste0("MCC vs n (vs Reference Graph, EBIC gamma=", gamma_ebic, ")"),
       x="Sample size n", y="MCC") +
  theme(legend.position="bottom")

p_edges <- ggplot(summ, aes(x=n, y=Edges_mean, color=Method, fill=Method)) +
  geom_ribbon(aes(ymin=Edges_mean-Edges_se, ymax=Edges_mean+Edges_se), alpha=0.12, color=NA) +
  geom_line(linewidth=1.0) + geom_point(size=2) +
  theme_bw() +
  labs(title="Average #Edges vs n", x="Sample size n", y="Edges (mean)") +
  theme(legend.position="bottom")

gridExtra::grid.arrange(p_tnr, p_mcc, p_edges, ncol=1)

#============================================================
# 8. Topology plots (largest non-full n)
#============================================================
n_topo <- max(sample_sizes)
cat(">>> Topology plot at n =", n_topo, "(largest non-full sample size)\n")

idx <- sample(1:n_full, n_topo, replace=FALSE)
X_sub <- X_raw[idx, , drop=FALSE]

sel_std <- select_by_ebic("Standard", X_sub, seq(0,1,length.out=grid_len), thr=thr_adj, gamma=gamma_ebic)
sel_wgt <- select_by_ebic("Weighted",  X_sub, seq(0,1,length.out=grid_len), thr=thr_adj, gamma=gamma_ebic)

fit_std <- sel_std$best$fit
fit_wgt <- sel_wgt$best$fit
met_std <- calc_tnr_mcc(fit_std$adj, ref_adj)
met_wgt <- calc_tnr_mcc(fit_wgt$adj, ref_adj)

par(mfrow=c(1,3), mar=c(1,1,3,1))
g_ref <- graph_from_adjacency_matrix(ref_adj, mode="undirected", diag=FALSE)
plot(g_ref, layout=layout_with_kk(g_ref),
     vertex.color="lightgrey", vertex.size=5, vertex.label=NA,
     main=paste0("Reference (full data)\nEdges=", edge_count(ref_adj)))

plot_net(fit_std$adj, ref_adj,
         paste0("Standard (EBIC)\n",
                "n=", n_topo, " edges=", edge_count(fit_std$adj),
                "\nTNR=", round(met_std["TNR"],3),
                " MCC=", round(met_std["MCC"],3)))

plot_net(fit_wgt$adj, ref_adj,
         paste0("Weighted (EBIC)\n",
                "n=", n_topo, " edges=", edge_count(fit_wgt$adj),
                "\nTNR=", round(met_wgt["TNR"],3),
                " MCC=", round(met_wgt["MCC"],3)))
par(mfrow=c(1,1))

cat("✅ Done.\n")
