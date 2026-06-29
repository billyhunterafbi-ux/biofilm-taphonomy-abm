# ============================================================
# BURGESS SHALE-STYLE TAPHONOMIC AGENT-BASED MODEL
# EXPERIMENT 1
# Is Lock-in Necessary for Exceptional Preservation?
# ============================================================

# ============================================================
# Setup
# ============================================================

# Use project root instead of absolute paths
base_dir <- getwd()
fig_dir  <- file.path(base_dir, "figures")
out_dir  <- file.path(base_dir, "outputs")

dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

library(ggplot2)
library(gridExtra)
library(reshape2)

source("hallucigenia_bodyplan.R")

set.seed(42)


# ============================================================
# Core functions
# ============================================================

laplacian_neumann_c <- function(Z) {
  N <- nrow(Z); M <- ncol(Z)
  
  Z_up    <- rbind(Z[2, ], Z[-N, ])
  Z_down  <- rbind(Z[-1, ], Z[N-1, ])
  Z_right <- cbind(Z[, -1], Z[, M-1])
  Z_left  <- cbind(Z[, 2], Z[, -M])
  
  Z_up + Z_down + Z_right + Z_left - 4 * Z
}


laplacian_mixed_o2 <- function(Z, bc_top = 1) {
  N <- nrow(Z); M <- ncol(Z)
  
  Z_up    <- rbind(rep(bc_top, M), Z[-N, ])
  Z_down  <- rbind(Z[-1, ], Z[N-1, ])
  Z_right <- cbind(Z[, -1], Z[, M-1])
  Z_left  <- cbind(Z[, 2], Z[, -M])
  
  Z_up + Z_down + Z_right + Z_left - 4 * Z
}


diffuse_o2 <- function(O2, D_O2, bc_top = 1) {
  
  n_sub <- max(1L, ceiling(D_O2 / 0.20))
  D_sub <- D_O2 / n_sub
  
  for (s in seq_len(n_sub)) {
    O2 <- O2 + D_sub * laplacian_mixed_o2(O2, bc_top)
    O2 <- pmax(pmin(O2, 1), 0)
    O2[1, ] <- bc_top
  }
  
  O2
}


anoxia_field <- function(O2, o2_crit = 0.2, sharpness = 20) {
  1 / (1 + exp(sharpness * (O2 - o2_crit)))
}


compute_fpi <- function(C_final, anox_final, body_mask, C_initial) {
  mean(C_final[body_mask] * anox_final[body_mask]) /
    (mean(C_initial[body_mask]) + 1e-9)
}


# ============================================================
# Body plan
# ============================================================

make_hallucigenia_offset <- function(N = 100, target_row = 10) {
  
  bp    <- make_hallucigenia(N)
  shift <- target_row - 50
  
  shift_matrix <- function(mat, fill = 0) {
    result <- matrix(fill, N, ncol(mat))
    
    for (r in seq_len(N)) {
      r_new <- r + shift
      if (r_new >= 1 && r_new <= N) result[r_new, ] <- mat[r, ]
    }
    
    result
  }
  
  bp$C             <- shift_matrix(bp$C, fill = 0)
  bp$k_tissue      <- shift_matrix(bp$k_tissue, fill = 0)
  bp$tissue        <- shift_matrix(bp$tissue, fill = 0)
  bp$tissue_labels <- shift_matrix(bp$tissue_labels, fill = "empty")
  
  bp
}


# ============================================================
# Model
# ============================================================

run_model <- function(
    k_tissue,
    body_mask,
    tissue_ref,
    k_bio        = 0.03,
    D_C          = 0.01,
    D_O2         = 0.10,
    k_resp_aero  = 0.015,
    k_resp_anaer = 0.005,
    k_sed        = 0.002,
    o2_crit      = 0.2
) {
  
  N     <- 100
  steps <- 120
  
  bp        <- make_hallucigenia_offset(N, target_row = 10)
  C         <- bp$C
  C_initial <- bp$C
  
  O2  <- matrix(1, N, N)
  phi <- matrix(runif(N * N, 0.2, 1), N, N)
  B   <- matrix(0, N, N)
  
  O2[1, ] <- 1
  
  mean_anox       <- numeric(steps)
  tissue_survival <- matrix(0, steps, 4)
  
  for (t in seq_len(steps)) {
    
    anox <- anoxia_field(O2, o2_crit)
    
    # Transport
    C <- C + D_C * laplacian_neumann_c(C)
    C <- pmax(pmin(C, 1), 0)
    
    # Diffusion
    O2 <- diffuse_o2(O2, D_O2)
    
    # Biofilm growth
    B <- B + anox * k_bio * C * phi * (1 - B)
    B <- pmax(pmin(B, 1), 0)
    
    # Tissue decay
    C <- C - (1 - anox) * k_tissue * C * (1 - 0.7 * B)
    C <- pmax(pmin(C, 1), 0)
    
    # Oxygen consumption
    O2 <- O2 -
      k_resp_aero * C * (1 - anox) -
      k_resp_anaer * B * anox -
      k_sed
    
    O2 <- pmax(pmin(O2, 1), 0)
    O2[1, ] <- 1
    
    # Diagnostics
    mean_anox[t] <- mean(anox[body_mask])
    
    for (k in 1:4) {
      tissue_survival[t, k] <- mean(C[tissue_ref == k])
    }
  }
  
  anox_final <- anoxia_field(O2, o2_crit)
  FPI        <- compute_fpi(C, anox_final, body_mask, C_initial)
  
  list(
    FPI             = FPI,
    C_final         = C,
    O2_final        = O2,
    anox_final      = anox_final,
    mean_anox       = mean_anox,
    tissue_survival = tissue_survival
  )
}


# ============================================================
# Experimental conditions
# ============================================================

conditions <- data.frame(
  condition   = c("Fast arrest", "Intermediate", "Slow arrest", "Oxic flush"),
  D_O2        = c(0.05, 0.10, 0.20, 0.50),
  k_resp_aero = c(0.100, 0.075, 0.050, 0.025)
)

n_runs        <- 20
FPI_threshold <- 0.10

results  <- data.frame()
rep_sims <- list()


# ============================================================
# Run experiment
# ============================================================

bp_base    <- make_hallucigenia_offset(100, 10)
k_tissue   <- bp_base$k_tissue
body_mask  <- bp_base$C > 0
tissue_ref <- bp_base$tissue

for (i in seq_len(nrow(conditions))) {
  
  cond_name  <- conditions$condition[i]
  D_O2_val   <- conditions$D_O2[i]
  k_resp_val <- conditions$k_resp_aero[i]
  
  FPI_values <- numeric(n_runs)
  
  for (j in seq_len(n_runs)) {
    
    set.seed(1000 + i * 100 + j)
    
    sim <- run_model(
      k_tissue    = k_tissue,
      body_mask   = body_mask,
      tissue_ref  = tissue_ref,
      D_O2        = D_O2_val,
      k_resp_aero = k_resp_val
    )
    
    FPI_values[j] <- sim$FPI
    
    if (j == 1) rep_sims[[cond_name]] <- sim
  }
  
  results <- rbind(results,
                   data.frame(condition = cond_name, FPI = FPI_values))
}

results$condition <- factor(results$condition,
                            levels = conditions$condition)


# ============================================================
# Summary
# ============================================================

summary_df <- aggregate(FPI ~ condition, results, mean)
summary_df$sd <- aggregate(FPI ~ condition, results, sd)$FPI
summary_df$probability <- aggregate(
  FPI ~ condition, results,
  function(x) mean(x > FPI_threshold)
)$FPI


# ============================================================
# Plot example
# ============================================================

p1 <- ggplot(results, aes(condition, FPI, fill = condition)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "none") +
  labs(
    title = "FPI by oxygen resupply regime (v10)",
    y = "Normalised FPI"
  )

ggsave(file.path(fig_dir, "exp1v10_fpi.png"),
       p1, width = 7, height = 5)
