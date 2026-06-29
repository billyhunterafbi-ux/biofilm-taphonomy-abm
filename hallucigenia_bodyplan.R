# ============================================================
# HALLUCIGENIA BODY PLAN CONSTRUCTOR
#
# Builds a spatial carcass mask and tissue-specific decay
# rate matrix for use in the taphonomy model (experiment1_v2.R)
#
# TISSUE TYPES AND DECAY RATES:
#
#   1. Sclerotised spines     k = 0.005   (most resistant)
#   2. Cuticularised trunk    k = 0.015
#   3. Lobopod appendages     k = 0.030
#   4. Soft internal tissue   k = 0.060   (least resistant)
#
# BODY PLAN (on 100x100 grid):
#
#   - Trunk: elongated ellipse, slightly curved
#   - Head bulb: small ellipse at anterior end
#   - Spines: 7 pairs, dorsal, angled outward
#   - Lobopods: 7 pairs, ventral, angled downward
#   - Pharynx: soft tissue zone near head
#
# USAGE:
#
#   source("hallucigenia_bodyplan.R")
#   bp <- make_hallucigenia(N = 100)
#
#   # plug into model:
#   C        <- bp$C          # initial carcass concentration
#   k_tissue <- bp$k_tissue   # per-cell decay rate
#
# ============================================================

# ------------------------------------------------------------
# HELPER: fill an ellipse on a matrix
# cx, cy  : centre (row, col)
# rx, ry  : semi-axes (rows, cols)
# angle   : rotation in radians (0 = axis-aligned)
# ------------------------------------------------------------

fill_ellipse <- function(mat, cx, cy, rx, ry, angle = 0, value = 1){
  
  N <- nrow(mat)
  M <- ncol(mat)
  
  for(i in 1:N){
    for(j in 1:M){
      
      # translate
      di <- i - cx
      dj <- j - cy
      
      # rotate
      di_r <-  di * cos(angle) + dj * sin(angle)
      dj_r <- -di * sin(angle) + dj * cos(angle)
      
      if((di_r / rx)^2 + (dj_r / ry)^2 <= 1){
        mat[i, j] <- value
      }
    }
  }
  
  mat
}

# ------------------------------------------------------------
# HELPER: draw a thin rectangle (spine or lobopod)
# cx, cy    : base centre (row, col)
# length    : length in cells
# width     : width in cells (1 or 2 for spines)
# angle     : angle from vertical in radians
# ------------------------------------------------------------

fill_rod <- function(mat, cx, cy, length, width = 1, angle = 0, value = 1){
  
  N <- nrow(mat)
  M <- ncol(mat)
  
  for(k in 0:length){
    for(w in -floor(width/2):floor(width/2)){
      
      # step along rod axis then offset perpendicular
      ri <- round(cx + k * cos(angle) + w * sin(angle))
      rj <- round(cy + k * sin(angle) - w * cos(angle))
      
      if(ri >= 1 && ri <= N && rj >= 1 && rj <= M){
        mat[ri, rj] <- value
      }
    }
  }
  
  mat
}

# ============================================================
# MAIN CONSTRUCTOR
# ============================================================

make_hallucigenia <- function(N = 100){
  
  # decay rate constants by tissue type
  k_spine    <- 0.005
  k_trunk    <- 0.015
  k_lobopod  <- 0.030
  k_soft     <- 0.060
  
  # initialise matrices
  C        <- matrix(0, N, N)   # carcass concentration
  k_tissue <- matrix(0, N, N)   # per-cell decay rate
  tissue   <- matrix(0, N, N)   # tissue type label (for plotting)
  
  # ----------------------------------------------------------
  # TRUNK
  # elongated ellipse, slightly tilted
  # centred at (50, 50), long axis along columns
  # ----------------------------------------------------------
  
  trunk_cx    <- 50
  trunk_cy    <- 50
  trunk_rx    <- 6      # narrow in row direction
  trunk_ry    <- 28     # long in col direction
  trunk_angle <- 0.08   # slight tilt
  
  C        <- fill_ellipse(C,        trunk_cx, trunk_cy, trunk_rx, trunk_ry, trunk_angle, value = 1)
  k_tissue <- fill_ellipse(k_tissue, trunk_cx, trunk_cy, trunk_rx, trunk_ry, trunk_angle, value = k_trunk)
  tissue   <- fill_ellipse(tissue,   trunk_cx, trunk_cy, trunk_rx, trunk_ry, trunk_angle, value = 2)
  
  # ----------------------------------------------------------
  # HEAD BULB
  # small ellipse at anterior (left) end
  # ----------------------------------------------------------
  
  head_cx <- 50
  head_cy <- 22
  
  C        <- fill_ellipse(C,        head_cx, head_cy, 5, 5, value = 1)
  k_tissue <- fill_ellipse(k_tissue, head_cx, head_cy, 5, 5, value = k_trunk)
  tissue   <- fill_ellipse(tissue,   head_cx, head_cy, 5, 5, value = 2)
  
  # ----------------------------------------------------------
  # PHARYNX / SOFT INTERNAL TISSUE
  # small soft zone just posterior to head
  # ----------------------------------------------------------
  
  pharynx_cx <- 50
  pharynx_cy <- 27
  
  C        <- fill_ellipse(C,        pharynx_cx, pharynx_cy, 3, 4, value = 1)
  k_tissue <- fill_ellipse(k_tissue, pharynx_cx, pharynx_cy, 3, 4, value = k_soft)
  tissue   <- fill_ellipse(tissue,   pharynx_cx, pharynx_cy, 3, 4, value = 4)
  
  # ----------------------------------------------------------
  # SPINES
  # 7 pairs, dorsal (above trunk centre row)
  # angled outward at ~40 degrees from vertical
  # evenly spaced along trunk
  # ----------------------------------------------------------
  
  # spine base positions along trunk (col positions)
  spine_cols <- seq(30, 72, length.out = 7)
  
  spine_angle_left  <-  pi * 0.28   # angled up-left
  spine_angle_right <- -pi * 0.28   # angled up-right
  spine_length      <- 9
  spine_width       <- 1
  
  for(sc in spine_cols){
    
    base_row <- trunk_cx - 5   # dorsal surface of trunk
    
    # left spine (dorsal left)
    C        <- fill_rod(C,        base_row, sc, spine_length, spine_width,
                         angle = spine_angle_left,  value = 1)
    k_tissue <- fill_rod(k_tissue, base_row, sc, spine_length, spine_width,
                         angle = spine_angle_left,  value = k_spine)
    tissue   <- fill_rod(tissue,   base_row, sc, spine_length, spine_width,
                         angle = spine_angle_left,  value = 1)
    
    # right spine (dorsal right)
    C        <- fill_rod(C,        base_row, sc, spine_length, spine_width,
                         angle = spine_angle_right, value = 1)
    k_tissue <- fill_rod(k_tissue, base_row, sc, spine_length, spine_width,
                         angle = spine_angle_right, value = k_spine)
    tissue   <- fill_rod(tissue,   base_row, sc, spine_length, spine_width,
                         angle = spine_angle_right, value = 1)
  }
  
  # ----------------------------------------------------------
  # LOBOPODS
  # 7 pairs, ventral (below trunk centre row)
  # shorter and broader than spines
  # angled downward and slightly outward
  # ----------------------------------------------------------
  
  lobopod_cols        <- spine_cols   # same spacing as spines
  lobopod_angle_left  <-  pi * 0.38
  lobopod_angle_right <- -pi * 0.38
  lobopod_length      <- 6
  lobopod_width       <- 2
  
  for(lc in lobopod_cols){
    
    base_row <- trunk_cx + 5   # ventral surface of trunk
    
    # left lobopod
    C        <- fill_rod(C,        base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_left,  value = 1)
    k_tissue <- fill_rod(k_tissue, base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_left,  value = k_lobopod)
    tissue   <- fill_rod(tissue,   base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_left,  value = 3)
    
    # right lobopod
    C        <- fill_rod(C,        base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_right, value = 1)
    k_tissue <- fill_rod(k_tissue, base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_right, value = k_lobopod)
    tissue   <- fill_rod(tissue,   base_row, lc, lobopod_length, lobopod_width,
                         angle = lobopod_angle_right, value = 3)
  }
  
  # ----------------------------------------------------------
  # SOFT INTERNAL TISSUE CORE
  # thin zone running along trunk interior
  # overwrites trunk k values in central strip
  # ----------------------------------------------------------
  
  for(i in (trunk_cx - 2):(trunk_cx + 2)){
    for(j in 28:72){
      if(C[i, j] > 0 && tissue[i, j] == 2){
        k_tissue[i, j] <- k_soft
        tissue[i, j]   <- 4
      }
    }
  }
  
  # ----------------------------------------------------------
  # TISSUE LABEL MATRIX FOR PLOTTING
  # 1 = spine, 2 = trunk cuticle, 3 = lobopod, 4 = soft
  # ----------------------------------------------------------
  
  tissue_labels <- matrix("empty", N, N)
  tissue_labels[tissue == 1] <- "spine"
  tissue_labels[tissue == 2] <- "trunk"
  tissue_labels[tissue == 3] <- "lobopod"
  tissue_labels[tissue == 4] <- "soft"
  
  return(list(
    C             = C,
    k_tissue      = k_tissue,
    tissue_labels = tissue_labels,
    tissue        = tissue
  ))
}

# ============================================================
# DIAGNOSTIC PLOT
# run this to check body plan geometry before plugging
# into the full model
# ============================================================

plot_bodyplan <- function(bp){
  
  library(ggplot2)
  library(reshape2)
  
  df <- melt(bp$tissue)
  colnames(df) <- c("row", "col", "tissue")
  
  df$tissue_name <- factor(
    df$tissue,
    levels = c(0, 1, 2, 3, 4),
    labels = c("empty", "spine", "trunk", "lobopod", "soft")
  )
  
  ggplot(df[df$tissue > 0, ], aes(col, row, fill = tissue_name)) +
    geom_tile() +
    scale_fill_manual(
      values = c(
        "spine"   = "#2c3e6b",
        "trunk"   = "#5b8fa8",
        "lobopod" = "#a8c5a0",
        "soft"    = "#e8a070"
      ),
      name = "Tissue type"
    ) +
    scale_y_reverse() +
    coord_equal() +
    theme_void() +
    labs(title = expression(italic("Hallucigenia")~"body plan — tissue type map"))
}

# ============================================================
# QUICK TEST
# ============================================================

bp <- make_hallucigenia(N = 100)
plot_bodyplan(bp)

cat("\nBody plan summary:\n")
cat("  Cells with carcass material:", sum(bp$C > 0), "\n")
cat("  Spine cells:  ", sum(bp$tissue == 1), "\n")
cat("  Trunk cells:  ", sum(bp$tissue == 2), "\n")
cat("  Lobopod cells:", sum(bp$tissue == 3), "\n")
cat("  Soft cells:   ", sum(bp$tissue == 4), "\n")

