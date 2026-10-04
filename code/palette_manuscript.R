# =============================================================================
# palette_manuscript.R — locked final colour and shape conventions
#
# Single source of truth for group / cohort / feature / marker / direction
# encodings used by all figure scripts. No other script should define its own
# group colours.
#
#   1 Group:     APTB == PTB / HC / LTBI / PNEU / LC / EPTB
#   2 Cohort:    IDC / ITC / PVC — grayscale + solid / dashed / dotted
#   3 Feature:   SOD2 / GOLPH3 / PCNX2 — ink #4D4D4D + solid shapes
#   4 Marker:    circRNA TB-index star / Sweeney3 open circle / RISK6 open diamond
#   5 Direction: Up / Down
#
# Usage:  source("code/palette_manuscript.R")
# =============================================================================
# ---- 1. Group ----
GROUP_COLS <- c(
  APTB = "#FF4B4B",
  PTB  = "#FF4B4B",  # ≡ APTB
  HC   = "#4DAF4A",
  LTBI = "#A259D6",
  PNEU = "#3B82F6",
  LC   = "#FF8C3B",
  EPTB = "#00ACC1"
)

IDC_GROUP_COLS <- GROUP_COLS[c("APTB", "HC", "LTBI", "PNEU", "LC")]
ITC_GROUP_COLS <- GROUP_COLS[c("APTB", "HC", "LC")]
PVC_GROUP_COLS <- GROUP_COLS[c("PTB", "EPTB", "HC")]

# ---- 2. Cohort: line identity, not hue ----
COHORT_COLS <- c(
  IDC = "#2B2B2B",
  ITC = "#6E6E6E",
  PVC = "#A0A0A0"
)
COHORT_LTY <- c(
  IDC = "solid",
  ITC = "dashed",
  PVC = "dotted"
)
COHORT_LWD <- c(
  IDC = 1.20,
  ITC = 1.05,
  PVC = 1.05
)

# ---- 3. Feature / genes ----
# Display / lock order matches manuscript: SOD2, GOLPH3, PCNX2.
# ggplot2 pch has no hexagon; GOLPH3 uses filled diamond (18) as stand-in.

FEATURE_COL <- "#4D4D4D"
GENE_COL <- FEATURE_COL
PANEL_GENES <- c("SOD2", "GOLPH3", "PCNX2")
GENE_COLS <- c(
  SOD2   = FEATURE_COL,
  GOLPH3 = FEATURE_COL,
  PCNX2  = FEATURE_COL
)
GENE_SHAPE <- c(
  SOD2   = 15,  # filled square
  GOLPH3 = 18,  # filled diamond ≈ hexagon stand-in in geom_point
  PCNX2  = 17   # filled triangle
)
GENE_GLYPH <- c(
  SOD2   = "\u25A0",  # black square
  GOLPH3 = "\u2B22",  # black hexagon (swatches / text legends)
  PCNX2  = "\u25B2"   # black triangle
)
GENE_LTY <- c(
  SOD2   = "solid",
  GOLPH3 = "solid",
  PCNX2  = "solid"
)
GENE_PATTERN <- c(
  SOD2   = "solid",
  GOLPH3 = "solid",
  PCNX2  = "solid"
)

# ---- 4. Marker / signatures ----
# ★ circRNA TB-index · ○ Sweeney3 · ◇ RISK6
MARKER_COL <- "#4D4D4D"
SIGNATURE_COL <- MARKER_COL
SIGNATURE_COLS <- c(
  "circRNA TB-index" = MARKER_COL,
  "Sweeney3"         = MARKER_COL,
  "RISK6"            = MARKER_COL
)
SIGNATURE_SHAPE <- c(
  "circRNA TB-index" = 8,  # star
  "Sweeney3"         = 1,  # open circle
  "RISK6"            = 5   # open diamond
)
SIGNATURE_LTY <- c(
  "circRNA TB-index" = "solid",
  "Sweeney3"         = "dashed",
  "RISK6"            = "dotted"
)
SIGNATURE_PATTERN <- c(
  "circRNA TB-index" = "solid",
  "Sweeney3"         = "stripe",
  "RISK6"            = "crosshatch"
)

# Back-compat alias
OTHER_COL <- FEATURE_COL

# ---- 5. Direction ----
DIR_COLS <- c(
  Up   = "#E64B35",
  Down = "#3C5488"
)

# ---- misc (functional / fallback; not clinical group hues) ----
FUNC_COLS <- c("#4D4D4D", "#6E6E6E", "#A0A0A0", "#C8C8C8",
               unname(GROUP_COLS[c("APTB", "HC", "LTBI", "PNEU")]))

roc_cols_idc <- function() {
  c(
    APTB_vs_HC    = unname(GROUP_COLS[["HC"]]),
    APTB_vs_LTBI  = unname(GROUP_COLS[["LTBI"]]),
    APTB_vs_PNEU  = unname(GROUP_COLS[["PNEU"]]),
    APTB_vs_LC    = unname(GROUP_COLS[["LC"]]),
    APTB_vs_nAPTB = unname(GROUP_COLS[["APTB"]])
  )
}
roc_cols_itc <- function() {
  c(
    APTB_vs_HC = unname(GROUP_COLS[["HC"]]),
    APTB_vs_LC = unname(GROUP_COLS[["LC"]])
  )
}
roc_cols_pvc <- function() {
  c(
    PTB_vs_HC  = unname(GROUP_COLS[["HC"]]),
    EPTB_vs_HC = unname(GROUP_COLS[["EPTB"]])
  )
}

VENN_CONTRAST_COLS <- c(
  "vs HC"   = unname(GROUP_COLS[["HC"]]),
  "vs LTBI" = unname(GROUP_COLS[["LTBI"]]),
  "vs PNEU" = unname(GROUP_COLS[["PNEU"]]),
  "vs LC"   = unname(GROUP_COLS[["LC"]])
)
