# Pesticide ecological risk screening tool (R Shiny app)
# Screens pesticide risk to waterways, pollinators and soil fauna at a chosen site.
# A. D. Canning, S. Lewis and D. J. Irvine. Version 1 (18 November 2026).
# Run with shiny::runApp(). Data files and Earth Engine key: see README.md.

# ---- Packages ----
library(shiny)
library(dplyr)
library(ggplot2)
LEGEND_BIG <- theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15, face = "bold"))
library(forcats)
library(ggpubr)
library(DT)
library(ggrepel)
library(jsonlite)
library(shinyjs)
library(leaflet)
HAS_DRAW <- requireNamespace("leaflet.extras", quietly = TRUE)
library(httr)
library(jose)
library(grid)
library(tidyr)
library(gridExtra)
library(patchwork)

# ---- Version and data ----
APP_VERSION <- "1 (18 November 2026)"

`%||%` <- function(x, y) if (is.null(x)) y else x

raw_data <- read.csv("Pesticides_APVMA.csv", header = TRUE, encoding = "UTF-8") %>%
  filter(MaximumPermittedIngredientApplicationRateumolPAIha > 0)
raw_data$pcode_int <- suppressWarnings(as.integer(trimws(as.character(raw_data$pcode))))
ALL_ACTIVES <- sort(unique(raw_data$PesticideactiveIngredient))
uses_raw <- tryCatch(if (file.exists("Pesticides_APVMA_uses.rds")) readRDS("Pesticides_APVMA_uses.rds") else NULL,
                     error = function(e) NULL)
name_key <- function(x) {
  w <- strsplit(gsub("[^A-Z0-9 ]", "", toupper(as.character(x))), " +")
  sing <- function(z) { z <- sub("OES$", "O", z); z <- sub("IES$", "Y", z); z <- sub("(CH|SH|X)ES$", "\\1", z)
                        sub("([^S])S$", "\\1", z) }
  vapply(w, function(z) paste(sing(z[nzchar(z)]), collapse = ""), "")
}
crop_base <- function(x) trimws(sub("\\s+-\\s.*$|\\s*\\(.*$|\\s*[:;,/].*$|\\s+-$", "", as.character(x)))
clean_name <- function(x) trimws(gsub("[?]", "'", gsub("[,:;?]+$", "", as.character(x))))
name_groups <- function(f, base = FALSE, w = NULL) {
  lv <- levels(f); nm <- if (base) crop_base(lv) else clean_name(lv); k <- name_key(nm)
  k[!nzchar(k)] <- NA
  w <- if (is.null(w)) tabulate(as.integer(f), length(lv)) else w
  ok <- !is.na(k)
  lab <- tapply(seq_along(lv)[ok], k[ok], function(i) { t <- tapply(w[i], nm[i], sum); names(t)[which.max(t)] })
  list(grp = factor(k[as.integer(f)]), label = lab)
}
HOST_LABEL <- character(0); PEST_LABEL <- character(0)
if (!is.null(uses_raw)) {
  g <- name_groups(uses_raw$host, base = TRUE); uses_raw$host_grp <- g$grp; HOST_LABEL <- g$label
  g <- name_groups(uses_raw$pest);              uses_raw$pest_grp <- g$grp; PEST_LABEL <- g$label
  rm(g)
}
host_label <- function(k) { x <- HOST_LABEL[k]; ifelse(is.na(x), k, x) }
pest_label <- function(k) { x <- PEST_LABEL[k]; ifelse(is.na(x), k, x) }

use_rates <- tryCatch(if (file.exists("Pesticides_APVMA_use_rates.rds")) readRDS("Pesticides_APVMA_use_rates.rds") else NULL,
                      error = function(e) NULL)
raw_data$rate_all_use_kg_ha <- raw_data$MaximumPermittedPesticideApplicationRatekgPAIha
raw_data$has_ref_rate <- FALSE
if (!is.null(use_rates)) {
  k <- match(raw_data$PesticideactiveIngredient, use_rates$ref$active)
  k_ok <- !is.na(k) & is.finite(raw_data$Molecularweightgmol) & raw_data$Molecularweightgmol > 0
  raw_data$MaximumPermittedPesticideApplicationRatekgPAIha[k_ok]    <- use_rates$ref$ref_kg_ai_ha[k[k_ok]]
  raw_data$MaximumPermittedIngredientApplicationRateumolPAIha[k_ok] <- use_rates$ref$ref_kg_ai_ha[k[k_ok]] * 1e9 /
    raw_data$Molecularweightgmol[k_ok]
  raw_data$has_ref_rate[k_ok] <- TRUE
  rm(k, k_ok)
  use_rates$host$host_grp <- name_groups(use_rates$host$host, base = TRUE)$grp
  use_rates$pest$host_grp <- name_groups(use_rates$pest$host, base = TRUE)$grp
  use_rates$pest$pest_grp <- name_groups(use_rates$pest$pest)$grp
  if (!is.null(use_rates$rows)) use_rates$rows$host_grp <- name_groups(use_rates$rows$host, base = TRUE)$grp
}
pesticides_r <- tryCatch({
  r <- read.csv("Pesticides_r.csv", header = TRUE, encoding = "UTF-8", check.names = FALSE)
  r[, c("PesticideactiveIngredient", "General human health issues",
        "Handling issues", "WHO_classification")]
}, error = function(e) NULL)

if (!is.null(pesticides_r)) {
  names(pesticides_r) <- c("PesticideactiveIngredient",
                           "General_human_health_issues",
                           "Handling_issues", "WHO_classification")
  raw_data <- merge(raw_data, pesticides_r,
                    by = "PesticideactiveIngredient", all.x = TRUE,
                    suffixes = c("", ".r"))
  if (!"WHO_classification" %in% names(raw_data)) {
    raw_data$WHO_classification <- NA_character_
  }
}

# ---- Model constants ----
PRACTICE_BENEFIT_DEFAULT <- 100
PRACTICE_MULT_FLOOR      <- 0.20
DEFAULT_DAYS_TO_RUNOFF <- 21
ARROW_NOTE <- "Circles and arrows show the shift from your baseline (Step 6)."
PATHWAY_TAGS <- list(A = c("Aquatic", "#0173b2"), P = c("Pollinator", "#b8860b"), S = c("Soil", "#8c564b"))
affects <- function(...) {
  lapply(c(...), function(k) {
    b <- PATHWAY_TAGS[[k]]
    tags$span(style = sprintf(paste0("display:inline-block; font-size:9.5px; font-weight:600; color:%s; ",
                                     "border:1px solid %s; border-radius:8px; padding:0 6px; margin-left:4px; ",
                                     "line-height:15px; vertical-align:middle; background:#fff;"), b[2], b[2]),
              title = paste("Changes", tolower(b[1]), "risk"), b[1])
  })
}
affects_legend <- function() {
  tags$div(style = "font-size:13px; color:#444; margin:0 0 10px 0; line-height:1.8;",
           "Tags show which risks each setting affects:", affects("A", "P", "S"))
}
SOIL_TWA_DAYS <- 14

# ---- Runoff timing and first-event depth ----
RUNOFF_EVENT_MM  <- 5
YEAR_PERCENTILE  <- 0.80
RUNOFF_HORIZON_D <- 365
RAIN_YEARS       <- 1985:2024
AMC_DRY_MM       <- 35.6
AMC_WET_MM       <- 53.3
CN_ROWCROP <- c(A = 67, B = 78, C = 85, D = 89)
hsg_from_texture <- function(clay, sand) {
  if (is.null(clay) || is.na(clay)) return("C")
  if (clay > 40) "D" else if (clay >= 20) "C" else if (clay >= 10 || (!is.na(sand) && sand < 90)) "B" else "A"
}
Q_REF_MM <- 25

# ---- Pollinator exposure ----
POLL_EXPE_UG_PER_KGHA  <- 2.4
VAPOUR_DEP_DEFAULT_PCT <- 1
VOL_FIT                <- c(a = 0.192, b = 0.397)
HVAP_J_MOL             <- 95000
VP_REF_K               <- 293.15
F_AIR_METHOD           <- c(Boom = 0.10, Aircraft = 0.25, Vertical = 0.10)
DRIFT_DIST_M <- c(2, 4, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500, 600, 800)
DRIFT_CURVES <- list(
  "Vertical|<=2m" = c(0.151, 0.05124, 0.02723, 0.01739, 0.01228, 0.00924, 0.0059, 0.004166, 0.002214, 0.001414, 0.0009982, 0.0007512, 0.0004797, 0.0003387, 0.00018, 0.0001149, NA, NA, NA, NA, NA),
  "Vertical|>2m non-foliated" = c(0.3961, 0.2352, 0.1734, 0.1397, 0.1181, 0.103, 0.04751, 0.02768, 0.01037, 0.005169, 0.003011, 0.001937, 0.0009652, 0.0005623, 0.0002107, 0.000105, NA, NA, NA, NA, NA),
  "Vertical|>2m foliated" = c(0.2583, 0.1106, 0.06731, 0.04733, 0.03601, 0.02881, 0.01687, 0.01112, 0.005216, 0.003048, 0.00201, 0.00143, 0.0008356, 0.0005508, 0.0002584, 0.000151, NA, NA, NA, NA, NA),
  "Boom|FINE" = c(0.2834, 0.1979, 0.1504, 0.1193, 0.0972, 0.0813, 0.0603, 0.0468, 0.0321, 0.0224, 0.0171, 0.01348, 0.0096, 0.00772, 0.00556, 0.00434, 0.003, 0.00233, 0.00188, 0.00147, NA),
  "Boom|MEDIUM" = c(0.07445, 0.05008, 0.03828, 0.031, 0.02658, 0.02346, 0.01937, 0.01626, 0.009975, 0.006124, 0.004549, 0.003653, 0.002696, 0.00219, 0.00156, 0.001284, 0.000964, 0.0006949, 0.0005623, 0.0004293, NA),
  "Boom|COARSE" = c(0.03161, 0.02328, 0.0189, 0.01632, 0.01465, 0.01334, 0.01139, 0.009653, 0.005573, 0.003248, 0.002493, 0.002063, 0.001596, 0.001338, 0.0009599, 0.0007956, 0.0006024, 0.0004283, 0.0003535, 0.0002765, NA),
  "Boom|VERY COARSE" = c(0.01773, 0.01394, 0.01172, 0.01043, 0.009631, 0.009011, 0.007902, 0.006766, 0.003847, 0.002238, 0.001752, 0.001472, 0.001155, 0.0009813, 0.0007091, 0.0005874, 0.0004479, 0.0003141, 0.0002612, 0.0002055, NA),
  "Boom|EXTREMELY COARSE" = c(0.01211, 0.009818, 0.008436, 0.007611, 0.007158, 0.006753, 0.005991, 0.005179, 0.002901, 0.001672, 0.001314, 0.001111, 0.0008784, 0.0007448, 0.0005361, 0.0004437, 0.0003401, 0.000238, 0.0001994, 0.000158, NA),
  "Boom|ULTRA COARSE" = c(0.00663, 0.00561, 0.00494, 0.00455, 0.00436, 0.00418, 0.00378, 0.00328, 0.00178, 0.001, 0.000798, 0.000677, 0.000541, 0.000462, 0.000337, 0.00028, 0.000215, 0.000149, 0.000126, 0.0001, NA),
  "Fixed-wing|FINE" = c(0.8427, 0.7944, 0.7436, 0.6958, 0.6516, 0.6099, 0.5377, 0.482, 0.3883, 0.3267, 0.2813, 0.2447, 0.188, 0.1498, 0.09757, 0.0719, 0.04852, 0.03809, 0.03211, 0.0281, 0.02276),
  "Fixed-wing|MEDIUM" = c(0.8386, 0.7026, 0.5916, 0.5063, 0.4399, 0.3868, 0.3115, 0.2626, 0.1923, 0.153, 0.1266, 0.1049, 0.07317, 0.05389, 0.03109, 0.02098, 0.01267, 0.0094, 0.007713, 0.006645, 0.005253),
  "Fixed-wing|COARSE" = c(0.5722, 0.4413, 0.3491, 0.2845, 0.2375, 0.2022, 0.1554, 0.1269, 0.08904, 0.06945, 0.05685, 0.04652, 0.03221, 0.02345, 0.01271, 0.008343, 0.004884, 0.003659, 0.003029, 0.002617, 0.002054),
  "Fixed-wing|VERY COARSE" = c(0.3767, 0.2761, 0.2099, 0.166, 0.1354, 0.1134, 0.08524, 0.0687, 0.04744, 0.03826, 0.03075, 0.02493, 0.01719, 0.01243, 0.006768, 0.004426, 0.002546, 0.001922, 0.001601, 0.001393, 0.001105),
  "Fixed-wing|EXTREMELY COARSE" = c(0.2803, 0.1991, 0.1479, 0.115, 0.09259, 0.07673, 0.05691, 0.04545, 0.03467, 0.02612, 0.02038, 0.0163, 0.01114, 0.008016, 0.004469, 0.00305, 0.001781, 0.001344, 0.001075, 0.0009184, 0.0007228),
  "Fixed-wing|ULTRA COARSE" = c(0.1823, 0.124, 0.08902, 0.06734, 0.05317, 0.04347, 0.0348, 0.03, 0.02075, 0.01461, 0.01094, 0.008668, 0.00587, 0.004195, 0.002335, 0.001575, 0.0008865, 0.0006791, 0.0005427, 0.0004628, 0.0003625),
  "Helicopter|FINE" = c(0.8441, 0.8251, 0.7984, 0.7662, 0.7341, 0.7013, 0.6393, 0.5805, 0.454, 0.3583, 0.2881, 0.2373, 0.1694, 0.1266, 0.0744, 0.05075, 0.03, 0.0213, 0.0173, 0.01521, 0.0123),
  "Helicopter|MEDIUM" = c(0.984, 0.9448, 0.8807, 0.7994, 0.7129, 0.6283, 0.4843, 0.3772, 0.2223, 0.1471, 0.1051, 0.08, 0.0513, 0.0365, 0.0193, 0.0124, 0.00688, 0.0047, 0.0037, 0.00317, 0.00218),
  "Helicopter|COARSE" = c(1.02, 0.8817, 0.721, 0.5765, 0.4609, 0.3717, 0.2517, 0.1804, 0.0958, 0.0606, 0.0432, 0.0331, 0.0208, 0.0142, 0.00722, 0.00459, 0.00252, 0.00172, 0.0014, 0.00128, 0.000888),
  "Helicopter|VERY COARSE" = c(0.9297, 0.6923, 0.5038, 0.3706, 0.2789, 0.2155, 0.1382, 0.0963, 0.0499, 0.0335, 0.02366, 0.0177, 0.0107, 0.00737, 0.00384, 0.00249, 0.00139, 0.000961, 0.000785, 0.000721, 0.000506),
  "Helicopter|EXTREMELY COARSE" = c(0.8194, 0.5578, 0.3821, 0.2698, 0.1971, 0.1488, 0.0925, 0.0635, 0.035, 0.0226, 0.0157, 0.01146, 0.00696, 0.00482, 0.00253, 0.00164, 0.000919, 0.00063, 0.000512, 0.000469, 0.000325),
  "Helicopter|ULTRA COARSE" = c(0.6543, 0.4003, 0.2547, 0.1702, 0.1194, 0.0876, 0.0522, 0.0359, 0.0208, 0.01269, 0.00848, 0.00597, 0.00359, 0.00247, 0.00127, 0.000818, 0.000451, 0.000307, 0.000248, 0.000227, 0.000159),
  "Boom-high|FINE" = c(0.796, 0.6198, 0.5087, 0.4316, 0.3733, 0.327, 0.2574, 0.2084, 0.1331, 0.09182, 0.06687, 0.05068, 0.03237, 0.02184, 0.009723, 0.005443, 0.003204, 0.002208, NA, NA, NA),
  "Boom-high|MEDIUM" = c(0.4338, 0.3011, 0.2285, 0.1817, 0.1483, 0.1235, 0.08911, 0.06749, 0.04027, 0.02707, 0.01984, 0.01521, 0.009908, 0.007214, 0.004271, 0.003167, 0.002324, 0.001854, NA, NA, NA),
  "Boom-high|COARSE" = c(0.1973, 0.1375, 0.1043, 0.08262, 0.06719, 0.05579, 0.04173, 0.03341, 0.02148, 0.01414, 0.01032, 0.007866, 0.005293, 0.003967, 0.002458, 0.001923, 0.001549, 0.001209, NA, NA, NA),
  "Boom-high|VERY COARSE" = c(0.1048, 0.07521, 0.05776, 0.04698, 0.03911, 0.03419, 0.02719, 0.02243, 0.01386, 0.008771, 0.006311, 0.004916, 0.00345, 0.002674, 0.001787, 0.00145, 0.001188, 0.0009288, NA, NA, NA),
  "Boom-high|EXTREMELY COARSE" = c(0.06925, 0.05071, 0.03954, 0.03323, 0.02879, 0.02559, 0.02084, 0.01736, 0.01043, 0.006357, 0.004605, 0.003634, 0.002602, 0.002047, 0.001403, 0.001155, 0.0009561, 0.0007495, NA, NA, NA),
  "Boom-high|ULTRA COARSE" = c(0.03763, 0.03006, 0.02517, 0.02194, 0.01943, 0.01754, 0.01444, 0.01203, 0.006771, 0.003855, 0.002824, 0.002249, 0.001634, 0.001299, 0.0009067, 0.0007545, 0.0006276, 0.0004879, NA, NA, NA)
)
DRIFT_METHODS <- c("Ground boom, 0.5 m or lower" = "Boom", "Ground boom, 1.0 m or lower" = "Boom-high",
                   "Fixed-wing aircraft" = "Fixed-wing", "Helicopter" = "Helicopter",
                   "Vertical sprayer, canopy 2 m or smaller" = "Vertical|<=2m",
                   "Vertical sprayer, canopy over 2 m, not fully foliated" = "Vertical|>2m non-foliated",
                   "Vertical sprayer, canopy over 2 m, fully foliated" = "Vertical|>2m foliated")
SPRAY_QUALITIES <- c("Fine (106-235 \u00b5m)" = "FINE", "Medium (236-340 \u00b5m)" = "MEDIUM",
                     "Coarse (341-403 \u00b5m)" = "COARSE", "Very coarse (404-502 \u00b5m)" = "VERY COARSE",
                     "Extremely coarse (503-665 \u00b5m)" = "EXTREMELY COARSE", "Ultra coarse (over 665 \u00b5m)" = "ULTRA COARSE")
SPRAY_QUALITY_SHORT <- c(FINE = "Fine", MEDIUM = "Medium", COARSE = "Coarse", `VERY COARSE` = "Very coarse",
                         `EXTREMELY COARSE` = "Extremely coarse", `ULTRA COARSE` = "Ultra coarse")
drift_group <- function(method) if (grepl("^Vertical", method)) "Vertical" else
  if (method %in% c("Fixed-wing", "Helicopter")) "Aircraft" else "Boom"
drift_deposition <- function(method, quality, dist_m) {
  key <- if (grepl("^Vertical", method)) method else paste0(method, "|", quality)
  y <- DRIFT_CURVES[[key]]
  if (is.null(y)) return(NA_real_)
  ok <- !is.na(y)
  exp(approx(log(DRIFT_DIST_M[ok]), log(y[ok]), xout = log(dist_m), rule = 2)$y)
}
vp_ref_pa <- function(d) {
  calc <- d$Henrys_law_Pam3mol * d$Solubility_mgL / d$Molecularweightgmol
  if ("VapourPressure20C_mPa" %in% names(d)) ifelse(is.na(d$VapourPressure20C_mPa), calc, d$VapourPressure20C_mPa / 1000) else calc
}
raw_data$VP_ref_Pa <- vp_ref_pa(raw_data)
pollinator_dose <- function(rate_kg, VP_ref_Pa, temp_c, method, quality, dist_m, dep_pct) {
  dep  <- drift_deposition(method, quality, dist_m)
  VP_T <- VP_ref_Pa * exp(-HVAP_J_MOL / 8.314 * (1 / (temp_c + 273.15) - 1 / VP_REF_K))
  f_vol <- ifelse(is.na(VP_T) | VP_T <= 0, 0,
                  pmin(1, 10^(VOL_FIT[["a"]] + VOL_FIT[["b"]] * log10(pmax(VP_T * 1000, 1e-12))) / 100))
  airborne <- pmin(1, F_AIR_METHOD[[drift_group(method)]] + f_vol)
  POLL_EXPE_UG_PER_KGHA * rate_kg * (dep + dep_pct / 100 * airborne)
}

# ---- Irrigation ----
DEFAULT_IRRIG_DEPTH_MM   <- 95
DEFAULT_IRRIG_INTERVAL_D <- 23
DEFAULT_IRRIG_FIRST_D    <- 12

IRRIGATION_PARAMS <- data.frame(
  type            = c("furrow", "furrow_alt", "furrow_surge", "furrow_reuse", "sprinkler", "drip"),
  label           = c("Furrow - conventional (open-end, every row)", "Furrow - alternate row", "Furrow - surge",
                      "Furrow - tailwater captured and reused", "Sprinkler - centre pivot or lateral move",
                      "Drip - surface or subsurface"),
  SC_gL           = c(0.70, 0.70, 0.70, 0.70, 0.15, 0.03),
  runoff_fraction = c(0.25, 0.15, 0.10, 0.01, 0.07, 0.01),
  RF_range        = c("0.04-0.70", "0.10-0.20", "0.05-0.15", "0-0.02", "0-0.15", "0-0.02"),
  SC_source       = c("Spencer et al. 1985", "as furrow", "as furrow", "as furrow", "Ebbert and Kim 1998", "assumed"),
  RF_source       = c("Irmak et al. 2011; Roth et al. 2013", "Irmak et al. 2011", "Irmak et al. 2011",
                      "Irmak et al. 2011; nominal overflow", "Evans n.d.", "Irmak et al. 2011"),
  stringsAsFactors = FALSE
)
IRRIGATION_STUDIES <- data.frame(
  system = c("furrow", "furrow", "furrow", "furrow", "furrow", "furrow", "furrow", "furrow", "furrow", "furrow",
             "furrow", "furrow_alt", "furrow_surge", "sprinkler", "furrow"),
  id     = c("irmak_everyrow", "roth2013", "burdekin_typical", "burdekin_improved", "mackay_flood", "reddy46",
             "calif_chico", "calif_davis", "iran_cont", "brazil_ne", "idaho_silt_sc", "roth2013_alt", "iran_surge",
             "pivot_upper", "acpm2026_t52"),
  label  = c("Every-row furrow, 10-35% runoff (Irmak et al. 2011): 0.23",
             "Australian cotton, tailwater 4-32% (Roth et al. 2013): 0.18",
             "Lower Burdekin sugarcane, typical (Water Quality Science and Agricultural Hub n.d.): 0.20",
             "Lower Burdekin sugarcane, improved irrigation, half the typical loss (same source): 0.10",
             "Mackay-Whitsunday sugarcane flood irrigation, 13 of 70 mm (Rohde et al. 2013): 0.19",
             "Uzbekistan cotton, 46 events, mean 39% (Reddy et al. 2013): 0.39",
             "California, Chico site (Long et al. 2010): 0.19",
             "California, Davis site, mean 43% (Long et al. 2010): 0.43",
             "Iran, continuous furrow, tailwater 26.75% (Radmanesh et al. 2023): 0.27",
             "North-east Brazil smallholdings, 20-70% (Araujo et al. 2019): 0.45",
             "Idaho silt loam, sediment 6.2 g/L (Lentz et al. 1996); runoff as central value",
             "Australian cotton, 2 m alternate-row example (Roth et al. 2013): 0.19",
             "Iran, surge furrow, tailwater 19.6% (Radmanesh et al. 2023): 0.20",
             "Centre pivot, upper end of 0-15% (Evans n.d.): 0.15",
             "Australian cotton, 25 of 120 mm as tailwater (CRDC and CottonInfo 2026, Table 5.2): 0.21"),
  runoff_fraction = c(0.225, 0.18, 0.20, 0.10, 0.19, 0.39, 0.19, 0.43, 0.27, 0.45, NA, 0.19, 0.20, 0.15, 0.21),
  SC_gL           = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6.2, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)
irrigation_source_choices <- function(type) {
  s <- IRRIGATION_STUDIES[IRRIGATION_STUDIES$system == type, , drop = FALSE]
  c(setNames("central", "Central value for this system"), setNames(s$id, s$label))
}
irrigation_values <- function(type, source_id = "central") {
  ir <- IRRIGATION_PARAMS[IRRIGATION_PARAMS$type == type, , drop = FALSE]
  if (nrow(ir) != 1) return(NULL)
  st <- IRRIGATION_STUDIES[IRRIGATION_STUDIES$system == type & IRRIGATION_STUDIES$id == (source_id %||% "central"), ,
                           drop = FALSE]
  list(SC = if (nrow(st) == 1 && !is.na(st$SC_gL)) st$SC_gL else ir$SC_gL,
       runoff_fraction = if (nrow(st) == 1 && !is.na(st$runoff_fraction)) st$runoff_fraction else ir$runoff_fraction,
       source = if (nrow(st) == 1) st$label else "central value for this system")
}
irrigation_label <- function(type) {
  if (identical(type, "rainfall")) return("Rainfall-driven")
  lab <- IRRIGATION_PARAMS$label[IRRIGATION_PARAMS$type == type]
  if (length(lab) == 1) lab else type
}

# ---- Runoff extraction and sediment ----
ER_MAX    <- 7.0
ER_fun    <- function(SC, Csoil_mgkg) pmin(ER_MAX, pmax(1, 2.16 * (SC * Csoil_mgkg / 1000)^(-0.17)))
Fdiss_fun <- function(Kd, SC) 1 / (1 + Kd * SC / 1000)
ERO_fun   <- function(Kd, SC, Csoil_mgkg, ERO_base = 30) {
  fd <- Fdiss_fun(Kd, SC)
  ERO_base * (fd + (1 - fd) * ER_fun(SC, Csoil_mgkg))
}
SC_RAINFALL_MAX <- 200
SC_RAINFALL_MIN <- 2

# ---- Decay, runoff timing and leaching ----
q10_day <- function(T) 2.58^((T - 20) / 10)

annual_cn_runoff <- function(rain, CN2) {
  if (is.null(rain) || is.na(CN2)) return(NA_real_)
  P <- pmax(0, rain$P_mm); P[is.na(P)] <- 0; n <- length(P)
  CN1 <- 4.2 * CN2 / (10 - 0.058 * CN2); CN3 <- 23 * CN2 / (10 + 0.13 * CN2)
  cs <- c(0, cumsum(P)); w5 <- cs[seq_len(n)] - cs[pmax(1, seq_len(n) - 5)]
  cn <- ifelse(w5 < AMC_DRY_MM, CN1, ifelse(w5 > AMC_WET_MM, CN3, CN2))
  S <- 25400 / cn - 254; Ia <- 0.2 * S
  Q <- ifelse(P > Ia, (P - Ia)^2 / (P + 0.8 * S), 0)
  yr <- format(rain$date, "%Y"); full <- names(which(table(yr) >= 365))
  mean(tapply(Q, yr, sum)[full])
}

runoff_timing_fn <- function(rain, app_date, CN2, irr = NULL, tmonth = rep(25, 12), years = RAIN_YEARS,
                             horizon = RUNOFF_HORIZON_D, q_event = RUNOFF_EVENT_MM, twa_days = SOIL_TWA_DAYS) {
  if (is.null(rain) || nrow(rain) < 365 || is.na(CN2)) return(NULL)
  P_all <- pmax(0, rain$P_mm); P_all[is.na(P_all)] <- 0
  names(P_all) <- as.character(rain$date)
  md <- format(as.Date(app_date), "%m-%d"); if (md == "02-29") md <- "02-28"
  CN1 <- 4.2 * CN2 / (10 - 0.058 * CN2); CN3 <- 23 * CN2 / (10 + 0.13 * CN2)
  off <- (-5):horizon
  I <- rep(0, length(off))
  if (!is.null(irr) && irr$depth_mm > 0 && is.finite(irr$first_offset) && irr$first_offset <= horizon) {
    io <- seq(irr$first_offset, horizon, by = max(1, irr$interval))
    I[off %in% io] <- irr$depth_mm
  }
  frac <- if (is.null(irr)) 0 else irr$frac
  idx <- which(off >= 1)
  E14 <- INF14 <- NULL
  per_year <- lapply(years, function(y) {
    t0 <- as.Date(paste0(y, "-", md)); dts <- t0 + off
    p <- unname(P_all[as.character(dts)])
    if (any(is.na(p))) return(NULL)
    cs <- c(0, cumsum(p + I))
    w5 <- cs[idx] - cs[idx - 5]
    cn <- ifelse(w5 < AMC_DRY_MM, CN1, ifelse(w5 > AMC_WET_MM, CN3, CN2))
    S  <- 25400 / cn - 254; Ia <- 0.2 * S; pp <- p[idx]
    Qr <- ifelse(pp > Ia, (pp - Ia)^2 / (pp + 0.8 * S), 0)
    Qi <- I[idx] * frac
    e   <- cumsum(q10_day(tmonth[as.integer(format(dts[idx], "%m"))]))
    inf <- cumsum(pmax(0, pp + I[idx] - Qr - Qi))
    tw  <- list(E = c(0, e[seq_len(twa_days)]), INF = c(0, inf[seq_len(twa_days)]))
    hit <- which(Qr + Qi >= q_event)
    row <- if (length(hit) == 0)
      data.frame(year = y, days = NA_real_, E = NA_real_, infil = NA_real_, Q_rain = 0, Q_irr = 0) else {
      k <- hit[1]
      data.frame(year = y, days = k, E = e[k], infil = if (k > 1) inf[k - 1] else 0, Q_rain = Qr[k], Q_irr = Qi[k])
    }
    list(row = row, tw = tw)
  })
  per_year <- Filter(Negate(is.null), per_year)
  if (length(per_year) == 0) return(NULL)
  pe <- do.call(rbind, lapply(per_year, `[[`, "row"))
  E14   <- do.call(rbind, lapply(per_year, function(z) z$tw$E))
  INF14 <- do.call(rbind, lapply(per_year, function(z) z$tw$INF))
  dd <- ifelse(is.na(pe$days), Inf, pe$days)
  qr <- pe$Q_rain[is.finite(dd) & pe$Q_rain > 0]
  list(per_year = pe, E14 = E14, INF14 = INF14, n_years = nrow(pe), median_days = stats::median(dd),
       p7 = mean(dd <= 7), p14 = mean(dd <= 14), p30 = mean(dd <= 30),
       share_irr = mean(is.finite(dd) & pe$Q_irr >= pe$Q_rain),
       median_Q_rain = if (length(qr)) stats::median(qr) else NA_real_,
       median_Q_event = stats::median((pe$Q_rain + pe$Q_irr)[is.finite(dd)]))
}

LEACH_STORAGE <- 1.6
leach_matrix <- function(Kd, bd, INF) exp(-INF / ((pmax(Kd, 0) * bd + LEACH_STORAGE) * 25))

first_event_loss <- function(C0, DT50_base, Kd, bd, MW, SC_rain, SC_irr, ev) {
  n <- length(C0); m <- nrow(ev)
  none <- matrix(is.na(ev$days), n, m, byrow = TRUE)
  E    <- matrix(ifelse(is.na(ev$E), 0, ev$E), n, m, byrow = TRUE)
  INF  <- matrix(ifelse(is.na(ev$infil), 0, ev$infil), n, m, byrow = TRUE)
  Cst  <- pmax(C0 * exp(-log(2) * E / DT50_base) * leach_matrix(Kd, bd, INF), C0 * 0.05)
  Cmg  <- Cst * MW / 1000
  Qr   <- matrix(ev$Q_rain, n, m, byrow = TRUE); Qi <- matrix(ev$Q_irr, n, m, byrow = TRUE)
  L    <- Cst * ERO_fun(Kd, SC_rain, Cmg) * Qr
  if (!is.na(SC_irr) && any(Qi > 0)) L <- L + Cst * ERO_fun(Kd, SC_irr, Cmg) * Qi
  L[none] <- 0
  apply(L, 1, function(v) if (anyNA(v)) NA_real_ else stats::quantile(v, YEAR_PERCENTILE, names = FALSE))
}

soil_twa <- function(C0, DT50_base, Kd, bd, E14, INF14) {
  n <- length(C0); k <- ncol(E14)
  w <- c(0.5, rep(1, k - 2), 0.5) / (k - 1)
  out <- vapply(seq_len(nrow(E14)), function(j) {
    E   <- matrix(E14[j, ], n, k, byrow = TRUE); INF <- matrix(INF14[j, ], n, k, byrow = TRUE)
    Ct  <- pmax(C0 * exp(-log(2) * E / DT50_base) * leach_matrix(Kd, bd, INF), C0 * 0.05)
    as.numeric(Ct %*% w)
  }, numeric(n))
  if (is.null(dim(out))) out <- matrix(out, nrow = n)
  apply(out, 1, function(v) if (anyNA(v)) NA_real_ else stats::quantile(v, YEAR_PERCENTILE, names = FALSE))
}

# ---- NDVI window and risk categories ----
NDVI_WINDOW_DAYS <- 365
NDVI_END_DATE    <- NA_character_

risk_category <- function(E, T) {
  dplyr::case_when(
    is.na(E) | is.na(T) ~ "No data",
    (E >= 0.66 & T >= 0.66) ~ "Severe",
    ((E >= 0.66 & T >= 0.33 & T < 0.66) | (T >= 0.66 & E >= 0.33 & E < 0.66)) ~ "High",
    ((E >= 0.66 & T < 0.33) | (T >= 0.66 & E < 0.33) |
       (E >= 0.33 & E < 0.66 & T >= 0.33 & T < 0.66)) ~ "Moderate",
    ((E >= 0.33 & E < 0.66 & T < 0.33) | (T >= 0.33 & T < 0.66 & E < 0.33)) ~ "Low",
    TRUE ~ "Negligible"
  )
}

# ---- User interface ----
ui <- fluidPage(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      * {
        box-sizing: border-box;
      }

      body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background-color: #f5f7fa;
        color: #2c3e50;
        line-height: 1.6;
      }

      /* Wizard nav: keep sticky bar above Shiny's default content padding */
      .container-fluid { padding-top: 0 !important; }

      /* Make wizard back button look clickable without Bootstrap primary colour */
      .btn-wizard-back {
        font-size: 16px; font-weight: 600; padding: 10px 28px;
        background: #f5f5f5; border: 2px solid #ccc; border-radius: 8px;
        min-width: 120px; cursor: pointer;
      }
      .btn-wizard-back:hover { background: #e8e8e8; border-color: #aaa; }

      /* Fetch loading spinner */
      @keyframes _spin { 0%{transform:rotate(0deg)} 100%{transform:rotate(360deg)} }
      .fetch-spinner-ring {
        display: inline-block; width: 18px; height: 18px;
        border: 3px solid #c8e6fa; border-top-color: #0173b2;
        border-radius: 50%; animation: _spin 0.75s linear infinite;
        vertical-align: middle; margin-right: 8px;
      }
      #location_fetch_overlay {
        display: none;
        background: #e8f4ff; border: 1px solid #90caf9; border-radius: 6px;
        padding: 12px 16px; margin: 10px 0; font-size: 13px; color: #0173b2;
        font-weight: 500;
      }

      .app-header {
        background: linear-gradient(135deg, #0173b2 0%, #029e73 100%);
        color: white;
        padding: 25px 20px;
        margin-bottom: 30px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      }
      
      .app-header h1 {
        margin: 0;
        font-size: 28px;
        font-weight: 600;
      }
      
      .content-section {
        background-color: white;
        padding: 25px;
        margin-bottom: 20px;
        border-radius: 6px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08);
      }
      
      .load-box { background:#f0f8ff; border:1px solid #cce0f5; border-radius:6px; padding:12px 16px 2px 16px;
                  font-size:13px; color:#33475b; }
      .load-box .form-group { margin:8px 0 6px 0; }

      .disclaimer-box {
        background-color: #fff8e1;
        border-left: 4px solid #f59e0b;
        padding: 20px;
        margin: 25px 0;
        border-radius: 4px;
      }
      
      .disclaimer-box h3 {
        margin-top: 0;
        color: #92400e;
      }
      
      .disclaimer-box p, .disclaimer-box ul {
        color: #78350f;
        margin: 10px 0;
      }
      
      .btn {
        border: none;
        border-radius: 4px;
        font-size: 14px;
        padding: 8px 16px;
        cursor: pointer;
        transition: all 0.3s ease;
        font-weight: 500;
      }
      
      .btn:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
      }
      
      .btn-primary, .btn-success {
        background-color: #0173b2;
        color: white;
      }
      
      .btn-primary:hover, .btn-success:hover {
        background-color: #014d86;
      }
      
      .btn-default {
        background-color: #e8eef5;
        color: #2c3e50;
      }
      
      .btn-default:hover {
        background-color: #d4dce8;
      }
      
      .btn-info {
        background-color: #029e73;
        color: white;
      }
      
      .btn-info:hover {
        background-color: #01745f;
      }
      
      .btn-warning {
        background-color: #cc78bc;
        color: white;
      }
      
      .btn-warning:hover {
        background-color: #b85fa8;
      }
      
      .btn-lg {
        padding: 12px 30px;
        font-size: 16px;
      }
      
      .form-group {
        margin-bottom: 20px;
      }
      
      .form-group label {
        font-weight: 500;
        color: #2c3e50;
        display: block;
        margin-bottom: 8px;
      }
      
      .form-group input,
      .form-group select {
        width: 100%;
        padding: 10px;
        border: 1px solid #d4dce8;
        border-radius: 4px;
        font-size: 14px;
        font-family: inherit;
      }
      
      .form-group input:focus,
      .form-group select:focus {
        outline: none;
        border-color: #0173b2;
        box-shadow: 0 0 0 3px rgba(1, 115, 178, 0.1);
      }
      
      .risk-card {
        background-color: white;
        border: 1px solid #e8eef5;
        border-radius: 6px;
        padding: 20px;
        margin: 10px 0;
        text-align: center;
      }
      
      .risk-card h4 {
        margin: 0 0 10px 0;
        font-size: 16px;
      }
      
      .risk-score {
        font-size: 32px;
        font-weight: bold;
        margin: 10px 0;
      }
      
      .risk-low {
        color: #029e73;
        border-top: 4px solid #029e73;
      }
      
      .risk-medium {
        color: #ee7733;
        border-top: 4px solid #ee7733;
      }
      
      .risk-high {
        color: #cc3311;
        border-top: 4px solid #cc3311;
      }
      
      .table {
        font-size: 14px;
      }
      
      .dataTables_wrapper {
        font-size: 13px;
      }
      
      /* Fix checkbox text overlap */
      .shiny-options-group {
        margin: 8px 0;
      }
      
      .shiny-options-group .checkbox {
        margin-bottom: 10px;
        display: block;
      }
      
      .shiny-options-group .checkbox label {
        display: flex;
        align-items: flex-start;
        gap: 10px;
        margin: 0;
        padding: 0;
        white-space: normal;
        font-weight: normal;
      }
      
      .shiny-options-group .checkbox input[type='checkbox'] {
        flex: 0 0 18px;
        height: 18px;
        margin-top: 2px;
        cursor: pointer;
      }
      
      .shiny-options-group .checkbox span {
        flex: 1;
        word-break: break-word;
        line-height: 1.4;
      }
      
      /* Collapsible table controls */
      td.dt-control {
        cursor: pointer;
        padding: 8px 12px !important;
        text-align: center;
        user-select: none;
        font-weight: bold;
        color: #0173b2;
      }
      
      td.dt-control:hover {
        background-color: #f0f8ff;
      }
      
      /* Risk table clickable rows */
      .risk-summary-table tbody tr {
        cursor: pointer;
      }
      
      .risk-summary-table tbody tr:hover {
        background-color: #f5f5f5 !important;
      }
      
      @media (max-width: 768px) {
        .app-header h1 {
          font-size: 22px;
        }
        .content-section {
          padding: 15px;
        }
      }
      /* ── Mobile-responsive wizard navigation bar ── */
      @media (max-width: 640px) {
        .wiz-outer {
          flex-wrap: wrap !important;
          padding: 8px 12px !important;
          row-gap: 6px;
        }
        .wiz-steps {
          order: -1 !important;
          width: 100% !important;
          justify-content: center !important;
        }
        .wiz-btn {
          min-width: 70px !important;
          padding: 8px 14px !important;
          font-size: 14px !important;
        }
      }
    ")),
    tags$style(HTML("
  /* Compact DataTable styling */
  .dataTable tbody tr {
    height: 35px !important;
  }
  
  .dataTable tbody td {
    padding: 4px 6px !important;
    line-height: 1.2 !important;
    white-space: normal !important;
    word-wrap: break-word !important;
    max-height: 35px !important;
    overflow: hidden !important;
    text-overflow: ellipsis !important;
    vertical-align: middle !important;
  }
  
  .dataTable thead th {
    padding: 6px 8px !important;
  }
  
  /* Risk summary table clickable first column */
  .risk-summary-table tbody td:first-child {
    cursor: pointer;
    font-weight: bold;
    color: #0173b2;
    text-decoration: underline;
  }
  
  .risk-summary-table tbody td:first-child:hover {
    background-color: #e6f2ff !important;
  }

  /* Pest group tile grid */
  .pest-tile-grid { display:grid; grid-template-columns:repeat(3,1fr); gap:14px; margin-bottom:20px; }
  .pest-tile { background:#fff; border:2px solid #d4dce8; border-radius:10px; padding:16px 14px;
    cursor:pointer; text-align:left; width:100%; min-height:110px; position:relative;
    transition:border-color 0.15s, background 0.15s, box-shadow 0.15s; }
  .pest-tile:hover { border-color:#01745f; box-shadow:0 2px 8px rgba(1,116,95,0.18); }
  .pest-tile-selected { border-color:#01745f !important; background:#f0fdf9 !important; }
  .pest-tile-selected::after { content:'\\2714'; position:absolute; top:8px; right:10px;
    color:#01745f; font-size:16px; font-weight:bold; }
  @media (max-width:680px) { .pest-tile-grid { grid-template-columns:repeat(2,1fr); } }
  @media (max-width:420px) { .pest-tile-grid { grid-template-columns:1fr; } }
  .info-tip { display:inline-block; width:15px; height:15px; line-height:15px; border-radius:50%; margin-left:5px;
              background:#0173b2; color:#fff; font:italic 700 11px Georgia, 'Times New Roman', serif; text-align:center;
              vertical-align:1px; cursor:pointer; user-select:none; }
  .info-tip:hover, .info-tip:focus { background:#01745f; outline:none; }
  .popover { font-size:12px; max-width:320px; }
  .ref-list { margin-top:22px; padding-top:10px; border-top:1px solid #e3e8ee; font-size:11px; color:#555; }
  .ref-list .ref-title { font-weight:700; color:#33475b; margin-bottom:4px; }
  .ref-list p { margin:0 0 4px 0; padding-left:14px; text-indent:-14px; line-height:1.45; }
  .src-line { font-size:11px; color:#777; margin:0 0 10px 0; }
"))

  ),
  shinyjs::useShinyjs(),
  tags$script(HTML("$(function() { $('body').popover({ selector: '[data-toggle=\\'popover\\']', trigger: 'hover focus', container: 'body' }); });")),

  div(class = "app-header",
      h1("Pesticide Risk Assessment Tool")
  ),

  uiOutput("ui_pages")
)

pest_group_regex <- list(
  grass_weeds     = "grass(es|y)?\\b|\\bryegrass|\\boat(e?s)?\\b|\\bbarley grass|\\bbrome|\\bfescue|\\bphalaris|\\bcanary grass|\\bcrowfoot|\\bcrabgrass|\\bliverseed|\\bpanic(e?s)?\\b|\\bpaspalum|\\bsetaria|\\bsilvergrass|\\bwindmill grass",
  broadleaf_weeds = "\\bbroadleaf|\\bthistle(e?s)?\\b|\\bpigweed|\\bbindweed|\\bragweed|\\bdock(e?s)?\\b|\\bdandelion|\\bclover(e?s)?\\b|\\bspurge(e?s)?\\b|\\blambsquarters|\\bnightshade|\\bmorning glory|\\bwild radish|\\bcharlock|\\bcapeweed|\\bmustard(e?s)?\\b|\\bsorrel(e?s)?\\b|\\bplantain(e?s)?\\b|\\bchickweed|\\bfleabane|\\bamaranth|\\bparthenium|\\bpaddy melon|\\bprickly acacia",
  sedges_woody    = "\\bsedge(e?s)?\\b|\\bnutgrass|\\bnutsedge|\\brush(e?s)?\\b|\\bwoody weed|\\bblackberr|\\blantana|\\bbrush\\b|\\bregrowth|\\bvine(e?s)?\\b|\\bprickly pear|\\bmesquite|\\bmimosa|\\bacacia|\\bwattle|\\bstump",
  sucking_insects = "\\baphid(e?s)?\\b|\\bmite(e?s)?\\b|\\bwhitefl|\\bscale(e?s)?\\b|\\bmealybug|\\bmealy bug|\\bthrip|\\bleafhopper|\\bplanthopper|\\bpsyllid(e?s)?\\b|\\bmirid(e?s)?\\b|\\blygus|\\bcotton bug|\\bstink bug|\\bgreen vegetable bug|\\bspotting bug|\\bcicada",
  caterpillars    = "\\bcaterpillar|\\barmyworm|\\barmy worm|\\bbollworm|\\bbudworm|\\bcutworm|\\bmoth(e?s)?\\b|\\blooper(e?s)?\\b|\\bheliothis|\\bhelicoverpa|\\bsemilooper|\\bleaf roller|\\bleafroller|\\bearworm|\\bborer(e?s)?\\b|\\bdiamondback|\\blarva",
  beetles_grubs   = "\\bbeetle(e?s)?\\b|\\bweevil(e?s)?\\b|\\bgrub(e?s)?\\b|\\bcanegrub|\\bgrasshopper|\\blocust(e?s)?\\b|\\bcricket(e?s)?\\b|\\bfruit fl(y|ies)\\b|\\bleafminer|\\bleaf miner|\\bant(e?s)?\\b",
  soil_pests      = "\\bwireworm|\\bnematode|\\btermite|\\brootworm|\\bsoil insect|\\bmillipede|\\bspringtail|\\bsymphylid|\\bslug(e?s)?\\b|\\bsnail(e?s)?\\b|\\bearth mite|\\bredlegged",
  leaf_diseases   = "\\brust(e?s)?\\b|\\bpowdery mildew|\\bdowny mildew|\\bbotrytis|\\bgrey mould|\\banthracnose|\\bblight(e?s)?\\b|\\balternaria|\\bcercospora|\\bsclerotinia|\\bsmut(e?s)?\\b|\\bleaf spot|\\bscab\\b|\\bblack spot|\\bsigatoka|\\bseptoria|\\bnet blotch|\\bscald\\b",
  root_diseases   = "\\bphytophthora|\\bpythium|\\bfusarium|\\brhizoctonia|\\broot rot|\\bcollar rot|\\bdamping.off|\\bdieback|\\bverticillium|\\bcrown rot|\\bsoil.borne|\\btake.all"
)

# ---- References ----
REFS <- list(
  apvma_pubcris = "APVMA 2026a. Public Chemical Registration Information System Search (PubCRIS). Australian Pesticides and Veterinary Medicines Authority, Canberra, Australia. https://portal.apvma.gov.au/pubcris (data accessed 22 October 2026).",
  apvma_labels  = "APVMA 2026b. Product labels (e-labels). Australian Pesticides and Veterinary Medicines Authority, Canberra, Australia. https://elabels.apvma.gov.au (accessed 21 October 2026).",
  apvma_sdram   = "APVMA 2019. Spray drift risk assessment manual, Stage one (July 2019). Australian Pesticides and Veterinary Medicines Authority, Canberra, Australia.",
  apvma_sdrat   = "APVMA 2023. Spray drift risk assessment tool (SDRAT), version 1.1 (June 2023). Microsoft Excel workbook. Australian Pesticides and Veterinary Medicines Authority, Canberra, Australia.",
  asabe         = "ASABE 2009. S572.1 Spray nozzle classification by droplet spectra. American Society of Agricultural and Biological Engineers, St Joseph, MI, USA.",
  araujo        = "Araujo, D.F., Costa, R.N. and Mateos, L. 2019. Pros and cons of furrow irrigation on smallholdings in northeast Brazil. Agricultural Water Management 221, 25-33.",
  crdc          = "CRDC and CottonInfo 2026. Australian Cotton Production Manual 2026-27. Cotton Research and Development Corporation, Narrabri, Australia. ISSN 2206-0081.",
  ebbert        = "Ebbert, J.C. and Kim, M.H. 1998. Relation between irrigation method, sediment yields, and losses of pesticides and nitrogen. Journal of Environmental Quality 27(2), 372-380.",
  evans_pivot   = "Evans, R.G. n.d. Center pivot irrigation. Northern Plains Agricultural Research Laboratory, United States Department of Agriculture, Agricultural Research Service, Sidney, MT, USA. https://www.ars.usda.gov/ARSUserFiles/21563/center%20pivot%20design%202.pdf",
  funk          = "Funk, C., Peterson, P., Landsfeld, M., Pedreros, D., Verdin, J., Shukla, S., Husak, G., Rowland, J., Harrison, L., Hoell, A. and Michaelsen, J. 2015. The climate hazards infrared precipitation with stations: a new environmental record for monitoring extremes. Scientific Data 2, 150066.",
  gentil        = "Gentil-Sergent, C., Basset-Mens, C., Gaab, J., Mottes, C., Melero, C. and Fantke, P. 2021. Quantifying pesticide emission fractions for tropical conditions. Chemosphere 275, 130014.",
  grundy        = "Grundy, M.J., Viscarra Rossel, R.A., Searle, R.D., Wilson, P.L., Chen, C. and Gregory, L.J. 2015. Soil and Landscape Grid of Australia. Soil Research 53(8), 835-844.",
  hengl         = "Hengl, T., Parente, L., Ho, Y.-F. and Simoes, R. 2023. OpenLandMap Open Land Data services, version 0.2. OpenGeoHub Foundation, Wageningen, the Netherlands. https://doi.org/10.5281/zenodo.10522799",
  houbraken     = "Houbraken, M., van den Berg, F., Butler Ellis, C.M., Dekeyser, D., Nuyttens, D., De Schampheleire, M. and Spanoghe, P. 2016. Volatilisation of pesticides under field conditions: inverse modelling and pesticide fate models. Pest Management Science 72(7), 1309-1321.",
  fleming       = "Fleming, D.E., Spencer, G.D. and Krutz, L.J. 2025. Pesticide runoff from conventional tillage, minimum tillage, and no-tillage cropping systems: Meta-analysis of published North American data. Journal of Environmental Quality 54(5), 913-924.",
  irmak         = "Irmak, S., Odhiambo, L.O., Kranz, W.L. and Eisenhauer, D.E. 2011. Irrigation efficiency and uniformity, and crop water use efficiency. Extension Circular EC732, University of Nebraska-Lincoln, Lincoln, NE, USA.",
  klein         = "Klein, M., Klein, J., Flade, J., Grossmann, D., T\u00fcrkowsky, D., O'Connor, I., Spycher, S., Reichenberger, S., Sittig, S., Multsch, S. and Thomas, K. 2023. Risk mitigation measures for pesticide runoff: How effective are they? Pest Management Science 79(12), 4897-4905.",
  lentz         = "Lentz, R.D., Sojka, R.E. and Carter, D.L. 1996. Furrow irrigation water-quality effects on soil loss and infiltration. Soil Science Society of America Journal 60(1), 238-245.",
  linke         = "Linke, S., Lehner, B., Ouellet Dallaire, C., Ariwi, J., Grill, G., Anand, M., Beames, P., Burchard-Levine, V., Maxwell, S., Moidu, H., Tan, F. and Thieme, M. 2019. Global hydro-environmental sub-basin and river reach characteristics at high spatial resolution. Scientific Data 6(1), 283.",
  long          = "Long, R.F., Hanson, B.R., Fulton, A.E. and Weston, D.P. 2010. Mitigation techniques reduce sediment in runoff from furrow-irrigated cropland. California Agriculture 64(3), 135-140.",
  munoz         = "Mu\u00f1oz-Sabater, J., Dutra, E., Agust\u00ed-Panareda, A., Albergel, C., Arduini, G., Balsamo, G., Boussetta, S., Choulga, M., Harrigan, S., Hersbach, H., Martens, B., Miralles, D.G., Piles, M., Rodr\u00edguez-Fern\u00e1ndez, N.J., Zsoter, E., Buontempo, C. and Th\u00e9paut, J.-N. 2021. ERA5-Land: a state-of-the-art global reanalysis dataset for land applications. Earth System Science Data 13(9), 4349-4383.",
  nrcs2004      = "USDA-NRCS 2004. National Engineering Handbook, Part 630 Hydrology, Chapter 10: Estimation of direct runoff from storm rainfall. United States Department of Agriculture, Natural Resources Conservation Service, Washington, DC, USA.",
  nrcs2009      = "USDA-NRCS 2009. National Engineering Handbook, Part 630 Hydrology, Chapter 7: Hydrologic soil groups. United States Department of Agriculture, Natural Resources Conservation Service, Washington, DC, USA.",
  poggio        = "Poggio, L., de Sousa, L.M., Batjes, N.H., Heuvelink, G.B.M., Kempen, B., Ribeiro, E. and Rossiter, D. 2021. SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. SOIL 7(1), 217-240.",
  radmanesh     = "Radmanesh, M., Ahmadi, S.H. and Sepaskhah, A.R. 2023. Measurement and simulation of irrigation performance in continuous and surge furrow irrigation using WinSRFR and SIRMOD models. Scientific Reports 13, 5768.",
  reddy         = "Reddy, J.M., Jumaboev, K., Matyakubov, B. and Eshmuratov, D. 2013. Evaluation of furrow irrigation practices in Fergana Valley of Uzbekistan. Agricultural Water Management 117, 133-144.",
  rohde         = "Rohde, K., McDuffie, K. and Agnew, J. 2013. Paddock to sub-catchment scale water quality monitoring of sugarcane management practices. Interim report 2012/13 wet season, Mackay Whitsunday region. Department of Natural Resources and Mines, Queensland Government, for Reef Catchments (Mackay Whitsunday Isaac) Limited, Mackay, Australia.",
  roth          = "Roth, G., Harris, G., Gillies, M., Montgomery, J. and Wigginton, D. 2013. Water-use efficiency and productivity trends in Australian irrigated cotton: a review. Crop & Pasture Science 64(12), 1033-1048.",
  silburn2013   = "Silburn, D.M., Foley, J.L. and deVoil, R.C. 2013. Managing runoff of herbicides under rainfall and furrow irrigation with wheel traffic and banded spraying. Agriculture, Ecosystems & Environment 180, 40-53.",
  spencer       = "Spencer, W.F., Cliath, M.M., Blair, J.W. and LeMert, R.A. 1985. Transport of pesticides from irrigated fields in surface runoff and tile drain waters. Conservation Research Report 31, United States Department of Agriculture, Agricultural Research Service, Washington, DC, USA.",
  usepa         = "US EPA, PMRA and CDPR 2014. Guidance for assessing pesticide risks to bees. United States Environmental Protection Agency, Health Canada Pest Management Regulatory Agency and California Department of Pesticide Regulation, Washington, DC, USA.",
  wqsah_burdekin = "Water Quality Science and Agricultural Hub n.d. Lower Burdekin. https://www.agwaterscience.com/sugarcane-region/lower-burdekin (accessed 28 October 2026).",
  silburn2024   = "Silburn, D.M. and Crawford, M.H. 2024. Pesticide extraction from soil into runoff in North American and Australian croplands. Soil Research 62(3), SR23148.",
  ppdb          = "Lewis, K.A., Tzilivakis, J., Warner, D.J. and Green, A. 2016. An international database for pesticide risk assessments and management. Human and Ecological Risk Assessment 22(4), 1050-1064.",
  anzg          = "ANZG 2018. Australian and New Zealand Guidelines for Fresh and Marine Water Quality. Australian and New Zealand Governments and Australian state and territory governments, Canberra, Australia. https://www.waterquality.gov.au/anz-guidelines",
  zhang         = "Zhang, X., Liu, X., Zhang, M., Dahlgren, R.A. and Eitzel, M. 2010. A review of vegetated buffers and a meta-analysis of their mitigation efficacy in reducing nonpoint source pollution. Journal of Environmental Quality 39(1), 76-84."
)
ref_list <- function(keys, title = "References") {
  r <- sort(unlist(REFS[unique(keys)]))
  tags$div(class = "ref-list", tags$div(class = "ref-title", title), lapply(r, function(x) tags$p(x)))
}
BASELINE_SOURCES <- tags$div(
  tags$div(class = "ref-list",
    tags$div(class = "ref-title", "Sources"),
    tags$p(tags$b("Rate and timing: "), "runoff timing from daily rainfall and the runoff curve number (USDA-NRCS 2004, 2009)."),
    tags$p(tags$b("Spraying and drift: "), "drift deposition from the APVMA standard drift curves (APVMA 2019, 2023); ",
           "bee contact dose (US EPA, PMRA and CDPR 2014); airborne spray and vapour (Gentil-Sergent et al. 2021; ",
           "Houbraken et al. 2016); droplet size classes (ASABE 2009)."),
    tags$p(tags$b("Irrigation: "), "runoff fractions (Irmak et al. 2011; Roth et al. 2013; Rohde et al. 2013; Reddy et al. 2013; ",
           "Long et al. 2010; Araujo et al. 2019; Radmanesh et al. 2023; Evans n.d.; Water Quality Science and Agricultural Hub n.d.; ",
           "CRDC and CottonInfo 2026); sediment in irrigation runoff (Spencer et al. 1985; Lentz et al. 1996; Ebbert and Kim 1998)."),
    tags$p(tags$b("Runoff controls: "), "filter strips (Zhang et al. 2010); tillage (Fleming et al. 2025); micro-dams (Klein et al. 2023).")),
  ref_list(c("nrcs2004", "nrcs2009", "apvma_sdram", "apvma_sdrat", "usepa", "gentil", "houbraken", "asabe",
             "irmak", "roth", "rohde", "reddy", "long", "araujo", "radmanesh", "evans_pivot", "wqsah_burdekin", "crdc",
             "spencer", "lentz", "ebbert", "zhang", "klein", "fleming")))

RESULTS_SOURCES <- tags$details(style = "margin-top:16px;", tags$summary("Sources and references"),
  tags$div(class = "ref-list",
    tags$p(tags$b("Chemical properties: "), "PPDB (Lewis et al. 2016); aquatic guideline values from ANZG (2018), otherwise ",
           "estimated from the lowest acute EC50 or LC50 across fish, invertebrates, algae and plants."),
    tags$p(tags$b("Exposure: "), "runoff and enrichment (Silburn and Crawford 2024); drift and bee dose (APVMA 2019, 2023; ",
           "US EPA, PMRA and CDPR 2014); airborne spray and vapour (Gentil-Sergent et al. 2021; Houbraken et al. 2016)."),
    tags$p(tags$b("Runoff controls: "), "filter strips (Zhang et al. 2010); tillage (Fleming et al. 2025); micro-dams (Klein et al. 2023).")),
  ref_list(c("ppdb", "anzg", "silburn2024", "apvma_sdram", "apvma_sdrat", "usepa", "gentil", "houbraken", "zhang", "klein", "fleming")))

# ---- Helpers and saved-assessment aliases ----
conc_text <- function(amount, unit) {
  a <- suppressWarnings(as.numeric(amount)); u <- gsub("\\s*/\\s*", "/", trimws(as.character(unit)))
  ifelse(is.na(a) | is.na(u) | !nzchar(u), "", paste(vapply(a, function(x) format(x, drop0trailing = TRUE, trim = TRUE), ""), u))
}
conc_label <- function(amount, unit) { x <- conc_text(amount, unit); ifelse(nzchar(x), paste0(" (", x, ")"), "") }

info_tip <- function(text) tags$span(class = "info-tip", tabindex = "0", role = "button", `data-toggle` = "popover",
                                     `data-trigger` = "hover focus", `data-placement` = "top", `data-container` = "body",
                                     `data-content` = text, "i")

ACTIVE_RENAMES <- c("Sulphur" = "Sulfur", "2 2 Dpa" = "2,2-DPA", "2 4 Db" = "2,4-DB",
                    "1 3 Dichloropropene" = "1,3-Dichloropropene", "1 Methylcyclopropene" = "1-Methylcyclopropene",
                    "2 4 D Isobutyl" = "2,4-D isobutyl ester", "Lambda Cyhalothrin" = "Lambda-cyhalothrin",
                    "Metalaxyl M" = "Metalaxyl-M")
TILLAGE_ALIASES <- c(Ridge_Till = "Minimum_Till", Chisel = "Minimum_Till", Conservation_Till = "Minimum_Till")
PRACTICES_KEPT  <- "Micro_Dams"
ACTIVE_MERGES <- c("S Indoxacarb" = "Indoxacarb")
PEST_GROUP_ALIASES <- list(chewing_insects = c("caterpillars", "beetles_grubs"), sedges_rushes = "sedges_woody",
                           fungal_disease = c("leaf_diseases", "root_diseases"))
upgrade_pest_groups <- function(g) unique(unlist(lapply(g, function(x) PEST_GROUP_ALIASES[[x]] %||% x)))
PEST_GROUP_CLASS <- list(grass_weeds = "Herbicide", broadleaf_weeds = "Herbicide", sedges_woody = "Herbicide",
                         sucking_insects = c("Insecticide", "Miticide", "Mixed function pesticide"),
                         caterpillars = c("Insecticide", "Mixed function pesticide"),
                         beetles_grubs = c("Insecticide", "Mixed function pesticide"),
                         soil_pests = c("Insecticide", "Nematicide", "Miticide", "Mixed function pesticide"),
                         leaf_diseases = c("Fungicide", "Mixed function pesticide"),
                         root_diseases = c("Fungicide", "Mixed function pesticide"))
GROUP_ACTIVES <- local({
  ex <- lapply(pest_group_regex, function(rx)
    unique(raw_data$PesticideactiveIngredient[grepl(rx, raw_data$Example_pests_controlled, ignore.case = TRUE, perl = TRUE)]))
  if (is.null(uses_raw)) return(ex)
  cls <- tapply(raw_data$Class, raw_data$PesticideactiveIngredient, function(x) x[1])
  pa  <- unique(raw_data[!is.na(raw_data$pcode_int), c("pcode_int", "PesticideactiveIngredient")])
  pl  <- levels(uses_raw$pest); pi <- as.integer(uses_raw$pest)
  setNames(lapply(names(pest_group_regex), function(g) {
    ix <- which(grepl(pest_group_regex[[g]], pl, ignore.case = TRUE, perl = TRUE))
    pc <- unique(uses_raw$pcode[pi %in% ix])
    a  <- unique(pa$PesticideactiveIngredient[pa$pcode_int %in% pc])
    union(ex[[g]], a[cls[a] %in% PEST_GROUP_CLASS[[g]]])
  }), names(pest_group_regex))
})

# ---- Server ----
server <- function(input, output, session) {
  eps  <- 1e-9

  current_page <- reactiveVal("landing")

  selected_pest_groups <- reactiveVal(character(0))

  active_pest_regex <- reactive({
    groups <- selected_pest_groups()
    if (length(groups) == 0) return(NULL)
    paste(pest_group_regex[groups], collapse = "|")
  })

  # ---- Soil state ----
  sel_om <- reactiveVal("moderate")
  sel_ph <- reactiveVal(7.0)

  # ---- User overrides ----
  user_runoff_frac <- reactiveVal(NULL)
  rain_daily_data  <- reactiveVal(NULL)
  user_cn          <- reactiveVal(NULL)
  user_clay_pct    <- reactiveVal(NULL)
  user_temp_c      <- reactiveVal(NULL)
  user_slope_dg    <- reactiveVal(NULL)
  user_sand_pct    <- reactiveVal(NULL)
  user_silt_pct    <- reactiveVal(NULL)
  user_bd_gcm3     <- reactiveVal(NULL)
  user_soc_pct     <- reactiveVal(NULL)

  # ---- Site data state ----
  map_coords      <- reactiveVal(NULL)
  hydroatlas_data <- reactiveVal(NULL)
  soilgrids_data   <- reactiveVal(NULL)
  soil_source_name <- reactiveVal(NULL)
  soil_source_id   <- reactiveVal(NULL)
  SOIL_SOURCES <- list(slga = list(name = "the Soil and Landscape Grid of Australia", cite = "Grundy et al. 2015", ref = "grundy"),
                       soilgrids = list(name = "SoilGrids", cite = "Poggio et al. 2021", ref = "poggio"),
                       olm = list(name = "OpenLandMap", cite = "Hengl et al. 2023", ref = "hengl"))
  soil_src <- function() {
    id <- soil_source_id()
    if (is.null(id)) { nm <- soil_source_name() %||% ""
      id <- if (grepl("SLGA", nm)) "slga" else if (grepl("OpenLandMap", nm)) "olm" else "soilgrids" }
    SOIL_SOURCES[[id]]
  }
  rusle_data        <- reactiveVal(NULL)

  gee_token_cache <- new.env()
  gee_token_fn <- function(key_file) {
    if (!is.null(gee_token_cache$token) && identical(gee_token_cache$key, key_file) &&
        Sys.time() < gee_token_cache$expires) return(gee_token_cache$token)
    key_data <- jsonlite::read_json(key_file)
    claims <- jose::jwt_claim(
      iss   = key_data$client_email,
      scope = "https://www.googleapis.com/auth/earthengine",
      aud   = "https://oauth2.googleapis.com/token",
      exp   = as.integer(Sys.time()) + 3600,
      iat   = as.integer(Sys.time())
    )
    jwt  <- jose::jwt_encode_sig(claims, key_data$private_key)
    resp <- httr::POST(
      "https://oauth2.googleapis.com/token",
      body   = list(grant_type = "urn:ietf:params:oauth:grant-type:jwt-bearer",
                    assertion  = jwt),
      encode = "form", httr::timeout(30)
    )
    tok <- httr::content(resp)$access_token
    if (httr::http_error(resp) || is.null(tok)) stop("Earth Engine sign-in failed (check the service-account key)")
    gee_token_cache$token <- tok; gee_token_cache$key <- key_file; gee_token_cache$expires <- Sys.time() + 50 * 60
    tok
  }

  fetch_hydroatlas_fn <- function(lat, lon, key_file,
                                  asset = "WWF/HydroATLAS/v1/Basins/level12") {
    token      <- gee_token_fn(key_file)
    region <- if (!is.null(fetch_geom)) list(type = "Polygon", coordinates = list(fetch_geom)) else
      list(type = "Point", coordinates = list(lon, lat))
    point_json <- jsonlite::toJSON(region, auto_unbox = TRUE, digits = NA)
    resp <- httr::GET(
      paste0("https://earthengine.googleapis.com/v1/projects/earthengine-public/assets/",
             asset, ":listFeatures"),
      httr::add_headers(Authorization = paste("Bearer", token)),
      query = list(region = point_json), httr::timeout(60)
    )
    if (httr::http_error(resp)) stop(sprintf("HydroATLAS request failed (HTTP %d)", httr::status_code(resp)))
    res <- httr::content(resp)
    fs <- res$features
    if (length(fs) > 1) {
      pr <- lapply(fs, `[[`, "properties"); nm <- names(pr[[1]])
      avg <- lapply(setNames(nm, nm), function(k) {
        v <- lapply(pr, `[[`, k)
        if (all(vapply(v, function(z) is.numeric(z) && length(z) == 1, TRUE))) mean(unlist(v)) else v[[1]]
      })
      res$features <- list(list(properties = avg))
      res$n_subbasins <- length(fs)
    }
    res
  }

  parse_hydroatlas <- function(result) {
    props   <- result$features[[1]]$properties
    pre     <- props$pre_mm_syr %||% NA_real_
    run     <- props$run_mm_syr %||% NA_real_
    aet     <- props$aet_mm_syr %||% NA_real_
    tmp_raw <- props$tmp_dc_syr %||% NA_real_

    pre_monthly <- sapply(sprintf("pre_mm_s%02d", 1:12),
                          function(b) props[[b]] %||% NA_real_)
    tmp_monthly <- sapply(sprintf("tmp_dc_s%02d", 1:12), function(b) props[[b]] %||% NA_real_) / 10

    annual_P  <- if (!is.na(pre) && pre > 0) pre else sum(pre_monthly, na.rm = TRUE)
    F_index   <- if (annual_P > 0 && !all(is.na(pre_monthly)))
      sum(pre_monthly^2, na.rm = TRUE) / annual_P else NA_real_
    R_factor  <- if (is.na(F_index)) NA_real_ else if (F_index < 55) 0.7397 * F_index^1.847 else
      95.77 - 6.081 * F_index + 0.4770 * F_index^2

    list(
      pre_mm      = pre,
      run_mm      = run,
      aet_mm      = aet,
      temp_c      = if (!is.na(tmp_raw)) tmp_raw / 10 else NA_real_,
      tmp_monthly = unname(tmp_monthly),
      slp_dg      = (props$slp_dg_sav %||% NA_real_) / 10,
      cly_pc      = props$cly_pc_sav %||% NA_real_,
      snd_pc      = props$snd_pc_sav %||% NA_real_,
      slt_pc      = props$slt_pc_sav %||% NA_real_,
      ero_kh      = props$ero_kh_sav %||% NA_real_,
      gwt_cm      = props$gwt_cm_sav %||% NA_real_,
      cmi_yr      = (props$cmi_ix_syr %||% NA_real_) / 100,
      ari_ix      = props$ari_ix_sav %||% NA_real_,
      R_factor    = R_factor,
      runoff_frac = if (!is.na(pre) && pre > 0) pmax(0, pmin(0.95, run / pre)) else NA_real_,
      infil_frac  = if (!is.na(pre) && pre > 0) pmax(0, 1 - run / pre) else NA_real_
    )
  }

  fetch_soilgrids_fn <- function(lat, lon) {
    url <- sprintf(
      paste0("https://rest.isric.org/soilgrids/v2.0/properties/query",
             "?lon=%f&lat=%f",
             "&property=phh2o&property=soc&property=clay&property=sand",
             "&property=silt&property=bdod&property=cec",
             "&depth=15-30cm",
             "&value=mean"),
      lon, lat
    )
    last_err <- "No response received"
    for (attempt in seq_len(3)) {
      if (attempt > 1) Sys.sleep(2)
      tryCatch({
        resp <- httr::GET(url, httr::timeout(45))
        sc   <- httr::status_code(resp)
        if (sc == 200) return(httr::content(resp, as = "parsed"))
        last_err <- sprintf("HTTP %d%s", sc,
                            if (sc == 503) ": ISRIC server busy" else "")
      }, error = function(e) {
        last_err <<- conditionMessage(e)
      })
    }
    stop(paste0("SoilGrids unavailable after 3 attempts: ", last_err))
  }

  usda_texture <- function(clay, silt, sand) {
    if (any(is.na(c(clay, silt, sand)))) return(NA_character_)
    if (clay >= 40 && sand > 45)                              return("Sandy clay")
    if (clay >= 40 && silt > 40)                              return("Silty clay")
    if (clay >= 40)                                           return("Clay")
    if (clay >= 35 && sand > 45)                              return("Sandy clay")
    if (clay >= 27 && clay < 40 && sand < 20)                 return("Silty clay loam")
    if (clay >= 27 && clay < 40 && sand >= 45)                return("Sandy clay loam")
    if (clay >= 27 && clay < 40)                              return("Clay loam")
    if (clay >= 20 && clay < 35 && sand >= 45 && silt < 28)   return("Sandy clay loam")
    if (silt >= 80 && clay < 12)                              return("Silt")
    if (silt >= 50)                                           return("Silt loam")
    if (clay >= 7 && clay < 27 && silt >= 28 && sand < 52)    return("Loam")
    if (sand >= 85 && silt + 1.5 * clay < 15)                 return("Sand")
    if (sand >= 70 && silt + 2 * clay < 30)                   return("Loamy sand")
    if ((sand >= 52 || sand >= 43) && clay < 7 && silt < 50)  return("Sandy loam")
    return("Loam")
  }

  parse_soilgrids <- function(result) {
    layers <- result$properties$layers
    extract_val <- function(name) {
      layer <- Filter(function(l) l$name == name, layers)
      if (length(layer) == 0) return(NA_real_)
      dfactor <- layer[[1]]$unit_measure$d_factor %||% 1
      vals <- sapply(layer[[1]]$depths, function(d) {
        v <- d$values$mean
        if (is.null(v)) NA_real_ else as.numeric(v)
      })
      mean(vals, na.rm = TRUE) / dfactor
    }
    ph_val   <- extract_val("phh2o")
    soc_gkg  <- extract_val("soc")
    clay_raw <- extract_val("clay")
    sand_raw <- extract_val("sand")
    silt_raw <- extract_val("silt")
    bd_gcm3  <- extract_val("bdod")
    cec_val  <- extract_val("cec")
    soc_pct  <- soc_gkg / 10
    om_pct   <- soc_pct * 1.724
    clay_pct <- clay_raw
    sand_pct <- sand_raw
    silt_pct <- if (!is.na(silt_raw)) {
      silt_raw
    } else if (!is.na(clay_pct) && !is.na(sand_pct)) {
      pmax(0, 100 - clay_pct - sand_pct)
    } else NA_real_
    texture_class <- usda_texture(clay_pct, silt_pct, sand_pct)
    tex_sum <- sum(c(clay_pct, silt_pct, sand_pct), na.rm = FALSE)
    texture_quality <- if (is.na(tex_sum)) "unknown" else
      if (tex_sum < 50)   "poor"    else
        if (tex_sum < 80 || tex_sum > 110) "fair" else "good"
    om_class <- if (!is.na(om_pct)) {
      if (om_pct < 1) "low" else if (om_pct <= 3) "moderate" else "high"
    } else "moderate"
    list(ph              = ph_val,
         soc_pct         = soc_pct,
         om_pct          = om_pct,
         om_class        = om_class,
         clay_pct        = clay_pct,
         sand_pct        = sand_pct,
         silt_pct        = silt_pct,
         bd_gcm3         = bd_gcm3,
         cec_cmol        = cec_val,
         texture_class   = texture_class,
         texture_sum     = tex_sum,
         texture_quality = texture_quality)
  }

  is_australia <- function(lat, lon) {
    lat >= -44 && lat <= -10 && lon >= 113 && lon <= 154
  }

  # ---- Earth Engine helpers ----
  gee_point_value <- function(image_id, bands, lat, lon, scale, key_file,
                              token = NULL, is_collection = FALSE,
                              filter_prop = NULL, filter_val = NULL,
                              date_start = NULL, date_end = NULL, geom = fetch_geom) {
    if (is.null(token)) token <- gee_token_fn(key_file)
    vl <- list()
    vl[["0"]] <- list(functionInvocationValue = list(
      functionName = "Image.reduceRegion",
      arguments = list(
        image      = list(valueReference = "1"),
        reducer    = list(valueReference = "2"),
        geometry   = list(valueReference = "3"),
        scale      = list(constantValue  = scale),
        bestEffort = list(constantValue  = TRUE)
      )
    ))
    if (is_collection && !is.null(date_start) && !is.null(date_end)) {
      vl[["1"]] <- list(functionInvocationValue = list(
        functionName = "reduce.median",
        arguments    = list(collection = list(valueReference = "5"))
      ))
      vl[["5"]] <- list(functionInvocationValue = list(
        functionName = "Collection.filter",
        arguments    = list(
          collection = list(valueReference = "4"),
          filter     = list(valueReference = "6")
        )
      ))
      vl[["4"]] <- list(functionInvocationValue = list(
        functionName = "ImageCollection.load",
        arguments    = list(id = list(constantValue = image_id))
      ))
      vl[["6"]] <- list(functionInvocationValue = list(
        functionName = "Filter.dateRangeContains",
        arguments    = list(
          leftValue  = list(functionInvocationValue = list(functionName = "DateRange", arguments = list(
            start = list(constantValue = date_start), end = list(constantValue = date_end)))),
          rightField = list(constantValue = "system:time_start")
        )
      ))
    } else if (is_collection && !is.null(filter_prop)) {
      vl[["1"]] <- list(functionInvocationValue = list(
        functionName = "Collection.first",
        arguments    = list(collection = list(valueReference = "5"))
      ))
      vl[["5"]] <- list(functionInvocationValue = list(
        functionName = "Collection.filter",
        arguments    = list(
          collection = list(valueReference = "4"),
          filter     = list(valueReference = "6")
        )
      ))
      vl[["4"]] <- list(functionInvocationValue = list(
        functionName = "ImageCollection.load",
        arguments    = list(id = list(constantValue = image_id))
      ))
      vl[["6"]] <- list(functionInvocationValue = list(
        functionName = "Filter.equals",
        arguments    = list(
          leftField  = list(constantValue = filter_prop),
          rightValue = list(constantValue = filter_val)
        )
      ))
    } else if (is_collection) {
      vl[["1"]] <- list(functionInvocationValue = list(
        functionName = "Collection.first",
        arguments    = list(collection = list(valueReference = "4"))
      ))
      vl[["4"]] <- list(functionInvocationValue = list(
        functionName = "ImageCollection.load",
        arguments    = list(id = list(constantValue = image_id))
      ))
    } else {
      vl[["1"]] <- list(functionInvocationValue = list(
        functionName = "Image.load",
        arguments    = list(id = list(constantValue = image_id))
      ))
    }
    vl[["2"]] <- list(functionInvocationValue = list(
      functionName = "Reducer.mean"
    ))
    vl[["3"]] <- if (!is.null(geom)) list(functionInvocationValue = list(
      functionName = "GeometryConstructors.Polygon",
      arguments    = list(coordinates = list(constantValue = list(geom)))
    )) else list(functionInvocationValue = list(
      functionName = "GeometryConstructors.Point",
      arguments    = list(coordinates = list(constantValue = c(lon, lat)))
    ))
    body_json <- jsonlite::toJSON(
      list(expression = list(result = "0", values = vl)),
      auto_unbox = TRUE, digits = NA
    )
    resp <- httr::POST(
      "https://earthengine.googleapis.com/v1/projects/earthengine-public/value:compute",
      httr::add_headers(Authorization  = paste("Bearer", token),
                        `Content-Type` = "application/json"),
      body = body_json, encode = "raw", httr::timeout(30)
    )
    if (httr::http_error(resp))
      stop(sprintf("GEE HTTP %d: %s", httr::status_code(resp),
                   jsonlite::toJSON(httr::content(resp))))
    result <- httr::content(resp)
    if (!is.null(result$error))
      stop(paste("GEE value:compute:", result$error$message %||% "unknown error"))
    vals <- result$result
    if (!is.null(geom) && (is.null(vals) || all(vapply(bands, function(b) is.null(vals[[b]]), TRUE))))
      return(gee_point_value(image_id, bands, lat, lon, scale, key_file, token, is_collection, filter_prop, filter_val,
                             date_start, date_end, geom = NULL))
    if (is.null(vals)) stop("GEE returned no pixel values for this location")
    out <- list()
    for (b in bands) out[[b]] <- if (!is.null(vals[[b]])) as.numeric(vals[[b]]) else NA_real_
    out
  }

  gee_slope_deg <- function(lat, lon, token) {
    fi  <- function(name, args) list(functionInvocationValue = list(functionName = name, arguments = args))
    cst <- function(x) list(constantValue = x)
    img <- fi("Terrain.slope", list(input = fi("Image.select", list(
      input = fi("Image.load", list(id = cst("MERIT/DEM/v1_0_3"))), bandSelectors = cst(list("dem"))))))
    vl <- list("0" = fi("Image.reduceRegion", list(image = img,
      reducer = list(functionInvocationValue = list(functionName = "Reducer.mean")),
      geometry = fi("GeometryConstructors.Point", list(coordinates = cst(c(lon, lat)))),
      scale = cst(90), bestEffort = cst(TRUE))))
    body_json <- jsonlite::toJSON(list(expression = list(result = "0", values = vl)), auto_unbox = TRUE, digits = NA)
    resp <- httr::POST("https://earthengine.googleapis.com/v1/projects/earthengine-public/value:compute",
                       httr::add_headers(Authorization = paste("Bearer", token), `Content-Type` = "application/json"),
                       body = body_json, encode = "raw", httr::timeout(30))
    if (httr::http_error(resp)) stop(sprintf("GEE HTTP %d", httr::status_code(resp)))
    v <- httr::content(resp)$result$slope
    if (is.null(v)) NA_real_ else as.numeric(v)
  }

  # ---- Daily rainfall ----
  RAIN_DATASETS <- list(
    chirps = list(id = "UCSB-CHG/CHIRPS/DAILY", band = "precipitation", scale = 5566, mult = 1,
                  label = "CHIRPS daily rainfall (Funk et al. 2015)", name = "CHIRPS", cite = "Funk et al. 2015", ref = "funk"),
    era5   = list(id = "ECMWF/ERA5_LAND/DAILY_AGGR", band = "total_precipitation_sum", scale = 11132, mult = 1000,
                  label = "ERA5-Land daily rainfall (Mu\u00f1oz-Sabater et al. 2021)", name = "ERA5-Land",
                  cite = "Mu\u00f1oz-Sabater et al. 2021", ref = "munoz"))
  rain_source_info <- function(rain) {
    src <- attr(rain, "source") %||% ""
    if (grepl("ERA5", src)) RAIN_DATASETS$era5 else RAIN_DATASETS$chirps
  }
  gee_daily_rain <- function(lat, lon, key_file, dataset = "chirps", start = "1985-01-01", end = "2025-01-01",
                             token = NULL) {
    ds <- RAIN_DATASETS[[dataset]]
    if (is.null(token)) token <- gee_token_fn(key_file)
    fi  <- function(name, args) list(functionInvocationValue = list(functionName = name, arguments = args))
    cst <- function(x) list(constantValue = x)
    coll <- fi("Collection.map", list(
      baseAlgorithm = list(functionDefinitionValue = list(argumentNames = list("img"), body = "1")),
      collection = fi("Collection.filter", list(
        collection = fi("ImageCollection.load", list(id = cst(ds$id))),
        filter = fi("Filter.dateRangeContains", list(
          leftValue  = fi("DateRange", list(start = cst(start), end = cst(end))),
          rightField = cst("system:time_start")))))))
    vl <- list(
      "0" = fi("ImageCollection.getRegion", list(collection = coll,
              geometry = fi("GeometryConstructors.Point", list(coordinates = cst(c(lon, lat)))),
              scale = cst(ds$scale))),
      "1" = fi("Image.select", list(input = list(argumentReference = "img"), bandSelectors = cst(list(ds$band)))))
    body_json <- jsonlite::toJSON(list(expression = list(result = "0", values = vl)), auto_unbox = TRUE, digits = NA)
    resp <- httr::POST(
      "https://earthengine.googleapis.com/v1/projects/earthengine-public/value:compute",
      httr::add_headers(Authorization = paste("Bearer", token), `Content-Type` = "application/json"),
      body = body_json, encode = "raw", httr::timeout(120))
    if (httr::http_error(resp))
      stop(sprintf("GEE HTTP %d: %s", httr::status_code(resp), substr(httr::content(resp, "text"), 1, 200)))
    rows <- httr::content(resp)$result
    if (length(rows) < 2) return(NULL)
    hdr <- unlist(rows[[1]]); body <- rows[-1]
    ti <- match("time", hdr); bi <- match(ds$band, hdr)
    tms <- vapply(body, function(r) as.numeric(r[[ti]]), 0)
    val <- vapply(body, function(r) if (is.null(r[[bi]])) NA_real_ else as.numeric(r[[bi]]), 0)
    data.frame(date = as.Date(as.POSIXct(tms / 1000, origin = "1970-01-01", tz = "UTC")), P_mm = val * ds$mult)
  }
  fetch_daily_rain_fn <- function(lat, lon, key_file, years = RAIN_YEARS) {
    token <- gee_token_fn(key_file)
    out <- NULL; src <- NULL
    if (abs(lat) <= 50) {
      out <- tryCatch(gee_daily_rain(lat, lon, key_file, "chirps", sprintf("%d-01-01", min(years)),
                                     sprintf("%d-01-01", max(years) + 1), token), error = function(e) NULL)
      if (!is.null(out) && sum(!is.na(out$P_mm)) > 3650) src <- RAIN_DATASETS$chirps$label else out <- NULL
    }
    if (is.null(out)) {
      starts <- seq(min(years), max(years), by = 10)
      out <- do.call(rbind, lapply(starts, function(y)
        gee_daily_rain(lat, lon, key_file, "era5", sprintf("%d-01-01", y),
                       sprintf("%d-01-01", min(y + 10, max(years) + 1)), token)))
      src <- RAIN_DATASETS$era5$label
    }
    if (is.null(out) || nrow(out) < 365) stop("no daily rainfall returned for this location")
    out <- out[order(out$date), ]
    attr(out, "source") <- src
    out
  }

  # ---- Soil: CSIRO SLGA ----
  fetch_slga_fn <- function(lat, lon, key_file) {
    token <- gee_token_fn(key_file)
    props <- c(CLY = "CLY", SLT = "SLT", SND = "SND",
               BDW = "BDW", SOC = "SOC", pHc = "pHc", ECE = "ECE")
    out <- list()
    for (code in names(props)) {
      band <- paste0(code, "_015_030_EV")
      out[[band]] <- tryCatch(
        gee_point_value("CSIRO/SLGA", band, lat, lon,
                        scale = 90, key_file = key_file, token = token,
                        is_collection = TRUE, filter_prop = "attribute_code",
                        filter_val = code)[[band]],
        error = function(e) NA_real_
      )
    }
    out
  }

  parse_slga <- function(result) {
    gb <- function(band) {
      v <- result[[band]]
      if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
    }
    if (all(vapply(result, function(v) is.null(v) || is.na(v), TRUE))) stop("SLGA returned no values here")
    clay_pct <- gb("CLY_015_030_EV")
    silt_pct <- gb("SLT_015_030_EV")
    sand_pct <- gb("SND_015_030_EV")
    bd_gcm3  <- gb("BDW_015_030_EV")
    soc_pct  <- gb("SOC_015_030_EV")
    ph_val   <- gb("pHc_015_030_EV")
    if (is.na(silt_pct) && !is.na(clay_pct) && !is.na(sand_pct))
      silt_pct <- pmax(0, 100 - clay_pct - sand_pct)
    tex_sum <- sum(c(clay_pct, silt_pct, sand_pct), na.rm = FALSE)
    tq      <- if (is.na(tex_sum)) "unknown" else
      if (tex_sum < 50) "poor" else if (tex_sum < 80 || tex_sum > 110) "fair" else "good"
    om_pct  <- if (!is.na(soc_pct)) soc_pct * 1.724 else NA_real_
    om_class <- if (!is.na(om_pct)) {
      if (om_pct < 1) "low" else if (om_pct <= 3) "moderate" else "high"
    } else "moderate"
    list(ph = ph_val, soc_pct = soc_pct, om_pct = om_pct, om_class = om_class,
         clay_pct = clay_pct, sand_pct = sand_pct, silt_pct = silt_pct,
         bd_gcm3 = bd_gcm3, cec_cmol = NA_real_,
         texture_class   = usda_texture(clay_pct, silt_pct, sand_pct),
         texture_sum     = tex_sum,
         texture_quality = tq)
  }

  # ---- Soil: OpenLandMap ----
  fetch_olm_fn <- function(lat, lon, key_file) {
    token <- gee_token_fn(key_file)
    gv <- function(id) gee_point_value(id, c("b30"), lat, lon, 250, key_file, token)
    list(
      clay = gv("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02"),
      sand = gv("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02"),
      ph   = gv("OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02"),
      oc   = gv("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02"),
      bd   = gv("OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02")
    )
  }

  parse_olm <- function(result) {
    get_b30 <- function(lst) {
      v <- lst[["b30"]]
      if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
    }
    if (all(is.na(vapply(result, get_b30, 0)))) stop("OpenLandMap returned no values here")
    clay_raw <- get_b30(result$clay)
    sand_raw <- get_b30(result$sand)
    ph_raw   <- get_b30(result$ph)
    oc_raw   <- get_b30(result$oc)
    bd_raw   <- get_b30(result$bd)
    clay_pct <- clay_raw
    sand_pct <- sand_raw
    silt_pct <- if (!is.na(clay_pct) && !is.na(sand_pct))
      pmax(0, 100 - clay_pct - sand_pct) else NA_real_
    ph_val   <- if (!is.na(ph_raw))  ph_raw / 10            else NA_real_
    soc_gkg  <- if (!is.na(oc_raw))  oc_raw * 5             else NA_real_
    soc_pct  <- if (!is.na(soc_gkg)) soc_gkg / 10           else NA_real_
    bd_gcm3  <- if (!is.na(bd_raw))  bd_raw * 10 / 1000     else NA_real_
    tex_sum  <- sum(c(clay_pct, silt_pct, sand_pct), na.rm = FALSE)
    tq       <- if (is.na(tex_sum)) "unknown" else
      if (tex_sum < 50) "poor" else if (tex_sum < 80 || tex_sum > 110) "fair" else "good"
    om_pct   <- if (!is.na(soc_pct)) soc_pct * 1.724 else NA_real_
    om_class <- if (!is.na(om_pct)) {
      if (om_pct < 1) "low" else if (om_pct <= 3) "moderate" else "high"
    } else "moderate"
    list(ph = ph_val, soc_pct = soc_pct, om_pct = om_pct, om_class = om_class,
         clay_pct = clay_pct, sand_pct = sand_pct, silt_pct = silt_pct,
         bd_gcm3 = bd_gcm3, cec_cmol = NA_real_,
         texture_class   = usda_texture(clay_pct, silt_pct, sand_pct),
         texture_sum     = tex_sum,
         texture_quality = tq)
  }

  add_field_draw <- function(m) {
    if (!HAS_DRAW) return(m)
    m <- leaflet.extras::addDrawToolbar(m, targetGroup = "field", singleFeature = TRUE,
           polylineOptions = FALSE, circleOptions = FALSE, markerOptions = FALSE, circleMarkerOptions = FALSE,
           polygonOptions = leaflet.extras::drawPolygonOptions(
             shapeOptions = leaflet.extras::drawShapeOptions(color = "#ffeb3b", fillOpacity = 0.15, weight = 2)),
           rectangleOptions = leaflet.extras::drawRectangleOptions(
             shapeOptions = leaflet.extras::drawShapeOptions(color = "#ffeb3b", fillOpacity = 0.15, weight = 2)),
           editOptions = leaflet.extras::editToolbarOptions())
    htmlwidgets::onRender(m, "function(el, x) { var map = this;
      var on = function() { Shiny.setInputValue('field_drawing', true); };
      var off = function() { setTimeout(function() { Shiny.setInputValue('field_drawing', false); }, 300); };
      map.on('draw:drawstart draw:editstart draw:deletestart', on);
      map.on('draw:drawstop draw:editstop draw:deletestop', off); }")
  }

  output$location_map <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery) %>%
      addTiles(urlTemplate = paste0("https://server.arcgisonline.com/ArcGIS/rest/services/Reference/",
                                    "World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}"),
               attribution = "Labels &copy; Esri, HERE, Garmin, &copy; OpenStreetMap contributors",
               options = tileOptions(opacity = 0.9)) %>%
      setView(lng = 134, lat = -26, zoom = 3) %>%
      add_field_draw()
  })

  soilgrids_status_msg  <- reactiveVal(list(type = "none", msg = ""))
  hydroatlas_status_msg <- reactiveVal(list(type = "none", msg = ""))

  run_soil_fetch <- function(lat, lon, source_pref = "auto") {
    au     <- is_australia(lat, lon)
    key_ok <- file.exists(gee_key_path)

    src_nice <- function(s) switch(s,
                                   soilgrids = "SoilGrids (ISRIC)",
                                   slga      = "CSIRO SLGA (90 m)",
                                   olm       = "OpenLandMap (250 m)",
                                   s
    )
    src_full <- function(s) switch(s,
                                   soilgrids = "SoilGrids (ISRIC, 250 m)",
                                   slga      = "CSIRO SLGA (Australia, 90 m)",
                                   olm       = "OpenLandMap (Global, 250 m)",
                                   s
    )

    pref_warn <- NULL
    order <- if (source_pref == "soilgrids") {
      "soilgrids"
    } else if (source_pref == "slga") {
      if (!key_ok) {
        pref_warn <- "SLGA selected but GEE key file not found - falling back to SoilGrids."
        c("soilgrids", "olm")
      } else if (!au) {
        pref_warn <- paste0("SLGA selected but this location appears to be outside Australia ",
                            "(bounds 10-44 deg S, 113-154 deg E) - ",
                            "falling back to SoilGrids.")
        c("soilgrids", "olm")
      } else {
        c("slga", "soilgrids", "olm")
      }
    } else if (source_pref == "olm") {
      if (!key_ok) {
        pref_warn <- "OpenLandMap selected but GEE key file not found - falling back to SoilGrids."
        "soilgrids"
      } else {
        c("olm", "soilgrids")
      }
    } else {
      if (au && key_ok) c("slga", "soilgrids", "olm") else
        if (key_ok)       c("soilgrids", "olm")          else "soilgrids"
    }

    if (!is.null(pref_warn))
      showNotification(pref_warn, type = "warning", duration = 8)

    n_src      <- length(order)
    all_errors <- character(0)

    for (i in seq_along(order)) {
      src   <- order[[i]]
      label <- src_nice(src)

      notif_id <- "soil_fetch_prog"
      showNotification(
        paste0("Soil: trying ", label,
               if (n_src > 1) paste0(" (", i, " of ", n_src, ")") else "",
               "..."),
        id = notif_id, duration = NULL, type = "message"
      )

      attempt <- tryCatch({
        sg <- switch(src,
                     soilgrids = parse_soilgrids(fetch_soilgrids_fn(lat, lon)),
                     slga      = parse_slga(fetch_slga_fn(lat, lon, gee_key_path)),
                     olm       = parse_olm(fetch_olm_fn(lat, lon, gee_key_path))
        )
        list(ok = TRUE, data = sg, src = src)
      }, error = function(e) {
        list(ok = FALSE, err = conditionMessage(e), src = src)
      })

      removeNotification(notif_id)

      if (isTRUE(attempt$ok)) {
        sg    <- attempt$data
        label_full <- src_full(attempt$src)
        soilgrids_data(sg)
        soil_source_name(label_full)
        soil_source_id(attempt$src)
        sel_om(sg$om_class)
        if (!is.na(sg$ph %||% NA_real_)) sel_ph(sg$ph)
        soilgrids_status_msg(list(
          type = "success",
          msg  = sprintf(
            "Soil loaded from %s: pH\u00a0%.1f, OC\u00a0%.1f%%, Clay\u00a0%.0f%%, Silt\u00a0%.0f%%",
            label_full,
            sg$ph       %||% NA_real_,
            sg$soc_pct  %||% NA_real_,
            sg$clay_pct %||% NA_real_,
            sg$silt_pct %||% NA_real_)
        ))
        return(invisible(NULL))
      }

      err_msg <- attempt$err
      all_errors <- c(all_errors, sprintf("%s: %s", label, err_msg))
      showNotification(paste0(label, " failed - ",
                              if (i < n_src) "trying next source..." else "no more sources."),
                       type = "warning", duration = 5)
    }

    tried_str <- paste(sapply(order, src_nice), collapse = " \u2192 ")
    soilgrids_status_msg(list(
      type = "error",
      msg  = paste0(
        "Soil data unavailable. Tried: ", tried_str, ". ",
        "Errors: ", paste(all_errors, collapse = "; "), ". ",
        "Calculations will use default values. ",
        "Try clicking the map again, change the source preference, or use the Retry button."
      )
    ))
  }

  # ---- RUSLE sediment concentration ----
  fetch_rusle_fn <- function(lat, lon, key_file) {
    token <- gee_token_fn(key_file)

    R_factor <- NA_real_

    merit_result <- gee_point_value("MERIT/Hydro/v1_0_1", "upa",
                                    lat, lon, scale = 90, key_file = key_file, token = token, geom = NULL)
    upa_km2   <- merit_result[["upa"]] %||% NA_real_
    slp_deg_m <- tryCatch(gee_slope_deg(lat, lon, token), error = function(e) NA_real_)
    theta_rad <- if (!is.na(slp_deg_m)) slp_deg_m * pi / 180 else NA_real_

    LS_factor <- if (is.na(upa_km2) || is.na(theta_rad)) NA_real_ else if (theta_rad <= 0) 0 else {
      lambda_m <- min(sqrt(upa_km2 * 1e6), 500)
      beta_ls  <- (sin(theta_rad) / 0.0896) / (3 * sin(theta_rad)^0.8 + 0.56)
      m_exp    <- beta_ls / (1 + beta_ls)
      (m_exp + 1) * (lambda_m / 22.13)^m_exp * (sin(theta_rad) / 0.0896)^1.3
    }

    ndvi_end_d <- if (is.na(NDVI_END_DATE)) Sys.Date() else as.Date(NDVI_END_DATE)
    ndvi_end   <- format(ndvi_end_d, "%Y-%m-%d")
    ndvi_start <- format(ndvi_end_d - NDVI_WINDOW_DAYS, "%Y-%m-%d")
    ndvi_result <- tryCatch({
      gee_point_value("LANDSAT/COMPOSITES/C02/T1_L2_32DAY_NDVI", "NDVI",
                      lat, lon, scale = 30, key_file = key_file, token = token,
                      is_collection = TRUE, date_start = ndvi_start, date_end = ndvi_end)
    }, error = function(e) list(NDVI = NA_real_))
    ndvi_raw <- ndvi_result[["NDVI"]] %||% NA_real_
    ndvi_val <- ndvi_raw

    C_factor  <- if (!is.na(ndvi_val) && ndvi_val > 0 && ndvi_val < 1) {
      pmax(0.001, pmin(0.999, exp(-2 * ndvi_val / (1 - ndvi_val))))
    } else 0.264

    list(R_factor = R_factor, LS_factor = LS_factor, C_factor = C_factor,
         upa_km2 = upa_km2, slp_deg = if (!is.na(theta_rad)) theta_rad * 180 / pi else NA_real_,
         ndvi = ndvi_val, ndvi_window = paste(ndvi_start, "to", ndvi_end))
  }

  gee_key_path <- Sys.getenv("PESTICIDE_GEE_KEY", "gee-service-account.json")

  # ---- Field location ----
  field_geom <- reactiveVal(NULL)
  fetch_geom <- NULL
  last_draw  <- reactiveVal(NULL)
  ring_area_ha <- function(ring) {
    lat0 <- mean(vapply(ring, `[`, 0, 2)) * pi / 180
    x <- vapply(ring, `[`, 0, 1) * 111320 * cos(lat0); y <- vapply(ring, `[`, 0, 2) * 110540
    abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2 / 1e4
  }
  ring_centroid <- function(ring) {
    r <- if (identical(ring[[1]], ring[[length(ring)]])) ring[-length(ring)] else ring
    c(lon = mean(vapply(r, `[`, 0, 1)), lat = mean(vapply(r, `[`, 0, 2)))
  }
  feature_ring <- function(f) {
    cc <- f$geometry$coordinates[[1]]
    ring <- lapply(cc, function(p) as.numeric(unlist(p))[1:2])
    if (!identical(ring[[1]], ring[[length(ring)]])) ring <- c(ring, ring[1])
    ring
  }

  fetch_site <- function(lat, lon, ring = NULL) {
    map_coords(list(lat = lat, lon = lon))
    field_geom(if (!is.null(ring)) list(ring = ring, area_ha = ring_area_ha(ring)) else NULL)
    fetch_geom <<- ring
    on.exit(fetch_geom <<- NULL, add = TRUE)
    user_runoff_frac(NULL); user_slope_dg(NULL); user_temp_c(NULL); user_cn(NULL)
    user_clay_pct(NULL);    user_sand_pct(NULL); user_silt_pct(NULL)
    user_bd_gcm3(NULL);     user_soc_pct(NULL)
    src_pref  <- input$soil_source_pref %||% "auto"
    src_label <- switch(src_pref,
                        auto      = "Auto (SLGA \u2192 SoilGrids \u2192 OLM)",
                        soilgrids = "SoilGrids (ISRIC)",
                        slga      = "CSIRO SLGA",
                        olm       = "OpenLandMap",
                        "Auto"
    )
    shinyjs::runjs(sprintf(
      "document.getElementById('location_fetch_overlay').style.display='block';
       var el = document.getElementById('fetch_overlay_src');
       if(el) el.innerText = 'Soil source: %s';",
      gsub("'", "\\\\'", src_label)
    ))
    soilgrids_data(NULL)
    hydroatlas_data(NULL)
    rain_daily_data(NULL)
    rusle_data(NULL)
    soilgrids_status_msg(list(type = "none", msg = ""))
    hydroatlas_status_msg(list(type = "none", msg = ""))

    # ---- Soil ----
    run_soil_fetch(lat, lon, src_pref)

    # ---- HydroATLAS ----
    if (!file.exists(gee_key_path)) {
      hydroatlas_status_msg(list(type = "error", msg = "GEE key file not found. Water balance using defaults."))
      shinyjs::runjs("document.getElementById('location_fetch_overlay').style.display='none';")
      return()
    }
    tryCatch({
      result <- fetch_hydroatlas_fn(lat, lon, gee_key_path)
      if (is.null(result$features) || length(result$features) == 0) {
        hydroatlas_status_msg(list(type = "error", msg = "No basin found here. Try a nearby point on a river or land area."))
      } else {
        parsed <- parse_hydroatlas(result)
        hydroatlas_data(parsed)
        hydroatlas_status_msg(list(type = "success",
                                   msg = sprintf("Water balance loaded: %.0f%% runoff, MAP %.0f mm, %.1f\u00b0C",
                                                 parsed$runoff_frac * 100, parsed$pre_mm %||% NA, parsed$temp_c %||% NA)))
      }
    }, error = function(e) {
      hydroatlas_status_msg(list(type = "error", msg = paste("Water balance unavailable:", conditionMessage(e))))
    })

    # ---- Daily rainfall ----
    tryCatch({
      rain <- fetch_daily_rain_fn(lat, lon, gee_key_path)
      rain_daily_data(rain)
    }, error = function(e) {
      rain_daily_data(NULL)
      showNotification(sprintf("Daily rainfall unavailable, so %d days to first runoff is assumed: %s",
                               as.integer(DEFAULT_DAYS_TO_RUNOFF), conditionMessage(e)), type = "warning", duration = 6)
    })

    # ---- RUSLE factors ----
    tryCatch({
      rusle_data(fetch_rusle_fn(lat, lon, gee_key_path))
    }, error = function(e) {
      rusle_data(NULL)
      showNotification(paste("Slope and ground cover unavailable, so a default sediment value is used:", conditionMessage(e)),
                       type = "warning", duration = 5)
    })

    shinyjs::runjs("document.getElementById('location_fetch_overlay').style.display='none';")
  }

  # ---- Location page ----
  observeEvent(input$location_map_click, {
    if (isTRUE(input$field_drawing)) return()
    if (!is.null(last_draw()) && difftime(Sys.time(), last_draw(), units = "secs") < 2) return()
    click <- input$location_map_click
    lat <- round(click$lat, 6)
    lon <- round(click$lng, 6)
    leafletProxy("location_map") %>%
      clearGroup("field") %>%
      clearMarkers() %>%
      addMarkers(lng = lon, lat = lat,
                 popup = sprintf("%.5f\u00b0, %.5f\u00b0", lat, lon))
    fetch_site(lat, lon)
  })
  draw_field <- function(f) {
    if (is.null(f) || is.null(f$geometry)) return()
    ring <- tryCatch(feature_ring(f), error = function(e) NULL)
    if (is.null(ring) || length(ring) < 4) return()
    last_draw(Sys.time())
    cen <- ring_centroid(ring)
    leafletProxy("location_map") %>% clearMarkers()
    fetch_site(round(cen[["lat"]], 6), round(cen[["lon"]], 6), ring)
  }
  observeEvent(input$location_map_draw_new_feature, draw_field(input$location_map_draw_new_feature))
  observeEvent(input$location_map_draw_edited_features, {
    fs <- input$location_map_draw_edited_features$features
    if (length(fs)) draw_field(fs[[1]])
  })
  observeEvent(input$location_map_draw_deleted_features, field_geom(NULL))

  observeEvent(input$retry_soil_fetch, {
    coords <- map_coords()
    if (is.null(coords)) {
      showNotification("Click the map first to set a location.", type = "warning", duration = 3)
      return()
    }
    src_pref  <- input$soil_source_pref %||% "auto"
    src_label <- switch(src_pref,
                        auto      = "Auto (SLGA \u2192 SoilGrids \u2192 OLM)",
                        soilgrids = "SoilGrids (ISRIC)",
                        slga      = "CSIRO SLGA",
                        olm       = "OpenLandMap",
                        "Auto"
    )
    soilgrids_data(NULL)
    soil_source_name(NULL)
    soilgrids_status_msg(list(type = "info",
                              msg = paste0("Retrying soil fetch from ", src_label, "...")))
    run_soil_fetch(coords$lat, coords$lon, src_pref)
  })

  output$fetch_status_panel <- renderUI({
    sg_s <- soilgrids_status_msg()
    ha_s <- hydroatlas_status_msg()
    if (sg_s$type == "none" && ha_s$type == "none") return(NULL)
    status_row <- function(s) {
      if (s$type == "none") return(NULL)
      col <- switch(s$type, "error"="#c0392b", "success"="#027a55", "info"="#0173b2")
      bg  <- switch(s$type, "error"="#fdecea", "success"="#e8faf3", "info"="#e8f4ff")
      div(style = sprintf("background:%s; border-left:3px solid %s; padding:7px 10px;
                            border-radius:3px; font-size:12px; color:%s; margin-bottom:6px;",
                          bg, col, col), s$msg)
    }
    tagList(if (sg_s$type == "error") status_row(sg_s), if (ha_s$type == "error") status_row(ha_s))
  })

  output$hydroatlas_results_panel <- renderUI({
    ha <- hydroatlas_data()
    if (is.null(ha)) return(NULL)
    rf   <- ha$runoff_frac
    inff <- ha$infil_frac
    tr <- function(lab, val, tip) tags$tr(tags$td(strong(lab), info_tip(tip), style = "padding:3px;"),
                                          tags$td(val, style = "padding:3px;"))
    big <- function(val, col, lab, tip) div(style = "text-align:center; padding:10px; background:#f0fdf9; border-radius:4px;",
                                            div(style = sprintf("font-size:28px; font-weight:bold; color:%s;", col), val),
                                            div(lab, info_tip(tip), style = "font-size:12px; color:#555;"))
    div(style = "margin-top: 16px; background: white; border: 1px solid #b2dfdb;
                  border-radius: 6px; padding: 16px;",
        h5("Water balance", style = "margin-top: 0; margin-bottom: 2px; color: #01745f;"),
        tags$p(class = "src-line", "Long-term sub-basin averages sourced from HydroATLAS (Linke et al. 2019)."),
        fluidRow(
          column(4, big(sprintf("%.0f%%", rf * 100), "#cc3311", "Runoff fraction",
                        sprintf("Share of rainfall that leaves the sub-basin as runoff (%.0f of %.0f mm/yr). Used for runoff only when no daily rainfall record is available.",
                                ha$run_mm, ha$pre_mm))),
          column(4, big(sprintf("%.0f%%", inff * 100), "#0173b2", "Infiltration + ET",
                        "Share of rainfall that soaks into the soil or returns to the air by evapotranspiration (100% minus the runoff fraction).")),
          column(4,
                 tags$table(style = "width:100%; font-size:12px;",
                            tr("Rainfall:", sprintf("%.0f mm/yr", ha$pre_mm), "Mean annual precipitation of the sub-basin."),
                            tr("Runoff:", sprintf("%.0f mm/yr", ha$run_mm), "Mean annual runoff from the sub-basin."),
                            tr("Actual ET:", sprintf("%.0f mm/yr", ha$aet_mm), "Mean annual actual evapotranspiration: evaporation from soil and water plus transpiration by plants."),
                            tr("CMI:", sprintf("%.2f", ha$cmi_yr %||% NA_real_), "Climate moisture index, from -1 (very dry) to +1 (very wet)."),
                            tr("Slope:", sprintf("%.1f\u00b0", ha$slp_dg %||% NA_real_), "Mean terrain slope of the sub-basin.")))
        )
    )
  })

  output$overwrite_panel <- renderUI({
    ha <- hydroatlas_data()
    sg <- soilgrids_data()
    if (is.null(ha)) return(NULL)

    rf_default   <- round(ha$runoff_frac, 3)
    temp_default <- if (!is.na(ha$temp_c)) round(ha$temp_c, 1) else 25
    slp_default  <- if (!is.na(ha$slp_dg)) round(ha$slp_dg, 1) else 2
    cly_default  <- if (!is.null(sg) && !is.na(sg$clay_pct)) round(sg$clay_pct) else
                    if (!is.na(ha$cly_pc)) round(ha$cly_pc) else 25L
    snd_default  <- if (!is.null(sg) && !is.na(sg$sand_pct)) round(sg$sand_pct) else
                    if (!is.null(ha) && !is.na(ha$snd_pc)) round(ha$snd_pc) else NA_integer_
    slt_default  <- if (!is.null(sg) && !is.na(sg$silt_pct)) round(sg$silt_pct) else
                    if (!is.null(ha) && !is.na(ha$slt_pc)) round(ha$slt_pc) else NA_integer_
    bd_default   <- if (!is.null(sg) && !is.na(sg$bd_gcm3)) round(sg$bd_gcm3, 2) else NA_real_
    soc_default  <- if (!is.null(sg) && !is.na(sg$soc_pct)) round(sg$soc_pct, 2) else NA_real_
    om_default   <- isolate(sel_om())
    ph_default   <- round(isolate(sel_ph()), 1)

    tagList(
      div(style = "margin-top:14px;",
          actionButton("toggle_overwrite", "Edit retrieved values",
                       class = "btn btn-default btn-sm",
                       style = "font-size:12px; color:#0173b2; border-color:#0173b2;"),
          shinyjs::hidden(
            div(id = "overwrite_inputs_div",
                style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;
                          padding:14px; margin-top:8px;",
                tags$b(style = "font-size:11px; color:#333;", "Water balance"),
                fluidRow(
                  column(6,
                         numericInput("ow_curve_number", "Curve number (30-98)",
                                      value = isolate(site_cn())$CN, min = 30, max = 98, step = 1)),
                  column(6,
                         numericInput("ow_slope_dg", "Slope (degrees)",
                                      value = slp_default, min = 0, max = 90, step = 0.5))
                ),
                tags$b(style = "font-size:11px; color:#333;", "Climate"),
                fluidRow(
                  column(6,
                         numericInput("ow_temp_c", "Temperature (\u00b0C)",
                                      value = temp_default, min = -10, max = 50, step = 0.5))
                ),
                tags$b(style = "font-size:11px; color:#333;", "Soil texture"),
                fluidRow(
                  column(4,
                         numericInput("ow_clay_pct", "Clay (%)",
                                      value = cly_default, min = 0, max = 100, step = 1)),
                  column(4,
                         numericInput("ow_sand_pct", "Sand (%)",
                                      value = snd_default, min = 0, max = 100, step = 1)),
                  column(4,
                         numericInput("ow_silt_pct", "Silt (%)",
                                      value = slt_default, min = 0, max = 100, step = 1))
                ),
                tags$b(style = "font-size:11px; color:#333;", "Soil chemistry & structure"),
                fluidRow(
                  column(6,
                         numericInput("ow_bd_gcm3", "Bulk density (g/cm\u00b3)",
                                      value = bd_default, min = 0.5, max = 2.0, step = 0.05)),
                  column(6,
                         numericInput("ow_soc_pct", "Soil organic carbon (%)",
                                      value = soc_default, min = 0, max = 20, step = 0.1))
                ),
                fluidRow(
                  column(6,
                         selectInput("ow_om", "Organic matter class",
                                     choices = c("Low" = "low", "Moderate" = "moderate", "High" = "high"),
                                     selected = om_default)),
                  column(6,
                         numericInput("ow_ph", "Soil pH",
                                      value = ph_default, min = 3, max = 9, step = 0.1))
                ),
                actionButton("reset_overwrite", "Reset to fetched values",
                             class = "btn btn-warning btn-xs",
                             style = "font-size:11px; margin-top:4px;")
            )
          )
      )
    )
  })

  observeEvent(input$toggle_overwrite, {
    shinyjs::toggle("overwrite_inputs_div")
  })

  observeEvent(input$ow_curve_number, {
    req(!is.na(input$ow_curve_number))
    v <- max(30, min(98, input$ow_curve_number))
    user_cn(if (isTRUE(all.equal(v, isolate(site_cn())$CN_table))) NULL else v)
  }, ignoreInit = TRUE)
  observeEvent(input$ow_slope_dg, {
    req(!is.na(input$ow_slope_dg))
    user_slope_dg(max(0, min(90, input$ow_slope_dg)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_temp_c, {
    req(!is.na(input$ow_temp_c))
    user_temp_c(max(-10, min(50, input$ow_temp_c)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_clay_pct, {
    req(!is.na(input$ow_clay_pct))
    user_clay_pct(max(0, min(100, input$ow_clay_pct)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_sand_pct, {
    req(!is.na(input$ow_sand_pct))
    user_sand_pct(max(0, min(100, input$ow_sand_pct)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_silt_pct, {
    req(!is.na(input$ow_silt_pct))
    user_silt_pct(max(0, min(100, input$ow_silt_pct)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_bd_gcm3, {
    req(!is.na(input$ow_bd_gcm3))
    user_bd_gcm3(max(0.5, min(2.0, input$ow_bd_gcm3)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_soc_pct, {
    req(!is.na(input$ow_soc_pct))
    user_soc_pct(max(0, min(20, input$ow_soc_pct)))
  }, ignoreInit = TRUE)
  observeEvent(input$ow_om, { sel_om(input$ow_om) }, ignoreInit = TRUE)
  observeEvent(input$ow_ph, {
    req(!is.na(input$ow_ph))
    sel_ph(max(3, min(9, input$ow_ph)))
  }, ignoreInit = TRUE)

  observeEvent(input$reset_overwrite, {
    user_runoff_frac(NULL); user_slope_dg(NULL); user_temp_c(NULL); user_cn(NULL)
    user_clay_pct(NULL);    user_sand_pct(NULL); user_silt_pct(NULL)
    user_bd_gcm3(NULL);     user_soc_pct(NULL)
    ha <- hydroatlas_data()
    sg <- soilgrids_data()
    if (!is.null(ha)) {
      updateNumericInput(session, "ow_curve_number", value = site_cn()$CN)
      if (!is.na(ha$temp_c)) updateNumericInput(session, "ow_temp_c",   value = round(ha$temp_c, 1))
      if (!is.na(ha$slp_dg)) updateNumericInput(session, "ow_slope_dg", value = round(ha$slp_dg, 1))
    }
    if (!is.null(sg)) {
      if (!is.na(sg$clay_pct)) updateNumericInput(session, "ow_clay_pct", value = round(sg$clay_pct))
      if (!is.na(sg$sand_pct)) updateNumericInput(session, "ow_sand_pct", value = round(sg$sand_pct))
      if (!is.na(sg$silt_pct)) updateNumericInput(session, "ow_silt_pct", value = round(sg$silt_pct))
      if (!is.na(sg$bd_gcm3))  updateNumericInput(session, "ow_bd_gcm3",  value = round(sg$bd_gcm3, 2))
      if (!is.na(sg$soc_pct))  updateNumericInput(session, "ow_soc_pct",  value = round(sg$soc_pct, 2))
    } else if (!is.null(ha)) {
      if (!is.na(ha$cly_pc)) updateNumericInput(session, "ow_clay_pct", value = round(ha$cly_pc))
      if (!is.na(ha$snd_pc)) updateNumericInput(session, "ow_sand_pct", value = round(ha$snd_pc))
      if (!is.na(ha$slt_pc)) updateNumericInput(session, "ow_silt_pct", value = round(ha$slt_pc))
    }
    sel_om(if (!is.null(sg) && !is.na(sg$om_class)) sg$om_class else "moderate")
    sel_ph(if (!is.null(sg) && !is.na(sg$ph)) sg$ph else 7.0)
    updateSelectInput(session, "ow_om", selected = sel_om())
    updateNumericInput(session, "ow_ph", value = round(sel_ph(), 1))
  })

  output$soilgrids_results_panel <- renderUI({
    sg <- soilgrids_data()
    s  <- soilgrids_status_msg()

    fv <- function(x, fmt) if (!is.null(x) && !is.na(x)) sprintf(fmt, x) else "n/a"

    status_div <- if (s$type %in% c("error", "info")) {
      col <- switch(s$type, "error"="#c0392b", "success"="#027a55", "info"="#0173b2", "#555")
      bg  <- switch(s$type, "error"="#fdecea", "success"="#e8faf3", "info"="#e8f4ff", "#f5f5f5")
      tagList(
        div(style = sprintf("background:%s; border-left:3px solid %s; padding:8px 10px;
                              border-radius:3px; font-size:12px; color:%s; margin-bottom:8px;",
                            bg, col, col), s$msg),
        if (s$type == "error" && !is.null(map_coords()))
          div(style = "margin-bottom:8px;",
              actionButton("retry_soil_fetch", "Retry soil data fetch",
                           class = "btn btn-warning btn-sm"))
        else NULL
      )
    } else NULL

    if (is.null(sg)) {
      tagList(
        status_div,
        div(style = "background:#f5f5f5; border:1px dashed #ccc; border-radius:4px;
                     padding:14px; text-align:center; color:#888; font-size:13px;",
            "Click the map to fetch site data.")
      )
    } else {
      td_label <- function(txt, tip = NULL)
        tags$td(strong(txt), if (!is.null(tip)) info_tip(tip),
                style = paste0("padding:3px 6px; white-space:nowrap;",
                               " color:#2e7d32;"))
      td_val <- function(txt, active = FALSE)
        tags$td(txt,
                style = paste0("padding:3px;",
                               if (active) " font-weight:600;" else ""))

      tagList(
        status_div,
        div(style = "background:white; border:1px solid #c8e6c9; border-radius:6px; padding:14px;",
            div(style = "display:flex; justify-content:space-between; align-items:baseline;",
                h5("Soil", style = "margin-top:0; margin-bottom:4px; color:#2e7d32;"),
                actionButton("retry_soil_fetch", "Refresh",
                             class = "btn btn-default btn-xs",
                             style = "font-size:11px; padding:2px 8px;")
            ),
            tags$p(class = "src-line", sprintf("Soil properties at 15-30 cm depth sourced from %s (%s).",
                                               soil_src()$name, soil_src()$cite)),
            fluidRow(
              column(6,
                     tags$table(style = "font-size:12px; width:100%;",
                                tags$tr(td_label("Soil pH:", "Soil pH (measured in water for SoilGrids and OpenLandMap, in calcium chloride for the Soil and Landscape Grid of Australia)."),
                                        td_val(fv(sg$ph, "%.1f"))),
                                tags$tr(td_label("Organic C:", "Soil organic carbon. More organic carbon binds more pesticide to the soil."),
                                        td_val(fv(sg$soc_pct, "%.2f%%"))),
                                tags$tr(td_label("Organic matter:", "Organic matter, estimated as organic carbon x 1.724. Low below 1%, moderate 1-3%, high above 3%."),
                                        td_val(if (!is.na(sg$om_pct))
                                          sprintf("%.1f%% (%s)", sg$om_pct, sg$om_class)
                                          else "n/a"))
                     )
              ),
              column(6,
                     tags$table(style = "font-size:12px; width:100%;",
                                tags$tr(td_label("Clay:", "Share of particles finer than 0.002 mm. Clay sets the hydrologic soil group used for the curve number."),
                                        td_val(fv(sg$clay_pct, "%.0f%%"))),
                                tags$tr(td_label("Silt:", "Share of particles 0.002-0.05 mm."),
                                        td_val(fv(sg$silt_pct, "%.0f%%"))),
                                tags$tr(td_label("Sand:", "Share of particles 0.05-2 mm."),
                                        td_val(fv(sg$sand_pct, "%.0f%%"))),
                                tags$tr(td_label("Bulk density:", "Dry mass of soil per unit volume. Converts the applied pesticide to a concentration in the surface soil."),
                                        td_val(fv(sg$bd_gcm3, "%.2f g/cm\u00b3"))),
                                tags$tr(td_label("Texture class:", "USDA soil texture class from the clay, silt and sand fractions."),
                                        td_val(if (!is.null(sg$texture_class) && !is.na(sg$texture_class))
                                          sg$texture_class else "n/a"))
                     )
              )
            ),
            {
              ts <- sg$texture_sum
              tq <- sg$texture_quality %||% "unknown"
              if (!is.na(ts)) {
                tc <- if (tq == "poor") "#c0392b" else if (tq == "fair") "#d97706" else "#888"
                bg <- if (tq == "poor") "background:#fdecea;" else
                  if (tq == "fair") "background:#fef3e2;" else ""
                note <- if (tq == "poor")
                  ": low, so values may be unreliable. Try another source."
                else if (tq == "fair" && ts > 110)
                  ": the mapped fractions disagree here, so texture is uncertain. Try another source."
                else ""
                tags$p(style = sprintf(
                  "font-size:10px; margin:6px 0 2px 0; padding:3px 6px; border-radius:3px; color:%s; %s",
                  tc, bg),
                  sprintf("Clay + Silt + Sand = %.0f%%", ts),
                  info_tip("Each fraction is mapped separately, so the three need not add to exactly 100%."), note)
              } else NULL
            },
            NULL
        )
      )
    }
  })

  output$ui_pages <- renderUI({
    switch(current_page(),
           "landing"         = page_landing(),
           "location"        = page_location(),
           "pest_filter"     = page_pest_filter(),
           "targets"         = page_targets(),
           "products"        = page_products(),
           "baseline"        = page_baseline(),
           page_risk_plots()
    )
  })

  # ---- Pages ----

  # ---- Management controls ----
  mgmt_controls <- function(p, s, open = TRUE, advanced = NULL, abouts = TRUE) {
    q   <- function(txt, ...) tags$label(class = "mg-q", txt, ...)
    sub <- function(txt) tags$label(class = "mg-sub", txt)
    about <- function(title, ...) if (abouts) tags$details(class = "mg-about", tags$summary(title), tags$div(...))
    section <- function(title, tags_, ...) tags$details(class = "mg-sec", open = if (open) NA,
                                                     tags$summary(tags$span(title), tags_), tags$div(class = "mg-body", ...))
    tagList(
      section("Rate and timing", affects("A", "P", "S"),
              q("Will you apply less than the full rate? Reduction (%)"),
              sliderInput(paste0(p, "rate_reduction"), label = NULL, min = 0, max = 90, value = s$rate_reduction,
                          step = 5, post = "%"),
              q("When will you spray?"),
              dateInput(paste0(p, "app_date"), label = NULL, value = s$app_date),
              uiOutput(paste0(p, "runoff_timing_panel"))),

      section("Spraying and drift", affects("P"),
              q("How will you apply it?"),
              selectInput(paste0(p, "app_method"), label = NULL, choices = DRIFT_METHODS, selected = s$app_method),
              conditionalPanel(
                condition = sprintf("!input.%sapp_method.startsWith('Vertical')", p),
                q("What droplet size will you spray?"),
                selectInput(paste0(p, "spray_quality"), label = NULL, choices = SPRAY_QUALITIES,
                            selected = s$spray_quality),
                tags$p(class = "mg-hint", "From your nozzle chart or the product label (median droplet diameter).")),
              q("How far downwind are the nearest flowering plants or beehives? (m)"),
              selectInput(paste0(p, "habitat_dist"), label = NULL, choices = DRIFT_DIST_M, selected = s$habitat_dist),
              about("About spray drift and pollinator exposure",
                    "Contact dose to bees in flowering vegetation downwind (APVMA method), from the APVMA drift curves at the ",
                    "20 km/h label maximum wind, plus an assumed share of airborne spray and vapour.")),

      section("Irrigation", affects("A", "S"),
              q("How is the field watered?"),
              selectInput(paste0(p, "runoff_source"), label = NULL,
                          choices = c("Rain only (no irrigation)" = "rainfall",
                                      setNames(IRRIGATION_PARAMS$type, IRRIGATION_PARAMS$label)),
                          selected = s$runoff_source),
              conditionalPanel(
                condition = sprintf("input.%srunoff_source != 'rainfall'", p),
                q("When is the first irrigation after spraying?"),
                dateInput(paste0(p, "irrig_first_date"), label = NULL, value = s$irrig_first_date),
                q("How many days between irrigations?"),
                numericInput(paste0(p, "irrig_interval"), label = NULL, value = s$irrig_interval,
                             min = 1, max = 120, step = 1),
                q("How much water per irrigation? (mm; 1 ML/ha = 100 mm)"),
                numericInput(paste0(p, "irrig_depth_mm"), label = NULL, value = s$irrig_depth_mm,
                             min = 0, max = 500, step = 5),
                sub("Runoff and sediment values from"),
                selectInput(paste0(p, "irrig_source"), label = NULL,
                            choices = irrigation_source_choices(if (s$runoff_source == "rainfall") "furrow" else
                                                                  s$runoff_source),
                            selected = s$irrig_source)),
              about("About irrigation assumptions",
                    "Irrigation runoff = water applied x runoff fraction, at a fixed sediment concentration. Central values:",
                    tags$br(),
                    lapply(seq_len(nrow(IRRIGATION_PARAMS)), function(i) {
                      r <- IRRIGATION_PARAMS[i, ]
                      tagList(tags$b(paste0(r$label, ": ")),
                              sprintf("runoff %.2f (%s), sediment %.2f g/L.", r$runoff_fraction, r$RF_range, r$SC_gL),
                              tags$br())
                    }), "Irrigation affects soil fauna only through leaching.")),

      section("Runoff controls", affects("A"),
              q("Is there a grassed or vegetated strip between the field and the nearest waterway? Width (m)"),
              sliderInput(paste0(p, "Vegetative_Filter_Strips"), label = NULL,
                          min = 0, max = 30, value = s$Vegetative_Filter_Strips, step = 5),
              about("About vegetative filter strips", "Reported pesticide removal by width (Zhang et al. 2010): 5 m 62%, 10 m 83%, 20 m 92%, 30 m 93%. Any reduction set under Model assumptions is applied."),
              q("How is the field tilled?"),
              selectInput(paste0(p, "tillage_type"), label = NULL,
                          choices = c("Conventional tillage"                    = "none",
                                      "Minimum tillage (chisel, disk or ridge)" = "Minimum_Till",
                                      "No-till"                                 = "No_Till"),
                          selected = s$tillage_type),
              about("About tillage practices",
                    "Pesticide loss in runoff compared with conventional tillage (Fleming et al. 2025): minimum tillage 15% lower; no-till no change."),
              q("Are micro-dams built between the crop rows?"),
              checkboxGroupInput(paste0(p, "land_practices"), label = NULL, selected = s$land_practices,
                                 choices = list("Micro-dams" = "Micro_Dams")),
              about("About micro-dams",
                    "Micro-dams reduced pesticide loss by 87% on average in two potato studies (Klein et al. 2023). Tillage, ",
                    "micro-dams and the filter strip combine multiplicatively; the total reduction is capped at 80%.")),
      if (!is.null(advanced))
        section("Model assumptions", affects("A", "P"),
                q("Want to be more conservative? Reduce the estimated benefits of applied practices by (%)"),
                sliderInput("practice_discount_pct", label = NULL, min = 0, max = 100,
                            value = 100 - advanced$practice_benefit_pct, step = 10, post = "%"),
                tags$p(class = "mg-hint", "0% applies the full reported reduction."),
                q("Share of airborne pesticide deposited on nearby flowering vegetation (%)"),
                numericInput("vapour_dep_pct", label = NULL, value = advanced$vapour_dep_pct, min = 0, max = 100,
                             step = 0.5),
                tags$p(class = "mg-hint", "No published value; 1% is an assumption."))
    )
  }
  MGMT_CSS <- tags$style(HTML("
    .mg-sec { background:#fff; border:1px solid #c8e6c9; border-radius:8px; margin-bottom:8px; }
    .mg-sec > summary { cursor:pointer; padding:8px 12px; font-weight:700; font-size:13px; color:#2e7d32;
                        list-style:none; display:flex; justify-content:space-between; align-items:center; }
    .mg-sec > summary::-webkit-details-marker { display:none; }
    .mg-sec > summary::before { content:'\\25B8'; margin-right:6px; color:#66bb6a; }
    .mg-sec[open] > summary::before { content:'\\25BE'; }
    .mg-sec > summary > span:first-child { flex:1; }
    .mg-body { padding:2px 12px 8px 12px; }
    .mg-q { font-weight:600; font-size:12px; color:#2c3e50; display:block; margin:6px 0 2px 0; }
    .mg-sub { font-size:11px; color:#555; display:block; margin:4px 0 2px 0; }
    .mg-hint { font-size:11px; color:#666; margin:-8px 0 8px 0; }
    .mg-about { margin:-6px 0 10px 0; }
    .mg-about > summary { font-size:11px; color:#0173b2; cursor:pointer; list-style:none; user-select:none; }
    .mg-about > div { font-size:11px; color:#444; background:#fff; border-radius:4px; padding:8px; margin-top:4px;
                      border:1px solid #ddd; line-height:1.6; }
  "))
  MGMT_JS <- tags$script(HTML("
    if (!window.mgBound) { window.mgBound = true;
      document.addEventListener('toggle', function(e) {
        if (e.target.tagName === 'DETAILS' && e.target.open) { $(window).trigger('resize'); $(e.target).trigger('shown'); }
      }, true); }"))

  # ---- Navigation bar ----
  wizard_nav <- function(step, back_id = NULL, next_id = NULL,
                         next_label = "Next", back_label = "Back") {
    step_info <- list(
      list(n = 1, label = "Disclaimer"),
      list(n = 2, label = "Location"),
      list(n = 3, label = "Pests"),
      list(n = 4, label = "Targets"),
      list(n = 5, label = "Products"),
      list(n = 6, label = "Baseline"),
      list(n = 7, label = "Results")
    )
    circles <- lapply(step_info, function(s) {
      state <- if (s$n < step) "done" else if (s$n == step) "active" else "pending"
      bg  <- switch(state, done = "#029e73", active = "#0173b2", pending = "#e0e0e0")
      col <- switch(state, done = "white",   active = "white",   pending = "#999")
      lbl_col <- if (state == "pending") "#aaa" else "#444"
      connector <- if (s$n < length(step_info))
        div(style = "width:28px; height:2px; background:#ddd; flex-shrink:0; margin-bottom:18px;")
      else NULL
      tagList(
        actionLink(paste0("wiz_goto_", s$n), title = paste0("Go to step ", s$n, ": ", s$label),
                   style = "text-decoration:none; cursor:pointer;", label =
        div(style = "display:flex; flex-direction:column; align-items:center; flex-shrink:0;",
            div(style = sprintf("width:34px; height:34px; border-radius:50%%; background:%s;
                                 color:%s; display:flex; align-items:center; justify-content:center;
                                 font-weight:700; font-size:15px; flex-shrink:0;", bg, col),
                as.character(s$n)),
            div(style = sprintf("font-size:10px; color:%s; margin-top:3px; white-space:nowrap;
                                 font-weight:%s;", lbl_col,
                                if (state == "active") "700" else "400"),
                s$label)
        )),
        connector
      )
    })
    div(style = "background:white; border-bottom:2px solid #e8e8e8; padding:14px 24px;
                  margin:-15px -15px 32px -15px; position:sticky; top:0; z-index:200;
                  box-shadow:0 2px 8px rgba(0,0,0,0.06);",
        div(class = "wiz-outer",
            style = "display:flex; align-items:center; justify-content:space-between;
                      max-width:1100px; margin:0 auto;",
            if (!is.null(back_id))
              actionButton(back_id, back_label,
                           class = "wiz-btn",
                           style = "font-size:16px; font-weight:600; padding:10px 28px;
                                     background:#f5f5f5; border:2px solid #ccc; border-radius:8px;
                                     min-width:120px; cursor:pointer;")
            else div(class = "wiz-btn", style = "min-width:120px;"),
            div(class = "wiz-steps", style = "display:flex; align-items:flex-start;", circles),
            if (!is.null(next_id))
              actionButton(next_id, next_label,
                           class = "btn btn-primary wiz-btn",
                           style = "font-size:16px; font-weight:600; padding:10px 28px;
                                     border-radius:8px; min-width:120px;")
            else div(class = "wiz-btn", style = "min-width:120px;")
        )
    )
  }

  # ---- Step 1: disclaimer ----
  page_landing <- function() {
    fluidPage(
      wizard_nav(1, back_id = NULL, next_id = "start_assessment", next_label = "I understand, proceed"),
      div(style = "max-width: 780px; margin: 0 auto;",
          h2("Pesticide Risk Assessment Tool",
             style = "color:#0173b2; margin-top:0; margin-bottom:6px;"),
          p(style = "color:#555; margin-bottom:28px;",
            "Compares the relative risk of pesticides to waterways, pollinators, soil fauna and human health."),
          div(class = "disclaimer-box",
              h4("Disclaimer", style = "margin-top:0; color:#92400e;"),
              tags$ul(
                tags$li("For information only; not a substitute for professional agronomic or environmental advice."),
                tags$li("Results are relative comparisons under simplified assumptions, not safety determinations."),
                tags$li("Tank mixes, metabolites and indirect ecological effects are not modelled."),
                tags$li("Registered uses and rates are read automatically from Australian product labels (APVMA). They apply in Australia only, may contain errors and may have changed since."),
                tags$li(strong("Whatever this tool shows, always follow the current product label and your local regulations.")),
                tags$li(strong("No warranty is given. We accept no responsibility for errors in the label information or for decisions made with this tool.")))
          ),
          div(class = "load-box",
              tags$b("Have a saved assessment?"),
              " Load the file from Export results and settings to return to its results (accepts the disclaimer).",
              fileInput("load_assessment", label = NULL, accept = c(".json", "application/json"),
                        buttonLabel = "Load file", placeholder = "pesticide_assessment_....json", width = "100%"))
      )
    )
  }

  observeEvent(input$start_assessment, { current_page("location") })
  observeEvent(input$back_to_landing,  { current_page("landing") })

  # ---- Step 2: location ----
  page_location <- function() {
    fluidPage(
      wizard_nav(2, back_id = "back_to_landing", next_id = "proceed_to_pest"),
      div(style = "max-width: 1000px; margin: 0 auto;",
          h2("Where is your field?", style = "margin-top:0; margin-bottom:6px; color:#01745f;"),
          p(style = "color:#555; margin-bottom:12px;",
            if (HAS_DRAW) "Click your field on the map, or outline it with the draw tool (on the map, below the zoom buttons) to average values over the field. "
            else "Click your field on the map. ",
            "Soil, water balance and daily rainfall load automatically; ", tags$b("then click Next.")),
          div(style = "background:#f0f8ff; border:1px solid #cce0f5; border-radius:6px;
                        padding:10px 14px; margin-bottom:14px; display:flex;
                        align-items:center; gap:14px; flex-wrap:wrap;",
              div(style = "flex:0 0 auto;",
                  tags$label(style = "font-weight:600; font-size:12px; color:#0173b2; white-space:nowrap;",
                             "Soil data source:")),
              div(style = "flex:1 1 260px; min-width:220px;",
                  selectInput("soil_source_pref", label = NULL,
                              choices = c(
                                "Auto (recommended)"                       = "auto",
                                "SoilGrids (global, 250 m)"                       = "soilgrids",
                                "Soil and Landscape Grid of Australia (90 m)"     = "slga",
                                "OpenLandMap (global, 250 m)"                     = "olm"
                              ),
                              selected = "auto",
                              width = "100%")),
          ),
          leafletOutput("location_map", height = "420px"),
          tags$div(id = "location_fetch_overlay",
                   tags$span(class = "fetch-spinner-ring"),
                   "Fetching site data (up to a minute)...",
                   tags$span(id = "fetch_overlay_src",
                             style = "display:block; font-size:11px; color:#5a9fd4; margin-top:4px; font-weight:400;")),
          uiOutput("fetch_status_panel"),
          uiOutput("hydroatlas_results_panel"),
          uiOutput("rainfall_results_panel"),
          uiOutput("soilgrids_results_panel"),
          uiOutput("overwrite_panel"),
          uiOutput("location_refs")
      )
    )
  }

  observeEvent(input$proceed_to_pest, { current_page("pest_filter") })

  output$location_refs <- renderUI({
    ha <- hydroatlas_data(); if (is.null(ha)) return(NULL)
    rain <- rain_daily_data()
    ref_list(c("linke", if (!is.null(rain)) c(rain_source_info(rain)$ref, "nrcs2004", "nrcs2009"),
               if (!is.null(soilgrids_data())) soil_src()$ref))
  })

  output$rainfall_results_panel <- renderUI({
    rain <- rain_daily_data(); ha <- hydroatlas_data()
    if (is.null(ha)) return(NULL)
    box <- "margin-top: 12px; background: white; border: 1px solid #b2dfdb; border-radius: 6px; padding: 14px;"
    if (is.null(rain))
      return(div(style = box, h5("Daily rainfall", style = "margin-top:0; color:#01745f;"),
                 p(style = "font-size:12px; color:#92400e;",
                   "Not loaded: ", DEFAULT_DAYS_TO_RUNOFF, " days to first runoff and the HydroATLAS runoff fraction are used.")))
    yrs <- format(range(rain$date), "%Y")
    yr  <- format(rain$date, "%Y"); full <- names(which(table(yr) >= 365))
    map <- mean(tapply(pmax(0, rain$P_mm), yr, sum, na.rm = TRUE)[full])
    cn  <- site_cn(); q <- annual_cn_runoff(rain, cn$CN)
    src <- rain_source_info(rain)
    tr <- function(lab, val, tip) tags$tr(tags$td(strong(lab), info_tip(tip), style = "padding:3px;"),
                                          tags$td(val, style = "padding:3px;"))
    div(style = box,
        h5("Daily rainfall", style = "margin-top:0; margin-bottom:2px; color:#01745f;"),
        tags$p(class = "src-line", sprintf("Long-term (%s-%s) daily rainfall sourced from %s (%s).", yrs[1], yrs[2], src$name, src$cite)),
        tags$table(style = "width:100%; font-size:12px;",
          tr("Mean annual rainfall:", sprintf("%.0f mm/yr (HydroATLAS %.0f)", map, ha$pre_mm %||% NA_real_),
             "Mean of the complete years in the daily record. The HydroATLAS sub-basin value is shown in brackets."),
          tr("Curve number:", sprintf("%g (row crop, soil group %s%s)", cn$CN, cn$hsg,
                                      if (!identical(cn$CN, cn$CN_table)) ", your value" else ""),
             "Runoff curve number for row crops in good condition on this hydrologic soil group, which is set from clay and sand content (USDA-NRCS 2004, 2009). Higher values give more runoff. It can be changed under Edit retrieved values."),
          tr("Mean annual runoff:", sprintf("%.0f mm/yr (HydroATLAS %.0f)", q, ha$run_mm %||% NA_real_),
             "Runoff calculated each day from the daily rainfall with this curve number, as a mean per year. The HydroATLAS sub-basin value is shown in brackets.")),
        if (is.finite(q) && isTRUE(ha$run_mm > 0) && q < 0.25 * ha$run_mm)
          tags$p(style = "font-size:12px; color:#92400e; margin:4px 0 0 0;",
                 sprintf("Daily curve-number runoff is well below the HydroATLAS value (%.0f vs %.0f mm/yr), so rain runoff events may be under-predicted here. Consider raising the curve number under Edit retrieved values.",
                         q, ha$run_mm)))
  })

  observe({
    req(current_page() == "location")
    if (!is.null(hydroatlas_data())) shinyjs::enable("proceed_to_pest")
    else                              shinyjs::disable("proceed_to_pest")
  })

  # ---- Step 3: pest groups ----
  page_pest_filter <- function() {
    fluidPage(
      wizard_nav(3, back_id = "back_from_pest", next_id = "proceed_to_risks"),
      div(style = "max-width: 860px; margin: 0 auto;",
          h2("What do you need to control?",
             style = "margin-top:0; color:#01745f;"),
          p(style = "color:#555; margin-bottom:20px;",
            "Select one or more groups, or none to include all chemicals."),
          div(style = "display:flex; justify-content:flex-end; margin-bottom:10px;",
              actionButton("pest_group_clear", "Clear selection",
                           class = "btn btn-default btn-sm")),
          uiOutput("pest_tile_grid")
      )
    )
  }

  # ---- Step 4: crop and target pests ----
  TGT_CSS <- tags$style(HTML("
    .tgt-card { background:#fff; border:1px solid #d6e4f0; border-radius:12px; padding:16px 18px; height:100%;
                box-shadow:0 1px 4px rgba(0,0,0,0.05); }
    .tgt-card-title { font-size:15px; font-weight:700; color:#01745f; margin-bottom:2px; }
    .tgt-hint { font-size:12px; color:#6b7a8d; margin:0 0 10px 0; }
    .tgt-banner { margin-top:16px; border-radius:12px; padding:14px 18px; display:flex; gap:18px; align-items:center; }
    .tgt-count { font-size:30px; font-weight:800; line-height:1; }
    .tgt-chip { display:inline-block; font-size:11px; padding:2px 9px; margin:2px 3px 2px 0; border-radius:10px;
                background:#fff; border:1px solid #a5d6a7; color:#2e7d32; }
    .tp-table { width:100%; border-collapse:separate; border-spacing:0; background:#fff; border:1px solid #d6e4f0;
                border-radius:12px; overflow:hidden; }
    .tp-table th { background:#f0f6fb; color:#0173b2; font-size:12px; font-weight:700; padding:10px 12px;
                   border-bottom:1px solid #d6e4f0; }
    .tp-table td { padding:8px 12px; vertical-align:middle; border-bottom:1px solid #eef3f8; font-size:13px; }
    .tp-table tr:last-child td { border-bottom:none; }
    .tp-table tr:hover td { background:#fafcfe; }
    .tp-table .form-group { margin-bottom:0; }
    .tp-sub { font-size:11px; color:#8a97a6; }
    .tp-def { color:#555; white-space:nowrap; }
    .tp-conv { font-size:11px; color:#2e7d32; margin-top:3px; }
    .tp-unit { min-width:64px; }
    .tp-chem { color:#01745f; text-decoration:none; border-bottom:1px dotted #01745f; }
    .tp-chem:hover { color:#0173b2; border-bottom-color:#0173b2; }
    .tp-label-link { font-size:11px; color:#0173b2; display:inline-block; margin-top:3px; }
    .lab-table { width:100%; font-size:12px; border-collapse:collapse; }
    .lab-table th { background:#f0f6fb; color:#0173b2; padding:6px 8px; text-align:left; }
    .lab-table td { padding:5px 8px; border-bottom:1px solid #eef3f8; vertical-align:top; }
    .lab-table tr.lab-sel td { background:#fff8e1; font-weight:600; }
    .lab-table tr.lab-sub td { border-bottom:1px solid #e3e8ee; padding:0 8px 6px 18px; background:#fbfcfd; }
    .lab-lines { margin:2px 0 0 0; padding-left:14px; font-size:11px; color:#444; font-weight:400; }
    .lab-note { color:#888; }
    .opt-callout { background:#fff8e1; border:1px solid #ffe082; border-left:4px solid #f9a825; border-radius:8px;
                   padding:10px 14px; margin-bottom:16px; font-size:14px; color:#5d4037; }
    .tp-legend { font-size:12px; color:#555; margin:-4px 0 12px 0; display:flex; flex-direction:column; gap:2px; }
  "))
  page_targets <- function() {
    fluidPage(
      TGT_CSS,
      wizard_nav(4, back_id = "back_from_targets", next_id = "proceed_from_targets"),
      div(style = "max-width: 980px; margin: 0 auto;",
          h2("What are you targeting, and where?", style = "margin-top:0; color:#01745f;"),
          div(class = "opt-callout", tags$b("Optional."),
              " Choose your crop and target pests to narrow the chemicals and use their label rates, or leave blank and ",
              tags$b("click Next to skip.")),
          fluidRow(
            column(6, div(class = "tgt-card",
                          div(class = "tgt-card-title", "Crop or situation"),
                          p(class = "tgt-hint", "Uses registered in Australia."),
                          uiOutput("tgt_host_ui"))),
            column(6, div(class = "tgt-card",
                          div(class = "tgt-card-title", "Target pests"),
                          p(class = "tgt-hint", "Pests in your Step 3 groups. ",
                            actionLink("tgt_toggle_all", uiOutput("tgt_toggle_label", inline = TRUE))),
                          uiOutput("tgt_pest_ui")))
          ),
          uiOutput("tgt_match_summary"),
          tags$p(class = "src-line", style = "margin-top:14px;",
                 "Registered crops, situations and pests sourced from PubCRIS (APVMA 2026a)."),
          ref_list("apvma_pubcris")
      )
    )
  }

  # ---- Step 5: products and rates ----
  page_products <- function() {
    fluidPage(
      TGT_CSS,
      wizard_nav(5, back_id = "back_from_products", next_id = "proceed_from_products"),
      div(style = "max-width: 1060px; margin: 0 auto;",
          h2("Which products, at what rates?", style = "margin-top:0; color:#01745f;"),
          div(class = "opt-callout", tags$b("Optional."),
              " Choose a product or enter your own rate, or leave blank and ", tags$b("click Next"),
              " to use the default rate shown."),
          div(class = "tp-legend",
              tags$span(tags$b("kg ai/ha"), " = kilograms of active ingredient per hectare (the same whichever product carries it)."),
              tags$span(tags$b("L/ha or kg/ha"), " with a product = litres or kilograms of that product per hectare, as on its label.")),
          div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:10px;",
              uiOutput("tgt_products_count"),
              actionButton("reset_rates", "Reset all to default rates", class = "btn btn-default btn-sm")),
          uiOutput("tgt_products_ui"),
          tags$p(class = "src-line", style = "margin-top:12px;",
                 "Label rates sourced from Australian product labels (APVMA 2026b) and registrations from PubCRIS (APVMA 2026a). ",
                 "Click a chemical for its details, or Label information for each product's label rates."),
          ref_list(c("apvma_pubcris", "apvma_labels")),
          tags$script(HTML("
            if (!window.tpBound) {
              window.tpBound = true;
              window.tpUpdate = function(row) {
                var m = row.data('units'); if (!m) return;
                var sel = row.find('select').val() || 'default';
                var e = m[sel] || m['default'];
                row.find('.tp-unit').text(e.u);
                var inp = row.find('input[type=number]'); inp.attr('placeholder', e.d);
                var v = parseFloat(inp.val());
                row.find('.tp-conv').text(isNaN(v) ? 'Default rate used' : '= ' + Number((v * e.f).toPrecision(3)) + ' kg ai/ha');
              };
              $(document).on('change', '.tp-row select', function() { tpUpdate($(this).closest('.tp-row')); });
              $(document).on('input change', '.tp-row input[type=number]', function() { tpUpdate($(this).closest('.tp-row')); });
            }"))
      )
    )
  }

  # ---- Step 6: baseline practice ----
  page_baseline <- function() {
    s <- isolate(baseline_settings())
    fluidPage(
      wizard_nav(6, back_id = "back_from_baseline", next_id = "proceed_from_baseline", next_label = "View risks"),
      div(style = "max-width: 760px; margin: 0 auto;",
          h2("How will you spray, and how is the field managed?", style = "margin-top:0; color:#01745f;"),
          p(style = "color:#555; margin-bottom:16px;",
            "Answer for your planned application and current practice. This is your starting point; on the results page you can try changes."),
          div(style = paste0("background: linear-gradient(160deg, #f1f8e9 0%, #e8f4fd 100%);",
                             "border: 2px solid #66bb6a; border-radius: 10px; padding: 14px 18px;"),
              MGMT_CSS, MGMT_JS, affects_legend(),
              mgmt_controls("b_", s, open = TRUE, abouts = FALSE)),
          BASELINE_SOURCES
      )
    )
  }

  # ---- Step 7: results ----
  page_risk_plots <- function() {
    fluidPage(
      wizard_nav(7, back_id = "back_from_risks", next_id = NULL),
      div(style = "display:flex; justify-content:flex-end; margin:-16px 0 20px 0;",
          actionButton("restart_from_p3", "New assessment, same field",
                       class = "btn btn-warning btn-sm")),
      fluidRow(
        column(3,
               style = "overflow-y: auto; max-height: 90vh; padding-right: 10px;",

               div(style = paste0(
                 "background: linear-gradient(160deg, #f1f8e9 0%, #e8f4fd 100%);",
                 "border: 2px solid #66bb6a; border-radius: 10px;",
                 "padding: 14px 14px 10px 14px; margin-bottom: 14px;"
               ),

               div(style = "margin-bottom:6px;",
                   tags$span(style = "font-size:14px; font-weight:700; color:#2e7d32; line-height:1.3;",
                             "Could management practices reduce your risk?")
               ),

               MGMT_CSS, MGMT_JS, affects_legend(),
               mgmt_controls("", isolate(scenario_init() %||% baseline_settings()), open = FALSE,
                             advanced = isolate(adv_init())),
               ),

               hr(),
               actionButton("clear_settings", "Clear all and restart",
                            class = "btn-warning btn-sm", style = "width:100%; margin-top:4px;")
        ),
        column(9,
               tabsetPanel(
                 tabPanel("Risk summary",
                          div(style = paste0("background:#e8f4fd; border:1px solid #90caf9; border-radius:6px;",
                                             "padding:8px 14px; margin-bottom:12px; font-size:13px; color:#0173b2;"),
                              "Click a row for chemical details, products and health information."
                          ),
                          DTOutput("risk_classification_table")
                 ),
                 tabPanel("Aquatic risk",
                          plotOutput("pesticidePlot", height = "700px"),
                          tags$p(style = "font-size:12px; color:#555; margin-top:4px;", ARROW_NOTE),
                          uiOutput("missing_note_aquatic"),
                          br(),
                          downloadButton("download_aquatic_plot", "Download PNG (600 DPI)")
                 ),
                 tabPanel("Pollinator risk",
                          plotOutput("terrestrial_risk_plot", height = "700px"),
                          tags$p(style = "font-size:12px; color:#555; margin-top:4px;", ARROW_NOTE),
                          uiOutput("missing_note_pollinator"),
                          br(),
                          downloadButton("download_above_plot", "Download PNG (600 DPI)")
                 ),
                 tabPanel("Soil fauna risk",
                          plotOutput("below_risk_plot", height = "700px"),
                          tags$p(style = "font-size:12px; color:#555; margin-top:4px;", ARROW_NOTE),
                          uiOutput("missing_note_soil"),
                          br(),
                          downloadButton("download_below_plot", "Download PNG (600 DPI)")
                 ),
                 tabPanel("Human risk",
                          plotOutput("HumanMatrix", height = "800px"),
                          br(),
                          downloadButton("download_human_plot", "Download PNG (600 DPI)")
                 )
               )
        )
      ),
      br(),
      fluidRow(
        column(12,
               style = "text-align: center;",
               downloadButton("download_pdf_report", "Download PDF report",
                              class = "btn-primary", style = "margin-right:12px;"),
               downloadButton("export_settings", "Export results and settings", class = "btn-info")
        )
      ),
      fluidRow(column(12, RESULTS_SOURCES))
    )
  }

  # ---- Navigation observers ----
  observeEvent(input$back_from_pest,  { current_page("location") })
  observeEvent(input$back_from_risks, { current_page("baseline") })
  observeEvent(input$back_from_baseline, {
    baseline_settings(read_settings("b_"))
    current_page("products")
  })
  observeEvent(input$back_from_targets,     { current_page("pest_filter") })
  observeEvent(input$proceed_from_targets,  { if (save_targets()) current_page("products") })
  observeEvent(input$back_from_products,    { save_products(); current_page("targets") })
  observeEvent(input$proceed_from_products, { save_products(); current_page("baseline") })

  WIZ_PAGES <- c("landing", "location", "pest_filter", "targets", "products", "baseline", "risk_plots")
  goto_step <- function(n) {
    cur <- current_page(); target <- WIZ_PAGES[n]
    if (identical(cur, target)) return(invisible(NULL))
    if (n >= 3 && is.null(hydroatlas_data())) {
      showNotification("Select a field location in Step 2 first.", type = "warning", duration = 4)
      return(invisible(NULL))
    }
    if (identical(cur, "baseline")) baseline_settings(read_settings("b_"))
    if (identical(cur, "targets") && !save_targets()) return(invisible(NULL))
    if (identical(cur, "products")) save_products()
    current_page(target)
  }
  lapply(seq_along(WIZ_PAGES), function(n)
    observeEvent(input[[paste0("wiz_goto_", n)]], goto_step(n), ignoreInit = TRUE))
  observeEvent(input$proceed_from_baseline, {
    baseline_settings(read_settings("b_"))
    scenario_init(NULL)
    current_page("risk_plots")
  })

  # ---- Reload a saved assessment ----
  as_settings <- function(x) {
    s <- default_settings()
    if (is.list(x)) for (n in intersect(names(x), names(s))) if (!is.null(x[[n]]) && length(x[[n]])) s[[n]] <- x[[n]]
    s$app_date <- as.Date(s$app_date); s$irrig_first_date <- as.Date(s$irrig_first_date)
    s$habitat_dist <- as.numeric(s$habitat_dist)
    s$land_practices <- intersect(as.character(unlist(x$land_practices %||% character(0))), PRACTICES_KEPT)
    if (s$tillage_type %in% names(TILLAGE_ALIASES)) s$tillage_type <- TILLAGE_ALIASES[[s$tillage_type]]
    s
  }
  na_fill <- function(x) if (is.list(x)) lapply(x, function(v) if (is.null(v) || (is.list(v) && !length(v))) NA_real_ else v) else x
  num1 <- function(x) if (is.null(x) || !length(x) || all(is.na(x))) NULL else as.numeric(x)[1]
  restore_assessment <- function(js) {
    if (is.null(js$map_lat) || is.null(js$map_lon)) stop("no field location in the file")
    ring <- js$field_outline$ring
    ring <- if (is.null(ring)) NULL else if (is.matrix(ring)) lapply(seq_len(nrow(ring)), function(i) as.numeric(ring[i, ]))
            else lapply(ring, as.numeric)
    st <- js$site
    if (!is.null(st$hydroatlas)) {
      map_coords(list(lat = js$map_lat, lon = js$map_lon))
      field_geom(if (!is.null(ring)) list(ring = ring, area_ha = ring_area_ha(ring)) else NULL)
      hydroatlas_data(na_fill(st$hydroatlas))
      soilgrids_data(if (!is.null(st$soil)) na_fill(st$soil) else NULL)
      soil_source_name(js$soil_source %||% NULL); soil_source_id(NULL)
      rusle_data(if (!is.null(st$rusle)) na_fill(st$rusle) else NULL)
      rain_daily_data(if (!is.null(st$rain)) structure(data.frame(date = as.Date(st$rain$date), P_mm = as.numeric(st$rain$P_mm)),
                                                       source = st$rain$source) else NULL)
      if (!is.null(st$soil_ph)) sel_ph(as.numeric(st$soil_ph))
      if (!is.null(st$soil_om)) sel_om(st$soil_om)
      hydroatlas_status_msg(list(type = "success", msg = "Site data loaded from the saved assessment."))
      soilgrids_status_msg(list(type = "none", msg = ""))
    } else {
      fetch_site(js$map_lat, js$map_lon, ring)
    }
    uo <- js$user_overrides
    user_runoff_frac(num1(uo$runoff_fraction)); user_slope_dg(num1(uo$slope_deg)); user_temp_c(num1(uo$temp_c))
    user_clay_pct(num1(uo$clay_pct)); user_sand_pct(num1(uo$sand_pct)); user_silt_pct(num1(uo$silt_pct))
    user_bd_gcm3(num1(uo$bd_gcm3)); user_soc_pct(num1(uo$soc_pct)); user_cn(num1(uo$curve_number))
    selected_pest_groups(upgrade_pest_groups(as.character(unlist(js$pest_groups %||% character(0)))))
    hs <- unique(name_key(crop_base(unlist(js$target_crops %||% character(0)))))
    ps <- unique(name_key(clean_name(unlist(js$target_pests %||% character(0)))))
    tgt_epoch(isolate(tgt_epoch()) + 1); rate_epoch(isolate(rate_epoch()) + 1)
    target_hosts(hs); target_pests(ps)
    m <- match_for(hs, ps); target_actives(if (m$refined) m$actives else NULL)
    pl <- js$product_plan; merged <- gone <- character(0)
    if (is.data.frame(pl) && nrow(pl)) {
      old <- pl$active %in% names(ACTIVE_RENAMES) & !pl$active %in% ALL_ACTIVES
      pl$active[old] <- unname(ACTIVE_RENAMES[pl$active[old]])
      mg <- pl$active %in% names(ACTIVE_MERGES) & !pl$active %in% ALL_ACTIVES
      if (any(mg)) {
        merged <- unique(pl$active[mg])
        pl$active[mg] <- unname(ACTIVE_MERGES[pl$active[mg]])
        pl <- pl[order(pl$active, -pl$kg_ha), , drop = FALSE]
        pl <- pl[!duplicated(pl$active), , drop = FALSE]
      }
      gone <- unique(pl$active[!pl$active %in% ALL_ACTIVES])
      pl <- pl[pl$active %in% ALL_ACTIVES, , drop = FALSE]
    }
    product_plan(if (is.data.frame(pl) && nrow(pl)) pl else NULL)
    baseline_settings(as_settings(js$baseline_settings))
    scenario_init(as_settings(js$scenario_settings %||% js$baseline_settings))
    adv_init(list(practice_benefit_pct = num1(js$practice_benefit_pct) %||% PRACTICE_BENEFIT_DEFAULT,
                  vapour_dep_pct = num1(js$airborne_deposition_pct) %||% VAPOUR_DEP_DEFAULT_PCT))
    old_p <- unique(as.character(unlist(c(js$baseline_settings$land_practices, js$scenario_settings$land_practices))))
    invisible(list(merged = merged, gone = gone, practices = setdiff(old_p, PRACTICES_KEPT)))
  }
  observeEvent(input$load_assessment, {
    f <- input$load_assessment; req(f)
    res <- tryCatch({
      js <- jsonlite::fromJSON(f$datapath, simplifyVector = TRUE)
      restore_assessment(js)
    }, error = function(e) conditionMessage(e))
    if (is.list(res)) {
      current_page("risk_plots")
      showNotification("Saved assessment loaded.", type = "message", duration = 4)
      for (a in res$merged)
        showNotification(sprintf("%s is now listed as %s. Check its rate at Step 5.", a, ACTIVE_MERGES[[a]]),
                         type = "warning", duration = NULL)
      if (length(res$gone))
        showNotification(paste0("No longer in the chemical data, so the saved rate was dropped: ",
                                paste(res$gone, collapse = ", "), "."), type = "warning", duration = NULL)
      if (length(res$practices))
        showNotification(paste0("No longer in the tool, so not applied: ",
                                paste(tolower(gsub("_", " ", res$practices)), collapse = ", "), "."),
                         type = "warning", duration = NULL)
    } else showNotification(paste("Could not load this file:", res), type = "error", duration = 8)
  })

  observeEvent(input$proceed_to_risks, { current_page("targets") })

  output$pest_tile_grid <- renderUI({
    selected <- selected_pest_groups()
    tile_defs <- list(
      list(id = "grass_weeds",     title = "Grass Weeds",
           desc = "Ryegrass, barnyard grass, wild oats, summer grasses"),
      list(id = "broadleaf_weeds", title = "Broadleaf Weeds",
           desc = "Thistles, pigweed, bindweed, capeweed, fleabane"),
      list(id = "sedges_woody",    title = "Sedges, Rushes & Woody Weeds",
           desc = "Nutgrass, sedges, vines, regrowth, woody weeds"),
      list(id = "sucking_insects", title = "Sucking Insects & Mites",
           desc = "Aphids, mites, whiteflies, thrips, mirids, scale"),
      list(id = "caterpillars",    title = "Caterpillars & Borers",
           desc = "Helicoverpa, armyworms, loopers, moths, borers"),
      list(id = "beetles_grubs",   title = "Beetles, Grubs & Other Insects",
           desc = "Beetles, weevils, canegrubs, grasshoppers, fruit fly, ants"),
      list(id = "soil_pests",      title = "Soil Pests & Nematodes",
           desc = "Wireworms, nematodes, termites, earth mites, snails"),
      list(id = "leaf_diseases",   title = "Leaf & Fruit Diseases",
           desc = "Rust, mildews, botrytis, blights, leaf spots"),
      list(id = "root_diseases",   title = "Root & Soil-borne Diseases",
           desc = "Phytophthora, pythium, fusarium, root rots")
    )
    tiles <- lapply(tile_defs, function(t) {
      cls <- if (t$id %in% selected) "pest-tile pest-tile-selected" else "pest-tile"
      actionButton(
        paste0("pest_tile_", t$id), label = NULL, class = cls,
        tagList(
          tags$span(style = "font-weight:700; font-size:14px; display:block; margin-bottom:4px;", t$title),
          tags$span(style = "font-size:12px; color:#6b7a8d; line-height:1.4;", t$desc)
        )
      )
    })
    div(class = "pest-tile-grid", tiles)
  })

  observeEvent(input$restart_from_p3, {
    selected_pest_groups(character(0))
    clear_targets()
    baseline_settings(default_settings())
    scenario_init(NULL)
    adv_init(list(practice_benefit_pct = PRACTICE_BENEFIT_DEFAULT, vapour_dep_pct = VAPOUR_DEP_DEFAULT_PCT))
    updateSliderInput(session, "practice_discount_pct", value = 100 - PRACTICE_BENEFIT_DEFAULT)
    updateNumericInput(session, "vapour_dep_pct", value = VAPOUR_DEP_DEFAULT_PCT)
    current_page("pest_filter")
  })

  # ---- Pest groups (Step 3) ----
  lapply(names(pest_group_regex), function(gid) {
    observeEvent(input[[paste0("pest_tile_", gid)]], {
      current <- selected_pest_groups()
      if (gid %in% current) selected_pest_groups(setdiff(current, gid))
      else selected_pest_groups(c(current, gid))
    }, ignoreInit = TRUE)
  })

  observeEvent(input$pest_group_clear, { selected_pest_groups(character(0)) })

  pesticide_selection <- reactive({
    vet <- if ("typedesc" %in% names(raw_data)) grepl("^VM", dplyr::coalesce(raw_data$typedesc, "")) else FALSE
    d <- raw_data[!vet, "PesticideactiveIngredient", drop = FALSE]
    groups <- selected_pest_groups()
    if (length(groups)) d <- d[d$PesticideactiveIngredient %in% unique(unlist(GROUP_ACTIVES[groups])), , drop = FALSE]
    d
  })

  # ---- Crop, pests, products and rates (Steps 4-5) ----
  target_actives <- reactiveVal(NULL)
  target_hosts   <- reactiveVal(character(0))
  target_pests   <- reactiveVal(character(0))
  product_plan   <- reactiveVal(NULL)
  tgt_show_all   <- reactiveVal(FALSE)
  tgt_epoch  <- reactiveVal(0)
  rate_epoch <- reactiveVal(0)
  tgt_id <- function(kind) paste0("tgt_", kind, "_", tgt_epoch())

  step3_actives  <- reactive(unique(pesticide_selection()$PesticideactiveIngredient))
  step3_products <- reactive(unique(na.omit(raw_data$pcode_int[raw_data$PesticideactiveIngredient %in% step3_actives()])))
  uses_step3 <- reactive({
    if (is.null(uses_raw)) return(NULL)
    uses_raw[uses_raw$pcode %in% step3_products(), , drop = FALSE]
  })
  title_case <- function(x) {
    x <- sub("-?\\s*Ref Label Life Stage\\s*$", "", x, ignore.case = TRUE)
    t <- tools::toTitleCase(tolower(x))
    has <- grepl(" - [^()]+$", t)
    t[has] <- vapply(which(has), function(i) {
      parts <- strsplit(t[i], " - ", fixed = TRUE)[[1]]; s <- parts[length(parts)]
      w <- strsplit(s, " ", fixed = TRUE)[[1]]; if (length(w) > 1) w[-1] <- tolower(w[-1])
      paste(c(parts[-length(parts)], paste(w, collapse = " ")), collapse = " - ")
    }, "")
    t
  }
  n_products <- function(n) sprintf("%d %s", as.integer(n), ifelse(as.integer(n) == 1, "product", "products"))

  output$tgt_host_ui <- renderUI({
    u <- uses_step3()
    if (is.null(u)) return(tags$p(style = "font-size:12px; color:#a00;",
      "Registration data (Pesticides_APVMA_uses.rds) not found, so crops and pests cannot be listed."))
    ph  <- unique(u[!is.na(u$host_grp), c("pcode", "host_grp")])
    cnt <- sort(table(droplevels(ph$host_grp)), decreasing = TRUE)
    selectizeInput(tgt_id("hosts"), label = NULL, multiple = TRUE, width = "100%",
                   choices = setNames(names(cnt), sprintf("%s (%s)", title_case(host_label(names(cnt))), n_products(cnt))),
                   selected = intersect(isolate(input[[tgt_id("hosts")]]) %||% isolate(target_hosts()), names(cnt)),
                   options = list(placeholder = "Any crop or situation (type to search)", plugins = list("remove_button")))
  })

  observeEvent(input$tgt_toggle_all, tgt_show_all(!tgt_show_all()))
  output$tgt_toggle_label <- renderUI(if (tgt_show_all()) "Show only my Step 3 groups" else "Show all registered pests")

  output$tgt_pest_ui <- renderUI({
    u <- uses_step3(); req(!is.null(u))
    hs <- input[[tgt_id("hosts")]] %||% character(0)
    if (length(hs)) u <- u[u$host_grp %in% hs, , drop = FALSE]
    pp  <- unique(u[!is.na(u$pest_grp), c("pcode", "pest_grp")])
    cnt <- sort(table(droplevels(pp$pest_grp)), decreasing = TRUE)
    nm  <- names(cnt)
    rx  <- active_pest_regex()
    if (!is.null(rx) && !tgt_show_all()) {
      keep <- grepl(rx, pest_label(nm), ignore.case = TRUE, perl = TRUE)
      if (any(keep)) { nm <- nm[keep]; cnt <- cnt[keep] }
    }
    selectizeInput(tgt_id("pests"), label = NULL, multiple = TRUE, width = "100%",
                   choices = setNames(nm, sprintf("%s (%s)", title_case(pest_label(nm)), n_products(cnt))),
                   selected = intersect(isolate(input[[tgt_id("pests")]]) %||% isolate(target_pests()), nm),
                   options = list(placeholder = "Any registered pest (type to search)", plugins = list("remove_button")))
  })

  match_for <- function(hs, ps) {
    acts <- step3_actives(); u <- uses_step3()
    hs <- hs %||% character(0); ps <- ps %||% character(0)
    if (is.null(u) || (length(hs) == 0 && length(ps) == 0))
      return(list(actives = sort(acts), pcodes = step3_products(), refined = FALSE))
    if (length(hs)) u <- u[u$host_grp %in% hs, , drop = FALSE]
    if (length(ps)) u <- u[u$pest_grp %in% ps, , drop = FALSE]
    pc <- unique(u$pcode)
    a  <- unique(raw_data$PesticideactiveIngredient[raw_data$pcode_int %in% pc])
    list(actives = sort(intersect(acts, a)), pcodes = pc, refined = TRUE)
  }
  target_match <- reactive(match_for(input[[tgt_id("hosts")]], input[[tgt_id("pests")]]))
  saved_match  <- reactive(match_for(target_hosts(), target_pests()))

  # ---- Default rates ----
  use_rate_rows <- function(hs, ps, pcodes, acts) {
    if (is.null(use_rates) || (!length(hs) && !length(ps))) return(NULL)
    take <- function(t) t[t$pcode %in% pcodes & t$active %in% acts, c("pcode", "active", "kg_ai_ha"), drop = FALSE]
    hr <- if (length(hs)) take(use_rates$host[use_rates$host$host_grp %in% hs, , drop = FALSE]) else NULL
    pr <- if (length(ps)) {
      t <- use_rates$pest[use_rates$pest$pest_grp %in% ps, , drop = FALSE]
      if (length(hs)) t <- t[t$host_grp %in% hs, , drop = FALSE]
      take(t)
    } else NULL
    if (!is.null(pr) && !is.null(hr)) hr <- hr[!hr$active %in% pr$active, , drop = FALSE]
    rbind(pr, hr)
  }
  ctx_rows <- reactive({
    m <- saved_match(); if (!m$refined) return(NULL)
    use_rate_rows(target_hosts(), target_pests(), m$pcodes, m$actives)
  })
  RATE_SOURCE <- c(label = "Label, this use", user = "Your rate", fallback = "Reference (no label rate for this use)",
                   reference = "Crop-use reference", all_use = "Label maximum, all uses")
  default_rate <- function(a) {
    i   <- match(a, raw_data$PesticideactiveIngredient)
    kg  <- raw_data$MaximumPermittedPesticideApplicationRatekgPAIha[i]
    src <- ifelse(raw_data$has_ref_rate[i], "reference", "all_use")
    r <- ctx_rows()
    if (!is.null(r)) {
      mx  <- tapply(r$kg_ai_ha, r$active, max)
      hit <- a %in% names(mx)
      kg[hit] <- mx[a[hit]]; src[hit] <- "label"
      src[!hit & a %in% saved_match()$actives] <- "fallback"
    }
    data.frame(active = a, kg = as.numeric(kg), source = src, stringsAsFactors = FALSE)
  }
  default_kg <- function(a) default_rate(a)$kg
  product_label_kg <- function(a, pcodes) {
    r <- ctx_rows(); if (is.null(r)) return(rep(NA_real_, length(pcodes)))
    r <- r[r$active == a, , drop = FALSE]
    vapply(pcodes, function(p) { x <- r$kg_ai_ha[r$pcode == p]; if (length(x)) max(x) else NA_real_ }, numeric(1))
  }
  rate_info <- reactive({
    ri <- default_rate(unique(raw_data$PesticideactiveIngredient))
    pl <- product_plan()
    if (!is.null(pl)) {
      u <- pl[!is.na(pl$kg_ha), , drop = FALSE]; h <- match(u$active, ri$active); ok <- !is.na(h)
      ri$kg[h[ok]] <- u$kg_ha[ok]; ri$source[h[ok]] <- "user"
    }
    ri
  })

  output$tgt_match_summary <- renderUI({
    m <- target_match(); n <- length(m$actives); ok <- n > 0
    col <- if (ok) "#2e7d32" else "#c62828"
    msg <- if (!m$refined) "chemicals from Step 3 (no crop or pest chosen)"
           else if (ok) sprintf("of %d chemicals registered for this use (%d products)",
                                length(step3_actives()), length(m$pcodes))
           else "chemicals registered for this combination. Widen the choice to continue."
    div(class = "tgt-banner", style = sprintf("background:%s; border:1px solid %s;", if (ok) "#f1f8e9" else "#ffebee",
                                              if (ok) "#c5e1a5" else "#ffcdd2"),
        div(class = "tgt-count", style = sprintf("color:%s;", col), n),
        div(div(style = sprintf("font-size:13px; color:%s; font-weight:600;", col), msg),
            if (ok && m$refined && n <= 40) div(style = "margin-top:6px;", lapply(m$actives, function(a) tags$span(class = "tgt-chip", a)))))
  })

  conc_factor <- function(cucode) {
    u <- tolower(gsub("\\s", "", as.character(cucode)))
    ifelse(u %in% c("g/l", "g/kg"), 1e-3, ifelse(u == "mg/l", 1e-6, ifelse(u %in% c("kg/l", "kg/kg"), 1, NA_real_)))
  }
  rate_unit <- function(cucode) {
    u <- tolower(gsub("\\s", "", as.character(cucode)))
    ifelse(u %in% c("g/l", "mg/l", "kg/l"), "L/ha", ifelse(u %in% c("g/kg", "kg/kg"), "kg/ha", NA_character_))
  }
  products_for <- function(a, pcodes, refined) {
    d <- raw_data[raw_data$PesticideactiveIngredient == a, c("pcode_int", "fpname", "camount", "cucode")]
    if (refined) d <- d[d$pcode_int %in% pcodes, , drop = FALSE]
    d <- d[!duplicated(d$pcode_int) & !is.na(conc_factor(d$cucode)) & !is.na(suppressWarnings(as.numeric(d$camount))), ]
    d[order(d$fpname), ]
  }
  tp_id <- function(kind, a) paste0("tp_", kind, "_", rate_epoch(), "_", match(a, ALL_ACTIVES))

  output$tgt_products_count <- renderUI({
    m <- saved_match()
    tags$span(style = "font-size:13px; color:#555;",
              sprintf("%d chemicals%s", length(m$actives),
                      if (m$refined) paste0(" registered for ", paste(title_case(c(host_label(target_hosts()), pest_label(target_pests()))), collapse = ", "),
                                            sprintf("; %d with a label rate for this use",
                                                    sum(default_rate(m$actives)$source == "label")))
                      else " from Step 3"))
  })

  output$tgt_products_ui <- renderUI({
    m <- saved_match(); acts <- m$actives
    if (!length(acts)) return(tags$p(style = "color:#c62828;", "No chemicals to show. Go back to Step 4 or Step 3."))
    pl <- isolate(product_plan())
    prev <- function(kind, a, col, dflt) {
      v <- isolate(input[[tp_id(kind, a)]])
      if (!is.null(v) && !(length(v) == 1 && is.na(v))) return(v)
      if (!is.null(pl) && a %in% pl$active) { x <- pl[[col]][match(a, pl$active)]; if (!is.na(x)) return(x) }
      dflt
    }
    rows <- lapply(acts, function(a) {
      d  <- products_for(a, m$pcodes, m$refined); dr <- default_rate(a); dk <- dr$kg
      f  <- suppressWarnings(as.numeric(d$camount)) * conc_factor(d$cucode)
      own <- product_label_kg(a, d$pcode_int)
      units <- c(list(default = list(u = "kg ai/ha", f = 1, d = signif(dk, 3))),
                 setNames(lapply(seq_len(nrow(d)), function(i)
                   list(u = rate_unit(d$cucode[i]), f = f[i],
                        d = signif((if (is.na(own[i])) dk else own[i]) / f[i], 3))), as.character(d$pcode_int)))
      ch <- c("Any product (rate in kg ai/ha)" = "default",
              setNames(as.character(d$pcode_int), paste0(d$fpname, conc_label(d$camount, d$cucode))))
      sel <- as.character(prev("prod", a, "pcode", "default")); if (!sel %in% ch) sel <- "default"
      rv  <- suppressWarnings(as.numeric(prev("rate", a, "product_rate", NA)))
      cu  <- units[[sel]]
      tags$tr(class = "tp-row", `data-units` = as.character(jsonlite::toJSON(units, auto_unbox = TRUE, digits = NA)),
              tags$td(style = "width:24%;",
                      tags$a(href = "#", class = "tp-chem", title = "Chemical details",
                             onclick = sprintf("Shiny.setInputValue('tp_chem_click', %s, {priority: 'event'}); return false;",
                                               jsonlite::toJSON(a, auto_unbox = TRUE)), tags$b(a)),
                      div(class = "tp-sub", n_products(nrow(d)))),
              tags$td(class = "tp-def", style = "width:13%;", sprintf("%s kg ai/ha", signif(dk, 3)),
                      div(class = "tp-sub", style = if (dr$source == "fallback") "color:#c62828;",
                          if (dr$source == "fallback") "No label rate for this use" else RATE_SOURCE[[dr$source]])),
              tags$td(style = "width:41%;", selectizeInput(tp_id("prod", a), NULL, choices = ch, selected = sel, width = "100%"),
                      tags$a(href = "#", class = "tp-label-link",
                             onclick = sprintf("Shiny.setInputValue('tp_label_click', {a: %s, p: $(this).closest('.tp-row').find('select').val(), n: Math.random()}, {priority: 'event'}); return false;",
                                               jsonlite::toJSON(a, auto_unbox = TRUE)), "Label information")),
              tags$td(style = "width:22%;",
                      div(class = "input-group input-group-sm",
                          tags$input(id = tp_id("rate", a), type = "number", class = "form-control", min = "0", step = "any",
                                     placeholder = cu$d, value = if (!is.na(rv)) rv else NULL),
                          tags$span(class = "input-group-addon input-group-text tp-unit", cu$u)),
                      div(class = "tp-conv", if (is.na(rv)) "Default rate used" else sprintf("= %s kg ai/ha", signif(rv * cu$f, 3)))))
    })
    tags$table(class = "tp-table",
               tags$thead(tags$tr(tags$th("Chemical"), tags$th("Default rate"), tags$th("Product"), tags$th("Your rate"))),
               tags$tbody(rows))
  })

  read_plan <- function() {
    m <- saved_match(); acts <- m$actives
    if (!length(acts)) return(NULL)
    do.call(rbind, lapply(acts, function(a) {
      pc <- input[[tp_id("prod", a)]] %||% "default"
      rt <- suppressWarnings(as.numeric(input[[tp_id("rate", a)]] %||% NA))
      if (length(rt) != 1) rt <- NA_real_
      row <- data.frame(active = a, include = TRUE, pcode = NA_integer_, product = NA_character_, conc = NA_character_,
                        unit = "kg ai/ha", product_rate = NA_real_, kg_ha = NA_real_, default_kg_ha = default_kg(a),
                        stringsAsFactors = FALSE)
      ok_rate <- !is.na(rt) && rt > 0
      if (nzchar(pc) && pc != "default") {
        d <- raw_data[raw_data$PesticideactiveIngredient == a & raw_data$pcode_int == as.integer(pc), ][1, ]
        row$pcode <- as.integer(pc); row$product <- d$fpname
        row$conc <- conc_text(d$camount, d$cucode); row$unit <- rate_unit(d$cucode)
        if (ok_rate) { row$product_rate <- rt; row$kg_ha <- rt * suppressWarnings(as.numeric(d$camount)) * conc_factor(d$cucode) }
      } else if (ok_rate) { row$product_rate <- rt; row$kg_ha <- rt }
      row
    }))
  }

  save_targets <- function() {
    m <- target_match()
    if (m$refined && !length(m$actives)) {
      showNotification("No chemicals are registered for this crop and pest combination. Widen the choice or skip.",
                       type = "warning", duration = 6)
      return(FALSE)
    }
    target_hosts(input[[tgt_id("hosts")]] %||% character(0)); target_pests(input[[tgt_id("pests")]] %||% character(0))
    target_actives(if (m$refined) m$actives else NULL)
    TRUE
  }
  save_products <- function() {
    pl <- read_plan()
    if (!is.null(pl)) {
      big <- pl$active[!is.na(pl$kg_ha) & is.finite(pl$default_kg_ha) & pl$kg_ha > 10 * pl$default_kg_ha]
      if (length(big)) showNotification(sprintf("Check the rate for %s: it is more than 10 times the default rate.",
                                                paste(big, collapse = ", ")), type = "warning", duration = 10)
    }
    product_plan(if (is.null(pl) || all(is.na(pl$kg_ha))) NULL else pl[!is.na(pl$kg_ha), , drop = FALSE])
    TRUE
  }
  clear_rates <- function() { product_plan(NULL); rate_epoch(isolate(rate_epoch()) + 1) }
  clear_targets <- function() {
    target_actives(NULL); target_hosts(character(0)); target_pests(character(0)); tgt_show_all(FALSE)
    tgt_epoch(isolate(tgt_epoch()) + 1); clear_rates()
  }
  observeEvent(input$reset_rates, clear_rates())

  selected_chemicals <- reactive({
    chems <- unique(pesticide_selection()$PesticideactiveIngredient)
    ta <- target_actives(); if (!is.null(ta)) chems <- intersect(chems, ta)
    if (length(chems) > 0) chems else NULL
  })

  # ---- Risk table ----
  table_order <- reactiveVal(NULL)
  output$risk_classification_table <- DT::renderDataTable({
    df <- data_mod() %>%
      distinct(PesticideactiveIngredient, .keep_all = TRUE) %>%
      mutate(
        Aquatic_Category = risk_category(ExposureScore, Tox_scaled),
        Aquatic_Default  = risk_category(ExposureScore_default, Tox_scaled),
        Above_Category   = risk_category(Volatilisation_scaled, FlyingEffect_scaled),
        Below_Category   = risk_category(SoilExp_scaled, BelowEffect_scaled),

        Human_Category = case_when(
          WHO_classification %in% c("Ia", "Ib") ~ "High",
          WHO_classification == "II"             ~ "Moderate",
          WHO_classification == "III"            ~ "Low",
          WHO_classification == "U"              ~ "Negligible",
          TRUE                                   ~ "Unclassified"
        )
      ) %>%
      arrange(
        desc(pmax(
          case_when(Aquatic_Category == "Severe" ~ 4, Aquatic_Category == "High" ~ 3,
                    Aquatic_Category == "Moderate" ~ 2, Aquatic_Category == "Low" ~ 1, TRUE ~ 0),
          case_when(Human_Category == "High" ~ 3, Human_Category == "Moderate" ~ 2,
                    Human_Category == "Low" ~ 1, Human_Category == "Negligible" ~ 0, TRUE ~ 0),
          case_when(Above_Category  == "Severe" ~ 4, Above_Category  == "High" ~ 3,
                    Above_Category  == "Moderate" ~ 2, Above_Category  == "Low" ~ 1, TRUE ~ 0),
          case_when(Below_Category  == "Severe" ~ 4, Below_Category  == "High" ~ 3,
                    Below_Category  == "Moderate" ~ 2, Below_Category  == "Low" ~ 1, TRUE ~ 0)
        )),
        PesticideactiveIngredient
      ) %>%
      mutate(`Rate applied (kg ai/ha)` = signif(AppliedRate_umol_ha * Molecularweightgmol / 1e9, 3)) %>%
      mutate(`Rate source` = unname(RATE_SOURCE[rate_info()$source[match(PesticideactiveIngredient, rate_info()$active)]])) %>%
      select(Pesticide = PesticideactiveIngredient,
             `Rate applied (kg ai/ha)`, `Rate source`,
             Aquatic = Aquatic_Category,
             `Aquatic (baseline)` = Aquatic_Default,
             `Pollinator` = Above_Category,
             `Soil fauna` = Below_Category,
             Human = Human_Category)
    table_order(df$Pesticide)
    if (identical(df$Aquatic, df$`Aquatic (baseline)`)) df <- df %>% dplyr::select(-`Aquatic (baseline)`)

    risk_colors <- c(
      "Severe"        = "#800026",
      "High"          = "#E31A1C",
      "Moderate"      = "#FDBF6F",
      "Low"           = "#B2DF8A",
      "Negligible"    = "#A6CEE3",
      "Unclassified"  = "#CCCCCC",
      "No data"       = "#CCCCCC"
    )

    dt <- datatable(
      df,
      rownames = FALSE,
      selection = 'single',
      options = list(
        paging = FALSE,
        dom = 'frti',
        autoWidth = FALSE
      ),
      class = 'risk-summary-table'
    )

    for (col in intersect(c("Aquatic", "Aquatic (baseline)", "Human", "Pollinator", "Soil fauna"), names(df))) {
      dt <- formatStyle(
        dt, columns = col,
        backgroundColor = styleEqual(names(risk_colors), unname(risk_colors)),
        color = styleEqual(names(risk_colors), ifelse(names(risk_colors) %in% c("Severe", "High"), "white", "#333333")),
        fontWeight = "bold"
      )
    }

    dt
  })

  observeEvent(input$risk_classification_table_rows_selected, {
    selected <- input$risk_classification_table_rows_selected
    if (length(selected) == 0) return()

    row_idx <- selected[1]
    ord <- table_order()
    if (is.null(ord) || row_idx > length(ord)) return()
    pesticide_name <- ord[row_idx]
    show_chemical_modal(pesticide_name)
  })

  missing_note <- function(x, y) renderUI({
    d <- data_mod() %>% distinct(PesticideactiveIngredient, .keep_all = TRUE)
    m <- sort(d$PesticideactiveIngredient[is.na(d[[x]]) | is.na(d[[y]])])
    if (length(m)) tags$p(style = "font-size:12px; color:#555;", paste0("Not shown (missing data): ", paste(m, collapse = ", "), "."))
  })
  output$missing_note_aquatic    <- missing_note("ExposureScore", "Tox_scaled")
  output$missing_note_pollinator <- missing_note("Volatilisation_scaled", "FlyingEffect_scaled")
  output$missing_note_soil       <- missing_note("SoilExp_scaled", "BelowEffect_scaled")

  # ---- Baseline and calculation helpers ----
  format_bracketed <- function(txt, fallback = "No information available.") {
    if (is.null(txt) || is.na(txt) || !nzchar(trimws(txt))) return(p(fallback))
    items <- regmatches(txt, gregexpr("(?<=\\[)[^\\]]+(?=\\])", txt, perl = TRUE))[[1]]
    if (length(items) == 0) return(p(txt))
    tags$ul(lapply(items, tags$li))
  }
  benefit_share   <- reactive({
    d <- input$practice_discount_pct
    if (is.null(d) || is.na(d)) adv_init()$practice_benefit_pct %||% PRACTICE_BENEFIT_DEFAULT else 100 - d
  })
  practice_shrink <- reactive(1 - benefit_share() / 100)

  # ---- Management settings ----
  num_or <- function(x, d) if (is.null(x) || length(x) == 0 || is.na(x)) d else x
  default_settings <- function() list(
    rate_reduction = 0, app_date = Sys.Date(), app_method = "Boom", spray_quality = "COARSE", habitat_dist = 2,
    runoff_source = "rainfall", irrig_source = "central", irrig_first_date = Sys.Date() + DEFAULT_IRRIG_FIRST_D,
    irrig_interval = DEFAULT_IRRIG_INTERVAL_D, irrig_depth_mm = DEFAULT_IRRIG_DEPTH_MM,
    Vegetative_Filter_Strips = 0, tillage_type = "none", land_practices = character(0))
  baseline_settings <- reactiveVal(default_settings())
  scenario_init <- reactiveVal(NULL)
  adv_init <- reactiveVal(list(practice_benefit_pct = PRACTICE_BENEFIT_DEFAULT, vapour_dep_pct = VAPOUR_DEP_DEFAULT_PCT))
  read_settings <- function(p) {
    s <- default_settings()
    for (id in names(s)) {
      v <- input[[paste0(p, id)]]
      if (!is.null(v) && length(v) > 0 && !(length(v) == 1 && is.na(v))) s[[id]] <- v
    }
    s$habitat_dist     <- as.numeric(s$habitat_dist)
    s$app_date         <- as.Date(s$app_date)
    s$irrig_first_date <- as.Date(s$irrig_first_date)
    s
  }
  scenario_settings <- reactive({
    if (is.null(input$app_date)) baseline_settings() else read_settings("")
  })
  settings_label <- function(s) {
    paste0(format(s$app_date, "%d %b %Y"), "; ", names(DRIFT_METHODS)[DRIFT_METHODS == s$app_method],
           if (!grepl("^Vertical", s$app_method)) paste0(", ", unname(SPRAY_QUALITY_SHORT[s$spray_quality]),
                                                          " spray") else "",
           ", habitat ", s$habitat_dist, " m; ", irrigation_label(s$runoff_source),
           if (s$runoff_source != "rainfall") sprintf(" (first irrigation %s, every %g days, %g mm)",
                                                      format(s$irrig_first_date, "%d %b"), s$irrig_interval,
                                                      s$irrig_depth_mm) else "",
           "; ", 100 - s$rate_reduction, "% of label rate; tillage ", gsub("_", " ", s$tillage_type),
           if (length(s$land_practices)) paste0("; ", paste(gsub("_", " ", s$land_practices), collapse = ", ")) else "",
           if (s$Vegetative_Filter_Strips > 0) paste0("; ", s$Vegetative_Filter_Strips, " m filter strip") else "")
  }

  for (pp in c("", "b_")) local({
    p <- pp
    observeEvent(input[[paste0(p, "runoff_source")]], {
      rs <- input[[paste0(p, "runoff_source")]]
      if (!is.null(rs) && rs != "rainfall") {
        ch  <- irrigation_source_choices(rs)
        cur <- input[[paste0(p, "irrig_source")]]
        updateSelectInput(session, paste0(p, "irrig_source"), choices = ch,
                          selected = if (!is.null(cur) && cur %in% ch) cur else "central")
      }
    }, ignoreInit = TRUE)
  })

  # ---- Runoff timing ----
  site_cn <- reactive({
    sg <- soilgrids_data(); ha <- hydroatlas_data()
    clay <- if (!is.null(user_clay_pct())) user_clay_pct() else
      if (!is.null(sg) && !is.na(sg$clay_pct)) sg$clay_pct else
      if (!is.null(ha) && !is.na(ha$cly_pc)) ha$cly_pc else NA_real_
    sand <- if (!is.null(user_sand_pct())) user_sand_pct() else
      if (!is.null(sg) && !is.na(sg$sand_pct)) sg$sand_pct else
      if (!is.null(ha) && !is.na(ha$snd_pc)) ha$snd_pc else NA_real_
    g <- hsg_from_texture(clay, sand)
    list(hsg = g, CN = if (!is.null(user_cn())) user_cn() else unname(CN_ROWCROP[g]),
         CN_table = unname(CN_ROWCROP[g]))
  })
  site_tmonth <- reactive({
    ha <- hydroatlas_data()
    t_ann <- if (!is.null(user_temp_c())) user_temp_c() else
      if (!is.null(ha) && !is.na(ha$temp_c)) ha$temp_c else 25
    tm <- if (!is.null(ha)) ha$tmp_monthly else NULL
    if (is.null(tm) || all(is.na(tm))) return(rep(t_ann, 12))
    tm[is.na(tm)] <- mean(tm, na.rm = TRUE)
    if (!is.null(user_temp_c())) tm <- tm + (user_temp_c() - mean(tm))
    tm
  })
  irr_schedule <- function(s) {
    if (is.null(s$runoff_source) || s$runoff_source == "rainfall") return(NULL)
    iv <- irrigation_values(s$runoff_source, s$irrig_source %||% "central")
    list(first_offset = as.numeric(s$irrig_first_date - s$app_date),
         interval = max(1, num_or(s$irrig_interval, DEFAULT_IRRIG_INTERVAL_D)),
         depth_mm = max(0, num_or(s$irrig_depth_mm, DEFAULT_IRRIG_DEPTH_MM)),
         frac = if (is.null(iv)) 0 else iv$runoff_fraction)
  }
  timing_for <- function(s) {
    rain <- rain_daily_data(); if (is.null(rain)) return(NULL)
    cn <- site_cn()
    tl <- runoff_timing_fn(rain, s$app_date, cn$CN, irr_schedule(s), site_tmonth())
    if (is.null(tl)) return(NULL)
    c(tl, list(source = attr(rain, "source") %||% "daily rainfall", hsg = cn$hsg, CN = cn$CN))
  }
  base_timing   <- reactive(timing_for(baseline_settings()))
  scen_timing   <- reactive(timing_for(scenario_settings()))
  b_live_timing <- reactive(timing_for(read_settings("b_")))
  timing_text_for <- function(s, tm) {
    pct <- function(x) sprintf("%.0f%%", 100 * x)
    if (is.null(tm))
      return(sprintf("No daily rainfall: %d days to first runoff assumed.", as.integer(DEFAULT_DAYS_TO_RUNOFF)))
    md <- if (is.finite(tm$median_days)) sprintf("about %d days", as.integer(round(tm$median_days))) else
      sprintf("more than %d days", RUNOFF_HORIZON_D)
    paste0("First runoff (at least ", RUNOFF_EVENT_MM, " mm): ", md, " after application (median of ",
           tm$n_years, " years, ", sub(" daily rainfall.*", "", tm$source), "; curve number ", tm$CN, "). Chance within 7, 14 and ",
           "30 days: ", pct(tm$p7), ", ", pct(tm$p14), ", ", pct(tm$p30), ".",
           if (!is.null(irr_schedule(s))) paste0(" From irrigation in ", pct(tm$share_irr), " of years.") else "")
  }
  timing_text <- reactive(timing_text_for(scenario_settings(), scen_timing()))
  timing_box <- function(txt) tags$div(style = paste0("font-size:11px; color:#444; background:#f7f9fb; ",
                                                      "border:1px solid #dde3ea; border-radius:4px; padding:6px; ",
                                                      "margin-bottom:6px; line-height:1.5;"), txt)
  output$runoff_timing_panel   <- renderUI(timing_box(timing_text()))
  output$b_runoff_timing_panel <- renderUI(timing_box(timing_text_for(read_settings("b_"), b_live_timing())))

  vfs_ratio <- function(w) {
    if (is.null(w) || is.na(w) || w <= 0) return(1)
    widths <- c(0, 5, 10, 20, 30)
    removal_percent <- c(0, 62, 83, 92, 93)
    predicted_removal <- approx(widths, removal_percent, xout = w, rule = 2)$y
    1 - (1 - practice_shrink()) * (predicted_removal / 100)
  }

  coef_means <- function(w) {
    vfs_val <- vfs_ratio(w)
    c(
      "No_Till"      = 1.00,
      "Minimum_Till" = 0.85,
      "Micro_Dams"   = 0.13,
      "Vegetative_Filter_Strips" = vfs_val
    )
  }

  # ---- Scores ----
  powlog_minmax01 <- function(x, power, log_min, log_max) {
    lx      <- log10(pmax(x, 1e-30))
    shifted <- pmax(lx - log_min, 0)
    span    <- log_max - log_min
    pmax(0, pmin(1, (shifted / span)^power))
  }

  AXIS_PCT_LO <- 1
  AXIS_PCT_HI <- 99

  axis_bounds <- local({
    raw_data <- raw_data %>% dplyr::distinct(PesticideactiveIngredient, .keep_all = TRUE)
    pct_range <- function(x) {
      x <- x[is.finite(x)]
      as.numeric(quantile(x, c(AXIS_PCT_LO, AXIS_PCT_HI) / 100, na.rm = TRUE))
    }

    tox_cols <- c("ETVsPC95FreshwaterumolML", "LC5072h", "LC5096h",
                  "LC507day", "Chronic21daysNOEC", "Chronic28daysNOEC")
    tox_ranges <- lapply(tox_cols, function(col) {
      vals <- log10(pmax(1 / pmax(raw_data[[col]], eps), 1e-30))
      pct_range(vals)
    })
    names(tox_ranges) <- tox_cols

    OC_b <- 0.0087; bd_b <- 1.5; ERO_base <- 30; SC_ref <- 5.0
    d <- raw_data %>%
      filter(MaximumPermittedIngredientApplicationRateumolPAIha > 0) %>%
      mutate(
        Kd_b       = pmax(KOCmLg * OC_b, eps),
        C_soil_0_b = MaximumPermittedIngredientApplicationRateumolPAIha / (25 * bd_b * 1000 * 10),
        C_soil_t_b = pmax(
          C_soil_0_b * exp(-log(2) * 21 / pmax(DT50Soil * 2.0, 0.1)),
          C_soil_0_b * 0.05
        ),
        ERO_b      = ERO_fun(Kd_b, SC_ref, C_soil_t_b * Molecularweightgmol / 1000, ERO_base),
        C_RO_b     = ERO_b * C_soil_t_b,
        k_b        = log(2) / pmax(DT50Soil * 2.0, 0.1),
        C_twa_b    = C_soil_0_b * pmax((1 - exp(-k_b * SOIL_TWA_DAYS)) / (k_b * SOIL_TWA_DAYS), 0.05),
        Dose_b     = pollinator_dose(MaximumPermittedPesticideApplicationRatekgPAIha,
                                     VP_ref_Pa, 20,
                                     "Boom", "COARSE", 2, VAPOUR_DEP_DEFAULT_PCT)
      )

    list(
      Tox     = tox_ranges,
      MinLD50 = pct_range(log10(pmax(
        1 / pmax(raw_data$Honeybees_Contact_acute_48hr_LD50_.ug_per_bee., eps), 1e-30))),
      Below   = pct_range(log10(pmax(
        1 / pmax(raw_data$Earthworms_Acute_14d_LC50_.mg.kg., eps), 1e-30))),
      Exposure = pct_range(log10(pmax(d$C_RO_b * Q_REF_MM, 1e-30))),
      Pollinator = pct_range(log10(pmax(d$Dose_b, 1e-30))),
      Csoil    = pct_range(log10(pmax(d$C_twa_b, 1e-30)))
    )
  })

  compute_scores <- function(s, tm) {
    till_sel  <- s$tillage_type %||% "none"
    tillage_p <- if (till_sel != "none") till_sel else character(0)
    other_p   <- s$land_practices %||% character(0)
    effective_practices <- c(tillage_p, other_p)

    if (!is.null(s$Vegetative_Filter_Strips) && s$Vegetative_Filter_Strips > 0) {
      effective_practices <- union(effective_practices, "Vegetative_Filter_Strips")
    }

    base_coef_means <- coef_means(s$Vegetative_Filter_Strips)

    mgmt_coef_means <- base_coef_means

    conservative_mgmt_coefs <- 1 - ((1 - mgmt_coef_means) * (1 - practice_shrink()))

    ha <- hydroatlas_data()

    tmonth      <- site_tmonth()
    temp_app    <- tmonth[as.integer(format(s$app_date, "%m"))]

    slp_site   <- if (!is.null(user_slope_dg())) user_slope_dg() else
      if (!is.null(ha) && !is.na(ha$slp_dg)) ha$slp_dg else 2.0
    slope_mult <- pmax(0.2, pmin(1.0, slp_site / 5))

    sg <- soilgrids_data()

    cly_site <- if (!is.null(user_clay_pct())) user_clay_pct() else
      if (!is.null(sg) && !is.na(sg$clay_pct)) sg$clay_pct else
      if (!is.null(ha) && !is.na(ha$cly_pc))   ha$cly_pc else 25

    sg_snd <- if (!is.null(user_sand_pct())) user_sand_pct() else
      if (!is.null(sg) && !is.na(sg$sand_pct)) sg$sand_pct else
      if (!is.null(ha) && !is.na(ha$snd_pc)) ha$snd_pc else
        pmax(0, 100 - cly_site - 30)
    sg_slt <- if (!is.null(user_silt_pct())) user_silt_pct() else
      if (!is.null(sg) && !is.na(sg$silt_pct)) sg$silt_pct else
      if (!is.null(ha) && !is.na(ha$slt_pc)) ha$slt_pc else
        pmax(0, 100 - cly_site - sg_snd)
    sg_cly <- cly_site

    bd_input <- if (!is.null(user_bd_gcm3())) {
      pmin(pmax(user_bd_gcm3(), 0.9), 1.85)
    } else if (!is.null(sg) && !is.na(sg$bd_gcm3)) {
      pmin(pmax(sg$bd_gcm3, 0.9), 1.85)
    } else {
      bd_base   <- if (cly_site > 35) 1.30 else if (cly_site > 15) 1.45 else 1.60
      bd_om_adj <- switch(sel_om(), "low"=0.08, "moderate"=0.00, "high"=-0.10)
      pmin(pmax(bd_base + bd_om_adj, 0.9), 1.85)
    }

    oc_input <- if (!is.null(user_soc_pct()) && user_soc_pct() > 0) {
      user_soc_pct() / 100
    } else if (!is.null(sg) && !is.na(sg$soc_pct) && sg$soc_pct > 0) {
      sg$soc_pct / 100
    } else {
      switch(sel_om(), "low"=0.005, "moderate"=0.0087, "high"=0.0145)
    }

    infil_mm <- if (cly_site > 35) 30 else if (cly_site > 15) 50 else 80

    is_irrigation     <- s$runoff_source != "rainfall"
    irr_vals          <- irrigation_values(s$runoff_source, s$irrig_source %||% "central")
    irrig_SC          <- if (!is.null(irr_vals)) irr_vals$SC else NA_real_

    filtered_data <- if (is.null(selected_chemicals()) || length(selected_chemicals()) == 0) {
      raw_data
    } else {
      raw_data %>% dplyr::filter(PesticideactiveIngredient %in% selected_chemicals())
    }

    ri  <- rate_info()
    ov  <- ri[ri$source %in% c("label", "user", "fallback"), , drop = FALSE]
    hit <- match(filtered_data$PesticideactiveIngredient, ov$active)
    k   <- !is.na(hit) & is.finite(filtered_data$Molecularweightgmol)
    if (any(k)) {
      filtered_data$MaximumPermittedPesticideApplicationRatekgPAIha[k]    <- ov$kg[hit[k]]
      filtered_data$MaximumPermittedIngredientApplicationRateumolPAIha[k] <- ov$kg[hit[k]] * 1e9 /
        filtered_data$Molecularweightgmol[k]
    }
    rate_scalar <- 1 - (s$rate_reduction %||% 0) / 100

    working_data <- filtered_data %>%
      mutate(EffectiveRate_umol_ha = MaximumPermittedIngredientApplicationRateumolPAIha * rate_scalar)

    vol_p <- intersect(effective_practices,
                       c("No_Till", "Minimum_Till", "Micro_Dams"))
    runoff_volume_mult <- if (length(vol_p) == 0) 1 else prod(conservative_mgmt_coefs[vol_p], na.rm = TRUE)

    combined_practice_mult <- max(PRACTICE_MULT_FLOOR, runoff_volume_mult * vfs_ratio(s$Vegetative_Filter_Strips))

    ha_rf <- {
      if (!is.null(user_runoff_frac())) user_runoff_frac()
      else {
        ha <- hydroatlas_data()
        if (!is.null(ha) && !is.na(ha$runoff_frac)) ha$runoff_frac else NA_real_
      }
    }

    SC_rainfall <- {
      rd <- rusle_data()
      if (!is.null(rd) && !is.na(rd$LS_factor) && !is.na(rd$C_factor)) {
        sand_f  <- pmax(0, pmin(100, sg_snd)) / 100
        silt_f  <- pmax(0, pmin(100, sg_slt)) / 100
        clay_f  <- pmax(0, pmin(100, sg_cly)) / 100
        OC_pct  <- pmax(0, oc_input) * 100
        K_f1    <- 0.2 + 0.3 * exp(-0.0256 * sand_f * 100 * (1 - silt_f))
        K_f2    <- if ((clay_f + silt_f) > 0) (silt_f / (clay_f + silt_f))^0.3 else 0
        K_f3    <- if (OC_pct > 0) 1 - 0.25 * OC_pct / (OC_pct + exp(3.72 - 2.95 * OC_pct)) else 1
        K_f4    <- 1 - 0.7 * (1 - sand_f) / ((1 - sand_f) + exp(-5.51 + 22.9 * (1 - sand_f)))
        K_williams <- K_f1 * K_f2 * K_f3 * K_f4 * 0.1317
        K_SI    <- pmax(0.001, pmin(0.15, K_williams))
        R_use   <- if (!is.null(ha) && !is.na(ha$R_factor)) ha$R_factor else
          if (!is.na(rd$R_factor)) rd$R_factor else 200
        A_rusle <- R_use * K_SI * rd$LS_factor * rd$C_factor
        run_mm  <- if (!is.null(ha) && !is.na(ha$run_mm) && ha$run_mm > 0) ha$run_mm else 150
        pmin(SC_RAINFALL_MAX, pmax(SC_RAINFALL_MIN, A_rusle / run_mm * 100))
      } else
        pmax(SC_RAINFALL_MIN, slope_mult * 5)
    }
    runoff_frac_rainfall <- if (!is.na(ha_rf)) ha_rf else pmin(0.95, 0.5 * slope_mult)

    irr_s       <- irr_schedule(s)
    Q_irr_event <- if (!is.null(irr_s)) irr_s$depth_mm * irr_s$frac else 0
    have_clim   <- !is.null(tm)
    eq_days <- function(n) cumsum(q10_day(tmonth[as.integer(format(s$app_date + seq_len(n), "%m"))]))
    ev <- if (have_clim) tm$per_year[, c("days", "E", "infil", "Q_rain", "Q_irr")] else {
      E21 <- eq_days(DEFAULT_DAYS_TO_RUNOFF)[DEFAULT_DAYS_TO_RUNOFF]
      d1 <- if (is_irrigation && is.finite(irr_s$first_offset)) {
        fo <- irr_s$first_offset
        as.integer(round(if (fo >= 1) fo else fo + ceiling((1 - fo) / irr_s$interval) * irr_s$interval))
      } else NA_integer_
      if (is_irrigation && Q_irr_event >= RUNOFF_EVENT_MM && !is.na(d1) && d1 <= RUNOFF_HORIZON_D) {
        data.frame(days = d1, E = eq_days(d1)[d1], infil = infil_mm * d1 / DEFAULT_DAYS_TO_RUNOFF, Q_rain = 0, Q_irr = Q_irr_event)
      } else
        data.frame(days = DEFAULT_DAYS_TO_RUNOFF, E = E21, infil = infil_mm,
                   Q_rain = Q_REF_MM * runoff_frac_rainfall / 0.25, Q_irr = 0)
    }
    E14   <- if (have_clim) tm$E14 else matrix(c(0, eq_days(SOIL_TWA_DAYS)), nrow = 1)
    INF14 <- if (have_clim) tm$INF14 else matrix(infil_mm * (0:SOIL_TWA_DAYS) / DEFAULT_DAYS_TO_RUNOFF, nrow = 1)

    out <- working_data %>%
      mutate(

        AppliedRate_umol_ha = EffectiveRate_umol_ha,

        C_soil_0 = AppliedRate_umol_ha / (25 * bd_input * 1000 * 10),

        Kd = pmax(KOCmLg * oc_input, eps),

        DT50_adj = pmax(DT50Soil * 2.0, 0.1),
        C_soil_t = soil_twa(C_soil_0, DT50_adj, Kd, bd_input, E14, INF14),

        Exposure_index = first_event_loss(C_soil_0, DT50_adj, Kd, bd_input, Molecularweightgmol,
                                          SC_rainfall, if (is_irrigation) irrig_SC else NA_real_, ev) *
          combined_practice_mult,

        ExposureScore = powlog_minmax01(Exposure_index, 1,
                                        axis_bounds$Exposure[1], axis_bounds$Exposure[2]),
        Tox_scaled    = powlog_minmax01(1 / pmax(.data[["ETVsPC95FreshwaterumolML"]], eps), 1,
                                        axis_bounds$Tox[["ETVsPC95FreshwaterumolML"]][1],
                                        axis_bounds$Tox[["ETVsPC95FreshwaterumolML"]][2]),

        Bee_LD50_ug_bee        = pmax(Honeybees_Contact_acute_48hr_LD50_.ug_per_bee., eps),
        FlyingEffect_scaled   = powlog_minmax01(1 / Bee_LD50_ug_bee, 1,
                                               axis_bounds$MinLD50[1], axis_bounds$MinLD50[2]),
        Pollinator_dose_ugbee = pollinator_dose(MaximumPermittedPesticideApplicationRatekgPAIha * rate_scalar,
                                                VP_ref_Pa, temp_app,
                                                s$app_method, s$spray_quality, s$habitat_dist,
                                                num_or(input$vapour_dep_pct, VAPOUR_DEP_DEFAULT_PCT)),
        Volatilisation_scaled = powlog_minmax01(Pollinator_dose_ugbee, 1,
                                               axis_bounds$Pollinator[1], axis_bounds$Pollinator[2]),

        BelowEffect_scaled = powlog_minmax01(1 / pmax(Earthworms_Acute_14d_LC50_.mg.kg., eps), 1,
                                            axis_bounds$Below[1], axis_bounds$Below[2]),
        SoilExp_scaled     = powlog_minmax01(C_soil_t, 1,
                                            axis_bounds$Csoil[1], axis_bounds$Csoil[2])
      )

    out
  }

  data_mod <- reactive({
    out  <- compute_scores(scenario_settings(), scen_timing())
    base <- compute_scores(baseline_settings(), base_timing())
    out$ExposureScore_default  <- base$ExposureScore
    out$Volatilisation_default <- base$Volatilisation_scaled
    out$SoilExp_default        <- base$SoilExp_scaled
    out
  })

  # ---- Risk matrix grid ----
  grid_df <- tibble(
    x_min = rep(c(0.0, 0.33, 0.66), each = 3),
    x_max = rep(c(0.33, 0.66, 1.05), each = 3),
    y_min = rep(c(0.0, 0.33, 0.66), times = 3),
    y_max = rep(c(0.33, 0.66, 1.05), times = 3),
    risk_level = c(
      "Negligible", "Low", "Moderate",
      "Low", "Moderate", "High",
      "Moderate", "High", "Severe"
    )
  )

  # ---- Risk plots ----
  RISK_LEVELS <- c("Negligible", "Low", "Moderate", "High", "Severe")
  RISK_PAL    <- c(Negligible = "#A6CEE3", Low = "#B2DF8A", Moderate = "#FDBF6F", High = "#E31A1C", Severe = "#800026")
  RISK_PLOTS  <- list(
    aquatic    = list(x = "ExposureScore", y = "Tox_scaled", x0 = "ExposureScore_default",
                      xlab = "Aquatic Exposure Score", ylab = "Aquatic Hazard Score"),
    pollinator = list(x = "Volatilisation_scaled", y = "FlyingEffect_scaled", x0 = "Volatilisation_default",
                      xlab = "Pollinator Exposure Score", ylab = "Pollinator Hazard Score"),
    soil       = list(x = "SoilExp_scaled", y = "BelowEffect_scaled", x0 = "SoilExp_default",
                      xlab = "Soil Fauna Exposure Score", ylab = "Soil Fauna Hazard Score"))
  baseline_arrows <- function(pd, x0, x1, y) {
    pa <- pd[!is.na(pd[[x0]]) & !is.na(pd[[x1]]) & abs(pd[[x0]] - pd[[x1]]) > 1e-6, , drop = FALSE]
    if (nrow(pa) == 0) return(NULL)
    list(
      geom_point(data = pa, aes(x = .data[[x0]], y = .data[[y]]), shape = 21, size = 3, colour = "grey35",
                 fill = "white"),
      geom_segment(data = pa, aes(x = .data[[x0]], xend = .data[[x1]], y = .data[[y]], yend = .data[[y]]),
                   colour = "grey20", linewidth = 0.4,
                   arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")))
  }
  risk_plot <- function(kind, dm, base_size = 14, label_size = 4, note = FALSE) {
    sp <- RISK_PLOTS[[kind]]
    pd <- dm %>% distinct(PesticideactiveIngredient, .keep_all = TRUE)
    pd <- pd[!is.na(pd[[sp$x]]) & !is.na(pd[[sp$y]]), , drop = FALSE]
    if (nrow(pd) == 0) return(ggplot() + annotate("text", x = 0, y = 0, label = "No data", size = 6) + theme_void())
    arrows <- baseline_arrows(pd, sp$x0, sp$x, sp$y)
    ggplot() +
      geom_rect(data = grid_df, aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max,
                                    fill = factor(risk_level, levels = RISK_LEVELS)), alpha = 0.4, colour = "grey70") +
      scale_fill_manual(values = RISK_PAL, name = "Risk level") +
      arrows +
      geom_point(data = pd, aes(x = .data[[sp$x]], y = .data[[sp$y]]), size = 3, colour = "black") +
      ggrepel::geom_text_repel(data = pd, aes(x = .data[[sp$x]], y = .data[[sp$y]], label = PesticideactiveIngredient),
                               size = label_size, box.padding = 0.4, point.padding = 0.2, max.overlaps = 30,
                               segment.size = 0.2) +
      labs(x = sp$xlab, y = sp$ylab, caption = if (note && !is.null(arrows)) ARROW_NOTE) +
      theme_pubr(base_size = base_size) + LEGEND_BIG +
      theme(plot.caption = element_text(hjust = 0, size = 9))
  }
  output$pesticidePlot         <- renderPlot(risk_plot("aquatic", data_mod(), base_size = 18, label_size = 5))
  output$terrestrial_risk_plot <- renderPlot(risk_plot("pollinator", data_mod()))
  output$below_risk_plot       <- renderPlot(risk_plot("soil", data_mod()))

  # ---- Pop-ups ----
  observeEvent(input$tp_chem_click, { req(input$tp_chem_click); show_chemical_modal(input$tp_chem_click) })
  no_em <- function(x) if (is.character(x)) gsub("\\s*\u2014\\s*", ", ", x) else x

  observeEvent(input$tp_label_click, {
    cl <- input$tp_label_click; req(cl$a)
    showModal(label_info_modal(cl$a, as.character(cl$p %||% "default")))
  })
  label_info_modal <- function(a, sel = "default") {
    m <- saved_match(); d <- products_for(a, m$pcodes, m$refined)
    f <- suppressWarnings(as.numeric(d$camount)) * conc_factor(d$cucode)
    use_txt <- if (m$refined) paste(title_case(c(host_label(target_hosts()), pest_label(target_pests()))), collapse = ", ") else NULL
    kg <- if (!is.null(ctx_rows())) product_label_kg(a, d$pcode_int) else {
      h <- if (!is.null(use_rates)) use_rates$host[use_rates$host$active == a, , drop = FALSE] else NULL
      vapply(d$pcode_int, function(p) { x <- h$kg_ai_ha[h$pcode == p]; if (length(x)) max(x) else NA_real_ }, numeric(1))
    }
    lr <- use_rates$rows
    lr <- if (is.null(lr)) NULL else lr[lr$active == a & lr$pcode %in% d$pcode_int &
                                          (if (length(target_hosts())) lr$host_grp %in% target_hosts() else TRUE), , drop = FALSE]
    lines_for <- function(pc) {
      x <- if (is.null(lr)) NULL else lr[lr$pcode == pc, , drop = FALSE]
      if (is.null(x) || !nrow(x)) return(NULL)
      x <- x[order(-x$kg_ai_ha), , drop = FALSE]; x <- x[!duplicated(paste(x$host, x$rate_text)), , drop = FALSE]
      tags$ul(class = "lab-lines", lapply(seq_len(min(nrow(x), 6)), function(j) tags$li(
        tags$b(title_case(as.character(x$host[j]))), ": ", x$rate_text[j],
        sprintf(" = %s kg ai/ha", signif(x$kg_ai_ha[j], 3)),
        if (isTRUE(x$per_row[j])) tags$span(class = "lab-note", sprintf(" (per row length, at %g m between rows)", x$row_spacing_m[j])))),
        if (nrow(x) > 6) tags$li(class = "lab-note", sprintf("and %d more lines on the label", nrow(x) - 6)))
    }
    ord <- order(d$pcode_int != suppressWarnings(as.integer(sel)), is.na(kg), d$fpname)
    rows <- lapply(ord, function(i) {
      tagList(
        tags$tr(class = if (identical(as.character(d$pcode_int[i]), sel)) "lab-sel",
                tags$td(d$fpname[i]), tags$td(conc_text(d$camount[i], d$cucode[i])),
                tags$td(if (is.na(kg[i])) "No rate read for this use" else sprintf("%s %s", signif(kg[i] / f[i], 3), rate_unit(d$cucode[i]))),
                tags$td(if (is.na(kg[i])) "" else signif(kg[i], 3)),
                tags$td(tags$a(href = sprintf("https://elabels.apvma.gov.au/%sELBL.pdf", d$pcode_int[i]), target = "_blank", "Label (PDF)"))),
        { li <- lines_for(d$pcode_int[i]); if (!is.null(li)) tags$tr(class = "lab-sub", tags$td(colspan = 5, li)) })
    })
    modalDialog(
      title = h3(a, style = "color:#0173b2; margin-top:0;"), size = "l", easyClose = TRUE,
      tags$p(style = "font-size:13px; color:#444;",
             if (!is.null(use_txt)) paste0("Highest label rate on each product label for ", use_txt, ".")
             else "Highest label rate on each product label for any crop or situation. Choose a crop in Step 4 for rates for that use."),
      tags$table(class = "lab-table",
                 tags$thead(tags$tr(tags$th("Product"), tags$th("Concentration"), tags$th("Label rate (product)"),
                                    tags$th("kg ai/ha"), tags$th("Label"))),
                 tags$tbody(rows)),
      tags$p(class = "src-line", style = "margin-top:10px;",
             "Rates were read automatically from the labels (APVMA 2026b) and may contain errors; rates per length of row are converted ",
             "with a typical row spacing for the crop. Always check the label before use."),
      ref_list("apvma_labels"),
      footer = modalButton("Close"))
  }

  show_chemical_modal <- function(pesticide_name) {
    pest_info <- raw_data %>% filter(PesticideactiveIngredient == pesticide_name) %>% slice(1)
    if (nrow(pest_info) == 0) return(invisible())
    sm <- tryCatch(saved_match(), error = function(e) NULL)
    prod <- raw_data %>% filter(PesticideactiveIngredient == pesticide_name)
    if (!is.null(sm) && isTRUE(sm$refined)) prod <- prod %>% filter(pcode_int %in% sm$pcodes)
    prod <- prod %>% distinct(fpname, .keep_all = TRUE) %>% arrange(fpname)
    n_prod <- nrow(prod); shown <- head(prod, 10)
    has <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(trimws(x)) && !grepl("^no data$", trimws(x), ignore.case = TRUE)
    head4 <- function(txt) h4(txt, style = "color: #0173b2; margin-top: 15px; margin-bottom: 5px;")
    row <- function(lab, val) if (has(val)) tags$tr(tags$td(strong(lab), style = "padding: 6px; width: 150px;"),
                                                    tags$td(val, style = "padding: 6px;"))
    showModal(modalDialog(
      title = h3(pest_info$PesticideactiveIngredient, style = "color: #0173b2; margin-top: 0;"),
      size = "l", easyClose = TRUE,
      div(style = "background-color: #f0f8ff; padding: 15px; margin-bottom: 10px; border-radius: 4px;",
          h4("Chemical information", style = "color: #0173b2; margin-top: 0;"),
          tags$table(row("Class:", pest_info$Class), row("Chemical group:", no_em(pest_info$Chemicalgroup)),
                     row("WHO class:", pest_info$WHO_classification %||% NA),
                     row("Aquatic guideline:", if (has(pest_info$ETV_source %||% NA) && !is.na(pest_info$ETVsPC95FreshwaterugL))
                       sprintf("%s \u00b5g/L (%s)", signif(pest_info$ETVsPC95FreshwaterugL, 3),
                               if (grepl("^ANZG", pest_info$ETV_source)) "ANZG default guideline value" else
                                 "estimated from the lowest acute EC50 or LC50") else NA))),
      if (has(pest_info$Description)) tagList(head4("Description"), p(no_em(pest_info$Description), style = "line-height: 1.6;")),
      head4(if (!is.null(sm) && isTRUE(sm$refined)) "Products registered for your use" else "Registered products"),
      if (n_prod == 0) p(em("No products listed.")) else
        tagList(tags$ul(lapply(seq_len(nrow(shown)), function(i)
                  tags$li(tags$strong(shown$fpname[i]), conc_label(shown$camount[i], shown$cucode[i])))),
                if (n_prod > 10) p(sprintf("and %d more.", n_prod - 10), style = "color:#555;")),
      if (has(pest_info$Example_pests_controlled)) tagList(head4("Pests and weeds controlled"),
                                                         p(no_em(pest_info$Example_pests_controlled), style = "line-height: 1.6;")),
      if (has(pest_info$Application_examples)) tagList(head4("Application examples"),
                                                     p(no_em(pest_info$Application_examples), style = "line-height: 1.6;")),
      if (has(pest_info$Handling_issues)) tagList(head4("Handling and safety"), format_bracketed(no_em(pest_info$Handling_issues))),
      if (has(pest_info$General_human_health_issues)) tagList(head4("Human health"),
                                                            format_bracketed(no_em(pest_info$General_human_health_issues))),
      footer = modalButton("Close")
    ))
  }

  # ---- Human health plot ----
  make_human_matrix <- function(dm) {
    base <- dm %>%
      distinct(PesticideactiveIngredient, .keep_all = TRUE) %>%
      select(Pesticide       = PesticideactiveIngredient,
             WHO_classification,
             Carcinogen      = Carcinogen.,
             Endocrine       = Endocrine_distrupter.,
             Reproduction    = Reproduction.development_effects.,
             Neurotoxicant   = Neurotoxicant.,
             Acetylcholinesterase = Acetyl_cholinesterase_inhibitor.,
             Respiratory     = Respiratory_tract_irritant.,
             `Skin Irritant` = Skin_irritant.,
             `Skin Sensitiser` = Skin_sensitiser.,
             `Eye Irritant`  = Eye_irritant.) %>%
      filter(if_any(Carcinogen:`Eye Irritant`, ~ !is.na(.x) & .x != "NA"))

    if (nrow(base) == 0)
      return(ggplot() + annotate("text", x=0, y=0, label="No data", size=6) + theme_void())

    cond_data <- base %>%
      pivot_longer(cols = Carcinogen:`Eye Irritant`,
                   names_to = "Condition", values_to = "Status") %>%
      mutate(
        Status = case_when(
          Status == "Yes"      ~ "Yes",
          Status == "Possible" ~ "Possible",
          Status == "No"       ~ "No",
          TRUE                 ~ "Unknown"
        ),
        Condition = fct_relevel(Condition,
                                "Carcinogen", "Endocrine", "Reproduction", "Neurotoxicant",
                                "Acetylcholinesterase", "Respiratory", "Skin Irritant",
                                "Skin Sensitiser", "Eye Irritant")
      )

    who_rank_map <- c("Ia"=1L,"Ib"=2L,"II"=3L,"III"=4L,"U"=5L)
    pesticide_order <- base %>%
      mutate(
        yes_count = rowSums(across(Carcinogen:`Eye Irritant`, ~ .x == "Yes"), na.rm = TRUE),
        who_rank  = coalesce(who_rank_map[WHO_classification], 6L)
      ) %>%
      arrange(desc(yes_count), who_rank, Pesticide) %>%
      pull(Pesticide)

    cond_data <- cond_data %>%
      mutate(Pesticide = factor(Pesticide, levels = pesticide_order))

    who_data <- base %>%
      mutate(
        Pesticide = factor(Pesticide, levels = pesticide_order),
        WHO_label = case_when(
          WHO_classification == "Ia"  ~ "Ia",
          WHO_classification == "Ib"  ~ "Ib",
          WHO_classification == "II"  ~ "II",
          WHO_classification == "III" ~ "III",
          WHO_classification == "U"   ~ "U",
          TRUE                        ~ "Not Listed"
        ),
        WHO_label = factor(WHO_label, levels = c("Ia","Ib","II","III","U","Not Listed")),
        x_col = "WHO"
      )

    p_who <- ggplot(who_data, aes(x = x_col, y = Pesticide, fill = WHO_label)) +
      geom_tile(color = "white", linewidth = 1) +
      scale_fill_manual(
        values = c("Ia"="#800026","Ib"="#E31A1C","II"="#FDBF6F",
                   "III"="#B2DF8A","U"="#A6CEE3","Not Listed"="#CCCCCC"),
        labels = c("Ia"="Ia: Extremely hazardous","Ib"="Ib: Highly hazardous",
                   "II"="II: Moderately hazardous","III"="III: Slightly hazardous",
                   "U"="U: Unlikely to present hazard","Not Listed"="Not listed"),
        drop = FALSE, name = "WHO class"
      ) +
      labs(x = NULL, y = NULL) +
      theme_pubr(base_size = 12) + LEGEND_BIG +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(), legend.position = "right")

    p_cond <- ggplot(cond_data, aes(x = Condition, y = Pesticide, fill = Status)) +
      geom_tile(color = "white", linewidth = 1) +
      scale_fill_manual(
        values = c("Yes"="#E31A1C","Possible"="#FDBF6F","No"="#A6CEE3","Unknown"="#CCCCCC"),
        breaks = c("Yes", "Possible", "No", "Unknown"), name = "Associated"
      ) +
      labs(x = "Health condition", y = "Pesticide",
           title = "Human health hazard profile") +
      theme_pubr(base_size = 12) + LEGEND_BIG +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            panel.grid = element_blank(), legend.position = "right")

    p_cond + p_who + patchwork::plot_layout(widths = c(9, 1))
  }
  output$HumanMatrix <- renderPlot(make_human_matrix(data_mod()))

  # ---- Downloads ----
  png_download <- function(kind, stem) downloadHandler(
    filename = function() paste0(stem, "_", Sys.Date(), ".png"),
    content  = function(file) ggsave(file, plot = risk_plot(kind, data_mod(), base_size = 16, label_size = 4.5, note = TRUE),
                                     width = 8, height = 8, dpi = 600, bg = "white"))
  output$download_aquatic_plot <- png_download("aquatic", "aquatic_risk")
  output$download_above_plot   <- png_download("pollinator", "pollinator_risk")
  output$download_below_plot   <- png_download("soil", "soil_fauna_risk")

  output$download_human_plot <- downloadHandler(
    filename = function() paste0("human_risk_", Sys.Date(), ".png"),
    content = function(file) {
      n_pest <- dplyr::n_distinct(data_mod()$PesticideactiveIngredient)
      h <- min(49, max(4, 1.5 + 0.3 * n_pest))
      ggsave(file, plot = make_human_matrix(data_mod()), width = 11,
             height = h, dpi = if (h > 20) 300 else 600, bg = "white")
    }
  )

  # ---- PDF report ----
  output$download_pdf_report <- downloadHandler(
    filename    = function() paste0("pesticide_risk_report_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content = function(file) {
      sg     <- soilgrids_data()
      ha     <- hydroatlas_data()
      coords <- map_coords()

      # ---- Helpers ----
      risk_cats <- function(d) {
        d %>%
          distinct(PesticideactiveIngredient, .keep_all = TRUE) %>%
          mutate(
            Aquatic = risk_category(ExposureScore, Tox_scaled),
            Human = case_when(
              WHO_classification %in% c("Ia", "Ib") ~ "High",
              WHO_classification == "II"             ~ "Moderate",
              WHO_classification == "III"            ~ "Low",
              WHO_classification == "U"              ~ "Negligible",
              TRUE                                   ~ "Unclassified"),
            `Pollinator` = risk_category(Volatilisation_scaled, FlyingEffect_scaled),
            `Soil fauna` = risk_category(SoilExp_scaled, BelowEffect_scaled)
          ) %>%
          mutate(`Rate (kg ai/ha)` = signif(AppliedRate_umol_ha * Molecularweightgmol / 1e9, 3),
                 `Rate source` = unname(RATE_SOURCE[rate_info()$source[match(PesticideactiveIngredient, rate_info()$active)]])) %>%
          select(Pesticide = PesticideactiveIngredient, `Rate (kg ai/ha)`, `Rate source`, Aquatic,
                 `Pollinator`, `Soil fauna`, Human)
      }
      risk_df  <- risk_cats(data_mod())

      COL_BLUE  <- "#0173b2"
      COL_GREEN <- "#029e73"
      COL_DARK  <- "#2c3e50"
      RISK_FILL <- c("Severe"="#800026","High"="#E31A1C","Moderate"="#FDBF6F",
                     "Low"="#B2DF8A","Negligible"="#A6CEE3","No data"="#CCCCCC","Unclassified"="#CCCCCC")
      RISK_TCOL <- c("Severe"="white","High"="white","Moderate"="#333333",
                     "Low"="#333333","Negligible"="#333333","No data"="#333333","Unclassified"="#333333")
      RPT_DATE  <- format(Sys.Date(), "%d %B %Y")
      ROWS_PP   <- 28L
      n_tab_pg  <- max(1L, as.integer(ceiling(nrow(risk_df) / ROWS_PP)))
      TOTAL_PG  <- 3L + n_tab_pg + 4L + 1L

      rpt_ttheme <- function(hdr_fill, stripe = "#f7f9fc", fs = 9, pad = c(4, 3)) {
        gridExtra::ttheme_default(
          colhead = list(
            bg_params = list(fill = hdr_fill, col = NA),
            fg_params = list(col = "white", fontface = "bold", fontsize = fs,
                             hjust = 0, x = 0.05)
          ),
          core = list(
            bg_params = list(fill = rep_len(c(stripe, "white"), 400), col = "#e5e5e5", lwd = 0.3),
            fg_params = list(fontsize = fs, hjust = 0, x = 0.05)
          ),
          padding = unit(pad, "mm")
        )
      }

      draw_chrome <- function(section_title, page_n) {
        grid.newpage()
        grid.rect(x = 0, y = 1, width = 1, height = 0.075,
                  just = c("left", "top"), gp = gpar(fill = COL_BLUE, col = NA))
        grid.text("PESTICIDE RISK ASSESSMENT REPORT",
                  x = 0.04, y = 0.977, just = c("left", "centre"),
                  gp = gpar(fontsize = 6.5, col = "#9bcfe8"))
        grid.text(section_title, x = 0.04, y = 0.957,
                  just = c("left", "centre"),
                  gp = gpar(fontsize = 12, fontface = "bold", col = "white"))
        grid.text(sprintf("Page %d of %d", page_n, TOTAL_PG),
                  x = 0.96, y = 0.966, just = c("right", "centre"),
                  gp = gpar(fontsize = 9, col = "white"))
        grid.rect(x = 0, y = 0.925, width = 1, height = 0.003,
                  just = c("left", "top"), gp = gpar(fill = COL_GREEN, col = NA))
        grid.rect(x = 0, y = 0, width = 1, height = 0.045,
                  just = c("left", "bottom"), gp = gpar(fill = "#f4f6f8", col = NA))
        grid.segments(0, 0.045, 1, 0.045, gp = gpar(col = "#d8d8d8", lwd = 0.4))
        grid.text(paste0("Pesticide Risk Assessment Tool, version ", APP_VERSION, "  |  ", RPT_DATE),
                  x = 0.04, y = 0.022, just = c("left", "centre"),
                  gp = gpar(fontsize = 7.5, col = "#888888"))
        grid.text("For informational purposes only",
                  x = 0.96, y = 0.022, just = c("right", "centre"),
                  gp = gpar(fontsize = 7, col = "#aaaaaa", fontface = "italic"))
        pushViewport(viewport(x = 0.5, y = 0.4835,
                              width = 0.90, height = 0.850,
                              just = c("centre", "centre")))
      }

      colour_risk_grob <- function(tg, df) {
        for (j in seq_len(ncol(df))[-1]) {
          for (i in seq_len(nrow(df))) {
            val <- as.character(df[i, j])
            if (!is.na(val) && val %in% names(RISK_FILL)) {
              bg_ix <- which(tg$layout$name == "core-bg" & tg$layout$t == i + 1 & tg$layout$l == j)
              if (length(bg_ix)) tg$grobs[[bg_ix]]$gp$fill <- RISK_FILL[[val]]
              fg_ix <- which(tg$layout$name == "core-fg" & tg$layout$t == i + 1 & tg$layout$l == j)
              if (length(fg_ix)) tg$grobs[[fg_ix]]$gp$col  <- RISK_TCOL[[val]]
            }
          }
        }
        tg
      }

      fmt <- function(x, f) if (!is.null(x) && !is.na(x)) sprintf(f, x) else "N/A"

      white_bg <- theme(plot.background  = element_rect(fill = "white", colour = NA),
                        panel.background = element_rect(fill = "white", colour = NA))
      pdf_plot <- function(kind) risk_plot(kind, data_mod(), base_size = 12, label_size = 3.5, note = TRUE) + white_bg
      p_aquatic <- pdf_plot("aquatic"); p_terrestrial <- pdf_plot("pollinator"); p_below <- pdf_plot("soil")

      p_human <- {
        pd <- data_mod() %>%
          distinct(PesticideactiveIngredient, .keep_all = TRUE) %>%
          select(Pesticide = PesticideactiveIngredient,
                 Carcinogen = Carcinogen., Endocrine = Endocrine_distrupter.,
                 Reproduction = Reproduction.development_effects.,
                 Neurotoxicant = Neurotoxicant.,
                 Acetylcholinesterase = Acetyl_cholinesterase_inhibitor.,
                 Respiratory = Respiratory_tract_irritant.,
                 `Skin Irritant` = Skin_irritant., `Skin Sensitiser` = Skin_sensitiser.,
                 `Eye Irritant` = Eye_irritant.) %>%
          tidyr::pivot_longer(cols = -Pesticide, names_to = "Condition", values_to = "Status") %>%
          mutate(
            Status    = case_when(Status == "Yes" ~ "Yes", Status == "Possible" ~ "Possible", Status == "No" ~ "No", TRUE ~ "Unknown"),
            Condition = fct_relevel(Condition, "Carcinogen","Endocrine","Reproduction",
                                    "Neurotoxicant","Acetylcholinesterase","Respiratory",
                                    "Skin Irritant","Skin Sensitiser","Eye Irritant")
          )
        ord <- pd %>% group_by(Pesticide) %>%
          summarise(n = sum(Status == "Yes"), .groups = "drop") %>% arrange(desc(n), Pesticide) %>% pull(Pesticide)
        pd <- pd %>% mutate(Pesticide = factor(Pesticide, levels = ord))
        ggplot(pd, aes(x = Condition, y = Pesticide, fill = Status)) +
          geom_tile(color = "white", linewidth = 0.8) +
          scale_fill_manual(values = c("Yes"="#E31A1C","Possible"="#FDBF6F","No"="#A6CEE3","Unknown"="#CCCCCC"),
                        breaks = c("Yes", "Possible", "No", "Unknown")) +
          labs(x = "Health condition", y = NULL, fill = "Associated") +
          theme_pubr(base_size = 11) + LEGEND_BIG +
          theme(axis.text.x    = element_text(angle = 45, hjust = 1),
                panel.grid     = element_blank(),
                legend.position = "right",
                plot.background  = element_rect(fill = "white", colour = NA),
                panel.background = element_rect(fill = "white", colour = NA))
      }

      if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(file, width = 8.27, height = 11.69, onefile = TRUE) else
        pdf(file, width = 8.27, height = 11.69)
      dv <- dev.cur()
      on.exit(if (dv %in% dev.list()) dev.off(dv), add = TRUE)

      # ---- Page 1: cover ----
      grid.newpage()
      grid.rect(x = 0, y = 1, width = 1, height = 0.42,
                just = c("left","top"), gp = gpar(fill = COL_BLUE, col = NA))
      grid.rect(x = 0, y = 0.58, width = 1, height = 0.010,
                just = c("left","top"), gp = gpar(fill = COL_GREEN, col = NA))
      grid.rect(x = 0, y = 0, width = 1, height = 0.58,
                just = c("left","bottom"), gp = gpar(fill = "white", col = NA))
      grid.text("Pesticide Risk",   x = 0.08, y = 0.80, just = c("left","centre"),
                gp = gpar(fontsize = 40, fontface = "bold", col = "white"))
      grid.text("Assessment Report", x = 0.08, y = 0.70, just = c("left","centre"),
                gp = gpar(fontsize = 40, fontface = "bold", col = "white"))
      grid.text("Environmental and human health risk screening",
                x = 0.08, y = 0.608, just = c("left","centre"),
                gp = gpar(fontsize = 13, col = "#b3e0ff"))
      grid.segments(x0 = 0.08, x1 = 0.92, y0 = 0.545, y1 = 0.545,
                    gp = gpar(col = "#e0e0e0", lwd = 0.8))
      grid.text("REPORT SUMMARY", x = 0.08, y = 0.530, just = c("left","centre"),
                gp = gpar(fontsize = 8, fontface = "bold", col = COL_BLUE,
                          letterSpacing = 1.5))
      info_items <- list(
        list(label = "Date",               value = RPT_DATE),
        list(label = "Chemicals assessed", value = as.character(nrow(risk_df))),
        list(label = "Runoff or irrigation",
             value = irrigation_label(scenario_settings()$runoff_source))
      )
      if (!is.null(coords))
        info_items <- c(info_items,
                        list(list(label = "Location",
                                  value = sprintf("%.4f\u00b0 lat,  %.4f\u00b0 lon",
                                                  coords$lat, coords$lon))))
      y0_info <- 0.490
      dy_info <- 0.058
      for (k in seq_along(info_items)) {
        yy <- y0_info - (k - 1) * dy_info
        grid.text(info_items[[k]]$label, x = 0.08, y = yy + 0.020,
                  just = c("left","centre"),
                  gp = gpar(fontsize = 8, col = "#888888"))
        grid.text(info_items[[k]]$value, x = 0.08, y = yy,
                  just = c("left","centre"),
                  gp = gpar(fontsize = 11, fontface = "bold", col = COL_DARK))
      }
      grid.rect(x = 0, y = 0, width = 1, height = 0.055,
                just = c("left","bottom"), gp = gpar(fill = "#f4f6f8", col = NA))
      grid.segments(0, 0.055, 1, 0.055, gp = gpar(col = "#dddddd", lwd = 0.5))
      grid.text(
        "For informational purposes only. Not a substitute for professional agronomic or environmental advice.",
        x = 0.5, y = 0.027, just = c("centre","centre"),
        gp = gpar(fontsize = 8, col = "#999999", fontface = "italic")
      )

      # ---- Page 2: site and soil ----
      draw_chrome("Site location and soil properties", 2)

      loc_df <- {
        rows <- c("Assessment date" = RPT_DATE)
        if (!is.null(coords))
          rows <- c(rows, "Latitude"  = sprintf("%.5f\u00b0", coords$lat),
                    "Longitude" = sprintf("%.5f\u00b0", coords$lon))
        data.frame(Parameter = names(rows), Value = unname(rows), stringsAsFactors = FALSE)
      }
      soil_df <- if (!is.null(sg)) {
        data.frame(
          Parameter = c("Soil pH (15-30 cm)", "Organic carbon (%)", "Organic matter (%)",
                        "Clay (%)", "Silt (%)", "Sand (%)",
                        "Bulk density (g/cm\u00b3)", "CEC (cmol/kg)",
                        "Texture class"),
          Value     = c(fmt(sg$ph,"%.1f"), fmt(sg$soc_pct,"%.2f"), fmt(sg$om_pct,"%.1f"),
                        fmt(sg$clay_pct,"%.1f"), fmt(sg$silt_pct,"%.1f"), fmt(sg$sand_pct,"%.1f"),
                        fmt(sg$bd_gcm3,"%.2f"), fmt(sg$cec_cmol,"%.1f"),
                        if (!is.null(sg$texture_class) && !is.na(sg$texture_class))
                          sg$texture_class else "N/A"),
          stringsAsFactors = FALSE)
      } else {
        data.frame(Parameter = "Soil data", Value = "Not fetched (click map to load)",
                   stringsAsFactors = FALSE)
      }

      th_blue  <- rpt_ttheme(COL_BLUE)
      th_green <- rpt_ttheme(COL_GREEN)
      th_dark  <- rpt_ttheme(COL_DARK)

      sec_label <- function(txt, col) {
        textGrob(txt, x = 0, hjust = 0,
                 gp = gpar(fontsize = 10, fontface = "bold", col = col))
      }

      grid.arrange(
        arrangeGrob(sec_label("Location", COL_DARK),
                    tableGrob(loc_df,  rows = NULL, theme = th_dark),
                    ncol = 1, heights = c(0.07, 0.93)),
        arrangeGrob(sec_label(paste0("Soil properties (", soil_source_name() %||% "SoilGrids 250 m", ")"), COL_BLUE),
                    tableGrob(soil_df, rows = NULL, theme = th_blue),
                    ncol = 1, heights = c(0.07, 0.93)),
        ncol = 2, widths = c(0.38, 0.62),
        newpage = FALSE
      )
      popViewport()

      # ---- Page 3: water balance and management ----
      draw_chrome("Water balance and management settings", 3)

      ov <- function(user, site, f) if (!is.null(user) && length(user) == 1 && !is.na(user)) paste(fmt(user, f), "(user)") else fmt(site, f)
      water_df <- if (!is.null(ha)) {
        wdf <- data.frame(
          Parameter = c("Mean annual rainfall (mm)", "Annual runoff (mm)",
                        "Evapotranspiration (mm)", "Mean temperature (\u00b0C)",
                        "Mean slope (\u00b0)", "Runoff fraction", "Infiltration and ET fraction"),
          Value     = c(fmt(ha$pre_mm,"%.0f"), fmt(ha$run_mm,"%.0f"), fmt(ha$aet_mm,"%.0f"),
                        ov(user_temp_c(), ha$temp_c, "%.1f"), ov(user_slope_dg(), ha$slp_dg, "%.1f"),
                        ov(user_runoff_frac(), ha$runoff_frac, "%.2f"), fmt(ha$infil_frac,"%.2f")),
          stringsAsFactors = FALSE)
        rd_pdf <- rusle_data()
        if (!is.null(rd_pdf)) {
          wdf <- rbind(wdf, data.frame(
            Parameter = c("RUSLE R / LS / C", "NDVI window (C factor)"),
            Value     = c(sprintf("%s / %s / %s", fmt(ha$R_factor, "%.0f"), fmt(rd_pdf$LS_factor, "%.3g"),
                                  fmt(rd_pdf$C_factor, "%.3f")),
                          rd_pdf$ndvi_window %||% "N/A"),
            stringsAsFactors = FALSE))
        }
        wdf
      } else {
        data.frame(Parameter = "Water balance", Value = "Not fetched (click map to load)",
                   stringsAsFactors = FALSE)
      }
      sc <- scenario_settings(); bs <- baseline_settings()
      mgmt_df <- data.frame(
        Setting = c("Runoff or irrigation type", "Irrigation assumptions",
                    "Vegetative filter strip", "Conservation practices", "Practice benefit applied",
                    "Application rate", "Application date", "Application method", "Runoff timing",
                    "Airborne deposition share", "Crop and pests (Step 4)", "Products and rates (Step 5)",
                    "Baseline (Step 6)"),
        Value   = c(
          irrigation_label(sc$runoff_source),
          {
            iv <- irrigation_values(sc$runoff_source, sc$irrig_source %||% "central")
            if (!is.null(iv))
              paste(strwrap(sprintf("Sediment %.2f g/L; runoff fraction %.2f; source: %s", iv$SC, iv$runoff_fraction,
                                    sub(":\\s*[0-9.]+$", "", iv$source)), 38), collapse = "\n")
            else "RUSLE sediment concentration; curve-number runoff"
          },
          paste0(sc$Vegetative_Filter_Strips, " m"),
          {
            all_p <- c(if (sc$tillage_type != "none") gsub("_", " ", sc$tillage_type),
                       if (length(sc$land_practices) > 0) gsub("_", " ", sc$land_practices))
            if (length(all_p) > 0) paste(all_p, collapse = ", ") else "None selected"
          },
          paste0(benefit_share(), "% of reported reduction"),
          paste(strwrap(paste0(100 - sc$rate_reduction, "% of the rate shown in the risk table"), 38),
                collapse = "\n"),
          paste(strwrap(paste(c(format(sc$app_date, "%d %b %Y"),
                  if (!is.null(irr_schedule(sc))) with(irr_schedule(sc), sprintf(
                    "irrigation every %g days from day %g, %g mm each", interval, first_offset, depth_mm))),
                collapse = "; "), 38), collapse = "\n"),
          paste(strwrap(paste0(names(DRIFT_METHODS)[DRIFT_METHODS == sc$app_method],
                               if (!grepl("^Vertical", sc$app_method))
                                 paste0(", ", unname(SPRAY_QUALITY_SHORT[sc$spray_quality]), " spray"),
                               "; flowering vegetation ", sc$habitat_dist, " m downwind"), 38), collapse = "\n"),
          paste(strwrap(timing_text(), 38), collapse = "\n"),
          paste0(num_or(input$vapour_dep_pct, VAPOUR_DEP_DEFAULT_PCT), "% (assumed)"),
          paste(strwrap(if (length(target_hosts()) + length(target_pests()) == 0) "Not refined" else
                  paste0("Crop or situation: ", if (length(target_hosts())) paste(title_case(host_label(target_hosts())), collapse = ", ") else "any",
                         "; pests: ", if (length(target_pests())) paste(title_case(pest_label(target_pests())), collapse = ", ") else "any"),
                38), collapse = "\n"),
          {
            pl <- product_plan(); u <- if (!is.null(pl)) pl[!is.na(pl$kg_ha), , drop = FALSE] else NULL
            ri <- rate_info(); sel <- selected_chemicals(); if (length(sel)) ri <- ri[ri$active %in% sel, , drop = FALSE]
            ns <- table(factor(ri$source[ri$source != "user"], names(RATE_SOURCE)))
            ns <- ns[ns > 0]
            paste(strwrap(paste(c(
              if (!is.null(u) && nrow(u)) paste(ifelse(is.na(u$product), sprintf("%s: %s kg ai/ha", u$active, signif(u$kg_ha, 3)),
                                        sprintf("%s: %s %s of %s (%s kg ai/ha)", u$active, signif(u$product_rate, 3), u$unit,
                                                u$product, signif(u$kg_ha, 3))), collapse = "; "),
              if (length(ns)) paste(sprintf("%s: %d %s", RATE_SOURCE[names(ns)], as.integer(ns),
                                            ifelse(as.integer(ns) == 1, "chemical", "chemicals")), collapse = "; ")),
              collapse = "; "), 38), collapse = "\n")
          },
          paste(strwrap(settings_label(bs), 38), collapse = "\n")
        ), stringsAsFactors = FALSE
      )

      top_table <- function(d, th) {
        tg <- tableGrob(d, rows = NULL, theme = th)
        tg$vp <- viewport(y = unit(1, "npc") - 0.5 * sum(tg$heights)); tg
      }
      grid.arrange(
        arrangeGrob(sec_label("Water balance (HydroATLAS)", COL_BLUE),
                    top_table(water_df, th_blue),
                    ncol = 1, heights = c(0.06, 0.94)),
        arrangeGrob(sec_label("Management settings", COL_GREEN),
                    top_table(mgmt_df, th_green),
                    ncol = 1, heights = c(0.06, 0.94)),
        ncol = 2, widths = c(0.50, 0.50),
        newpage = FALSE
      )
      popViewport()

      # ---- Risk summary pages ----
      key_df <- as.data.frame(setNames(as.list(rep("", 5)), RISK_LEVELS), stringsAsFactors = FALSE)
      key_grob <- tableGrob(key_df, rows = NULL,
                            theme = gridExtra::ttheme_default(
                              colhead = list(
                                bg_params = list(fill = unname(RISK_PAL), col = NA),
                                fg_params = list(col = c("#333333","#333333","#333333","white","white"),
                                                 fontsize = 9, fontface = "bold")),
                              core = list(bg_params = list(fill = "white", col = NA),
                                          fg_params = list(fontsize = 1)),
                              padding = unit(c(3, 2), "mm")))
      risk_pdf <- risk_df
      risk_pdf$`Rate source` <- sub("Reference (no label rate for this use)", "Reference (no label rate)",
                                    sub("Label maximum, all uses", "Label max, all uses", risk_pdf$`Rate source`, fixed = TRUE),
                                    fixed = TRUE)
      fit_table <- function(df, avail_mm = 185) {
        fs <- 8
        repeat {
          tg <- tableGrob(df, rows = NULL, theme = rpt_ttheme(COL_DARK, fs = fs, pad = c(2.5, 2.2)))
          w  <- sum(convertWidth(tg$widths, "mm", valueOnly = TRUE))
          if (w <= avail_mm || fs <= 5.5) break
          fs <- max(5.5, fs * avail_mm / w * 0.98)
        }
        tg$vp <- viewport(y = unit(1, "npc") - 0.5 * sum(tg$heights))
        tg
      }
      for (pg in seq_len(n_tab_pg)) {
        draw_chrome(if (n_tab_pg > 1) sprintf("Risk classification summary (%d of %d)", pg, n_tab_pg) else
                      "Risk classification summary", 3L + pg)
        part <- if (nrow(risk_pdf) == 0) data.frame(Pesticide = "No chemicals to show") else
          risk_pdf[((pg - 1L) * ROWS_PP + 1L):min(nrow(risk_pdf), pg * ROWS_PP), , drop = FALSE]
        grid.arrange(
          colour_risk_grob(fit_table(part), part),
          arrangeGrob(textGrob("Risk level key:", x = 0, hjust = 0,
                               gp = gpar(fontsize = 9, fontface = "bold", col = COL_DARK)),
                      key_grob, ncol = 2, widths = c(0.15, 0.85)),
          ncol = 1, heights = c(0.92, 0.08), newpage = FALSE)
        popViewport()
      }

      # ---- Risk plot pages ----
      plot_list <- list(
        list(p = p_aquatic,     title = "Aquatic risk"),
        list(p = p_terrestrial, title = "Pollinator risk"),
        list(p = p_below,       title = "Soil fauna risk"),
        list(p = p_human,       title = "Human health hazard profile")
      )
      for (k in seq_along(plot_list)) {
        draw_chrome(plot_list[[k]]$title, 3L + n_tab_pg + k)
        grid.draw(ggplotGrob(plot_list[[k]]$p))
        popViewport()
      }

      # ---- Disclaimer and references ----
      draw_chrome("Disclaimer and references", TOTAL_PG)
      disc <- c("For information only; not a substitute for professional agronomic or environmental advice.",
                "Results are relative comparisons under simplified assumptions, not safety determinations.",
                "Tank mixes, metabolites and indirect ecological effects are not modelled.",
                "Registered uses and rates are read automatically from Australian product labels (APVMA). They apply in Australia only, may contain errors and may have changed since.",
                "Whatever this report shows, always follow the current product label and your local regulations.",
                "No warranty is given. The authors accept no responsibility for errors in the label information or for decisions made with this tool.")
      st <- scenario_settings()
      iv_key <- {
        iv <- irrigation_values(st$runoff_source, st$irrig_source %||% "central")
        k <- c(irmak = "Irmak", roth = "Roth", rohde = "Rohde", reddy = "Reddy", long = "Long", araujo = "Araujo",
               radmanesh = "Radmanesh", evans_pivot = "Evans", wqsah_burdekin = "Water Quality Science", crdc = "CRDC", lentz = "Lentz")
        if (is.null(iv)) character(0) else names(k)[vapply(k, function(x) grepl(x, iv$source %||% "", fixed = TRUE), TRUE)]
      }
      keys <- c("apvma_pubcris", "apvma_labels", "ppdb", "anzg", "silburn2024", "linke", "funk", "munoz", "grundy", "poggio", "hengl",
                "nrcs2004", "nrcs2009", "apvma_sdram", "apvma_sdrat", "usepa", "gentil", "houbraken", "asabe", "spencer", "ebbert",
                "zhang", "klein", "fleming", iv_key)
      refs <- sort(unlist(REFS[intersect(unique(keys), names(REFS))]))
      wrap <- function(x, w) paste(strwrap(x, w), collapse = "\n")
      body <- c("Disclaimer", paste0("- ", vapply(disc, wrap, "", w = 125)), "", "References", vapply(refs, wrap, "", w = 135))
      n_lines <- vapply(strsplit(body, "\n"), length, 1L)
      fs <- 6.5; lh <- fs * 1.25 / 72 / (0.85 * 11.69)
      y <- 0.99
      for (i in seq_along(body)) {
        bold <- body[i] %in% c("Disclaimer", "References")
        grid.text(body[i], x = 0.06, y = y, just = c("left", "top"),
                  gp = gpar(fontsize = if (bold) 9 else fs, fontface = if (bold) "bold" else "plain", col = COL_DARK))
        y <- y - lh * (n_lines[i] + 0.35) * (if (bold) 1.4 else 1)
      }
      popViewport()

      dev.off(dv)
    }
  )

  output$export_settings <- downloadHandler(
    filename = function() paste0("pesticide_assessment_", Sys.Date(), ".json"),
    content = function(file) {
      sg <- soilgrids_data()
      settings <- list(
        exported           = format(Sys.time(), "%Y-%m-%d %H:%M"),
        format_version     = 2,
        app_version        = APP_VERSION,
        rate_tables_built  = if (!is.null(use_rates)) as.character(use_rates$built %||% NA) else NULL,
        field_outline      = field_geom(),
        map_lat            = if (!is.null(map_coords())) map_coords()$lat else NULL,
        map_lon            = if (!is.null(map_coords())) map_coords()$lon else NULL,
        soil_source        = soil_source_name(),
        soil_ph            = if (!is.null(sg)) sg$ph else NULL,
        soil_soc_pct       = if (!is.null(sg)) sg$soc_pct else NULL,
        soil_clay_pct      = if (!is.null(sg)) sg$clay_pct else NULL,
        soil_bd_gcm3       = if (!is.null(sg)) sg$bd_gcm3 else NULL,
        baseline_settings  = lapply(baseline_settings(), function(v) if (inherits(v, "Date")) as.character(v) else v),
        scenario_settings  = lapply(scenario_settings(), function(v) if (inherits(v, "Date")) as.character(v) else v),
        practice_benefit_pct = benefit_share(),
        airborne_deposition_pct = num_or(input$vapour_dep_pct, VAPOUR_DEP_DEFAULT_PCT),
        runoff_event_mm    = RUNOFF_EVENT_MM,
        runoff_timing      = if (!is.null(scen_timing())) with(scen_timing(), list(
          rainfall_source = source, curve_number = CN, soil_group = hsg, n_years = n_years,
          median_days = if (is.finite(median_days)) median_days else NULL,
          p_within_7d = p7, p_within_14d = p14, p_within_30d = p30)) else NULL,
        pest_groups        = as.list(selected_pest_groups()),
        target_crops       = as.list(target_hosts()),
        target_pests       = as.list(target_pests()),
        target_crop_names  = as.list(unname(host_label(target_hosts()))),
        target_pest_names  = as.list(unname(pest_label(target_pests()))),
        product_plan       = product_plan(),
        chemical_rates     = { ri <- rate_info(); sel <- selected_chemicals(); if (length(sel)) ri <- ri[ri$active %in% sel, , drop = FALSE]
                               data.frame(active = ri$active, kg_ai_ha = signif(ri$kg, 4),
                                          source = unname(RATE_SOURCE[ri$source])) },
        selected_chemicals = as.list(unname(selected_chemicals())),
        user_overrides     = list(runoff_fraction = user_runoff_frac(), slope_deg = user_slope_dg(),
                                  temp_c = user_temp_c(), clay_pct = user_clay_pct(),
                                  sand_pct = user_sand_pct(), silt_pct = user_silt_pct(),
                                  bd_gcm3 = user_bd_gcm3(), soc_pct = user_soc_pct(),
                                  curve_number = user_cn()),
        site = list(hydroatlas = hydroatlas_data(), soil = sg, rusle = rusle_data(), soil_ph = sel_ph(), soil_om = sel_om(),
                    rain = { r <- rain_daily_data()
                             if (is.null(r)) NULL else list(source = attr(r, "source"), date = as.character(r$date),
                                                            P_mm = round(r$P_mm, 3)) })
      )
      write(jsonlite::toJSON(settings, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA), file)
    }
  )

  observeEvent(input$clear_settings, {
    showModal(modalDialog(
      title = "Clear all and restart?",
      "This clears the field, its data and all choices.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_clear", "Clear all", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_clear, {
    removeModal()
    map_coords(NULL)
    hydroatlas_data(NULL)
    rain_daily_data(NULL)
    rusle_data(NULL)
    soilgrids_data(NULL); soil_source_name(NULL); soil_source_id(NULL)
    soilgrids_status_msg(list(type = "none", msg = ""))
    hydroatlas_status_msg(list(type = "none", msg = ""))
    sel_om("moderate")
    sel_ph(7.0)
    field_geom(NULL)
    user_runoff_frac(NULL); user_slope_dg(NULL); user_temp_c(NULL); user_cn(NULL)
    user_clay_pct(NULL);    user_sand_pct(NULL); user_silt_pct(NULL)
    user_bd_gcm3(NULL);     user_soc_pct(NULL)
    selected_pest_groups(character(0))
    clear_targets()
    baseline_settings(default_settings())
    scenario_init(NULL)
    adv_init(list(practice_benefit_pct = PRACTICE_BENEFIT_DEFAULT, vapour_dep_pct = VAPOUR_DEP_DEFAULT_PCT))
    updateSliderInput(session, "practice_discount_pct", value = 100 - PRACTICE_BENEFIT_DEFAULT)
    updateNumericInput(session, "vapour_dep_pct", value = VAPOUR_DEP_DEFAULT_PCT)
    current_page("landing")
    showNotification("All settings cleared.", type = "message", duration = 3)
  })
}

shinyApp(ui = ui, server = server)
