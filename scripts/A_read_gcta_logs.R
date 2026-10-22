# =============================================================================
# Script: A_read_gcta_logs.R
#
# Purpose:
#   Extract variance component estimates and heritability estimates from GCTA
#   REML log files, and calculate approximate 95% confidence intervals.
#
# Input:
#   - GCTA REML .log files in the working directory
#
# Outputs:
#   - gcta_reml_results.rds: A combined tibble containing sample sizes, 
#     variance component estimates, standard errors, and 95% confidence 
#     intervals for each log file.
#
# Methods:
#   - Parses sex, statin status, and LDL-C phenotype type from filenames.
#   - Extracts the number of observations and REML summary estimates from
#     GCTA log files.
#   - Calculates approximate 95% confidence intervals as estimate +/- 1.96 SE.
#   - Combines results from all log files into a single tibble.
#
# Notes:
#   - The script assumes that GCTA log files contain a summary result block
#     with estimates for V(G), V(e), Vp, and V(G)/Vp.
#   - Confidence intervals use a normal approximation and are not constrained
#     to biologically plausible parameter ranges.
#   - Run this script from the directory containing the GCTA log files.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

library(tidyverse)

# =============================================================================
# 1. READ GCTA REML LOG FILES
# =============================================================================

# Parse a single GCTA REML log file and extract variance component estimates,
# standard errors, sample size, and approximate 95% confidence intervals.

read_gcta_log <- function(file) {
  
  nm <- tools::file_path_sans_ext(basename(file))
  pieces <- str_split(nm, "_")[[1]]
  
  sex <- pieces[1]
  
  if (pieces[2] == "max") {
    # Format: sex_max_nostat_max_unadj_pheno
    statin_status <- pieces[3]
    type <- paste(pieces[4:5], collapse = "_")
  } else {
    # Format: sex_nostat_unadj_pheno
    statin_status <- pieces[2]
    type <- paste(pieces[3], collapse = "_")
  }
  
  txt <- readLines(file)

  obs_line <- grep("observations,", txt, value = TRUE)
  
  n <- str_match(obs_line, "^(\\d+) observations")[,2] |>
    as.numeric()
  
  start <- grep("^Summary result of REML analysis:", txt)
  end <- grep("^Sampling variance/covariance", txt)
  
  summary_block <- txt[(start + 2):(end - 2)]

  summary_block <- summary_block[nzchar(summary_block)]

  tmp <- str_split(summary_block, "\\s+", simplify = TRUE)
  

  df <- data.frame(
    source = tmp[,1],
    estimate = as.numeric(tmp[,2]),
    se = as.numeric(tmp[,3])
  )
  
# =============================================================================
# 2. CALCULATE 95% CONFIDENCE INTERVALS
# =============================================================================
  
# Calculate approximate two-sided 95% confidence intervals using a
# normal approximation to the sampling distribution of each estimate.
  
  df <- df |>
    mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se
    )
  
# =============================================================================
# 3. FORMAT OUTPUT
# =============================================================================
  
# Reshape the parsed estimates into one row per log file, with separate
# columns for each variance component, its standard error, and its
# confidence interval.
  
  out <- tibble(
    sex = sex,
    statin_status = statin_status,
    ldl_type = type,
    n = n,
    
    # Additive genetic variance
    VG = df$estimate[df$source == "V(G)"],
    VG_SE = df$se[df$source == "V(G)"],
    VG_lower = df$lower[df$source == "V(G)"],
    VG_upper = df$upper[df$source == "V(G)"],
    
    # Residual variance
    Ve = df$estimate[df$source == "V(e)"],
    Ve_SE = df$se[df$source == "V(e)"],
    Ve_lower = df$lower[df$source == "V(e)"],
    Ve_upper = df$upper[df$source == "V(e)"],
    
    # Total phenotypic variance
    Vp = df$estimate[df$source == "Vp"],
    Vp_SE = df$se[df$source == "Vp"],
    Vp_lower = df$lower[df$source == "Vp"],
    Vp_upper = df$upper[df$source == "Vp"],
    
    # SNP-based heritability (V(G) / Vp)
    h2 = df$estimate[df$source == "V(G)/Vp"],
    h2_SE = df$se[df$source == "V(G)/Vp"],
    h2_lower = df$lower[df$source == "V(G)/Vp"],
    h2_upper = df$upper[df$source == "V(G)/Vp"]
  )
  
  out
}

# =============================================================================
# 4. PROCESS ALL LOG FILES
# =============================================================================

# Identify all GCTA REML log files in the working directory.

files <- list.files(
  pattern = "\\.log$",
  full.names = TRUE
)

# Parse each log file and combine the results into a single tibble.
results <- map_dfr(files, read_gcta_log)

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================

# Save the combined REML results for use in downstream analyses.
saveRDS(results, "gcta_reml_results.rds")
