# =============================================================================
# LDL-C Sex-Stratified Z-Score Analysis
#
# Purpose:
# Compare female and male GWAS effect estimates at selected loci using a
# z-score difference test, with Bonferroni correction across tested loci.
#
# Method:
# Select lead variants using positional clumping of combined-GWAS results,
# with additional sex-specific candidates included by the selection procedure.
# Gene annotation uses the Ensembl REST API; rsID lookup uses NCBI dbSNP.
#
# Inputs:
#   GWAS_AoU_ADJLDL_FEMALE.rds
#   GWAS_AoU_ADJLDL_MALE.rds
#   GWAS_AoU_ADJLDL_COMBINED.rds
#
# Output:
#   lead_snp_summary_table_final.csv
#
# Notes:
#   Input variant IDs must use the format chr:position:ref:alt.
#   Genomic positions used for annotation must correspond to GRCh38.
#   Positional clumping uses physical distance and does not evaluate LD.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

# Install required packages if they are not already available.
required_packages <- c("tidyverse", "httr", "jsonlite")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(httr)
library(jsonlite)


# =============================================================================
# 1. LOAD DATA
# =============================================================================

female <- readRDS('GWAS_AoU_ADJLDL_FEMALE.rds')
male <- readRDS('GWAS_AoU_ADJLDL_MALE.rds')
combined <- readRDS('GWAS_AoU_ADJLDL_COMBINED.rds')

# =============================================================================
# 2. COMPUTE Z-SCORE DIFFERENCE TEST (all SNPs)
# =============================================================================

merged <- merge(female, male, by="ID", suffixes=c("_fem","_male"))
merged$Z_diff <- (merged$BETA_fem - merged$BETA_male) /
  sqrt(merged$SE_fem^2 + merged$SE_male^2)
merged$P_diff <- 2 * pnorm(-abs(merged$Z_diff))

# Apply a Bonferroni correction across all variants tested in the combined
# GWAS. This threshold is also used to identify significant variants in
# the sex-stratified GWAS results.
gws <- 0.05 / nrow(combined)
cat("GWS threshold:", gws, "\n")

# Convert -log10(p) values to p-values for significance filtering.
female$P <- 10^(-female$LOG10P)
male$P <- 10^(-male$LOG10P)
combined$P <- 10^(-combined$LOG10P)

# =============================================================================
# 3. CLUMPING: identify lead SNPs
# =============================================================================

# Parse chromosome and position from ID (format: chr19:44892652:C:G)
combined_gws <- combined |>
  dplyr::filter(P < gws) |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE) |>
  mutate(chr_num = as.numeric(str_remove(chr, "chr"))) |>
  arrange(P)

# Greedy positional clumping: process variants from smallest to largest
# p-value, retain the most significant remaining variant, and remove
# variants on the same chromosome within the specified physical window.
# This uses genomic distance only; linkage disequilibrium is not evaluated.
clump_greedy <- function(df, window_kb=1000) {
  lead_snps <- list()
  while(nrow(df) > 0) {
    lead <- df[1, ]
    lead_snps[[length(lead_snps) + 1]] <- lead
    df <- df |>
      dplyr::filter(!(chr == lead$chr & abs(pos - lead$pos) < window_kb * 1000))
  }
  bind_rows(lead_snps)
}

lead_snps_combined <- clump_greedy(combined_gws)
cat("Number of independent loci from combined GWAS:", nrow(lead_snps_combined), "\n")

# Find any sex-specific GWS SNPs not captured by combined leads
female_gws <- female |>
  dplyr::filter(P < gws) |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE)

male_gws <- male |>
  dplyr::filter(P < gws) |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE)

sex_specific_candidates <- bind_rows(female_gws, male_gws) |>
  distinct(ID, .keep_all=TRUE) |>
  dplyr::filter(!ID %in% lead_snps_combined$ID) |>
  dplyr::filter(!map_lgl(seq_len(n()), function(i) {
    any(chr[i] == lead_snps_combined$chr &
          abs(pos[i] - lead_snps_combined$pos) < 1e6)
  }))

cat("Additional sex-specific loci not in combined:", nrow(sex_specific_candidates), "\n")

lead_snp_ids <- unique(c(lead_snps_combined$ID, sex_specific_candidates$ID))
cat("Total independent loci for z-score testing:", length(lead_snp_ids), "\n")

# All GWS SNP IDs across any analysis 
gws_snps <- unique(c(
  female$ID[female$P < gws],
  male$ID[male$P < gws],
  combined$ID[combined$P < gws]
))
cat("Total GWS SNPs across all analyses:", length(gws_snps), "\n")

# =============================================================================
# 4. APPLY Z-SCORE TEST TO LEAD SNPs ONLY
# =============================================================================

# Apply Bonferroni correction across the selected variants available in
# both sex-stratified GWAS results.
merged_gws <- merged[merged$ID %in% lead_snp_ids, ]

# Check how many selected variants are available in both sex-stratified GWAS.
cat("Selected variants available in both GWAS:", nrow(merged_gws), "\n")

merged_gws$P_diff_adj  <- p.adjust(merged_gws$P_diff, method="bonferroni")
merged_gws$sex_diff_sig <- merged_gws$P_diff_adj < 0.05

cat("Z-diff results:\n")
print(table(merged_gws$sex_diff_sig))

# Inspect significant SNPs
sig_snps <- merged_gws[merged_gws$sex_diff_sig == TRUE, ]
sig_snps$direction <- case_when(
  10^(-sig_snps$LOG10P_fem) < gws & 10^(-sig_snps$LOG10P_male) >= gws ~ "female-specific",
  10^(-sig_snps$LOG10P_male) < gws & 10^(-sig_snps$LOG10P_fem) >= gws ~ "male-specific",
  sign(sig_snps$BETA_fem) != sign(sig_snps$BETA_male) ~ "opposite direction",
  abs(sig_snps$BETA_fem) > abs(sig_snps$BETA_male) ~ "larger effect in females",
  abs(sig_snps$BETA_male) > abs(sig_snps$BETA_fem) ~ "larger effect in males",
  TRUE ~ "check manually"
)

cat("Significant sex-differentiated SNPs:\n")
print(sig_snps |> dplyr::select(ID, BETA_fem, BETA_male, Z_diff, P_diff_adj, direction))

# =============================================================================
# 5. BUILD SUMMARY TABLE (base, before annotation)
# =============================================================================

combined_sub <- combined[combined$ID %in% lead_snp_ids,
                         c("ID","BETA","SE","LOG10P")]
names(combined_sub)[2:4] <- c("BETA_combined","SE_combined","LOG10P_combined")

summary_table <- merged_gws |>
  dplyr::select(ID, BETA_fem, SE_fem, LOG10P_fem,
                BETA_male, SE_male, LOG10P_male,
                Z_diff, P_diff, P_diff_adj, sex_diff_sig) |>
  left_join(combined_sub, by="ID")

# =============================================================================
# 6. GENE ANNOTATION 
# =============================================================================

# This endpoint uses the current Ensembl assembly (GRCh38).
ensembl_rest <- "https://rest.ensembl.org"

# Annotate each variant by querying the Ensembl REST API for protein-coding
# genes overlapping a genomic interval centered on the variant position.
# The genomic coordinates are assumed to use GRCh38.

# Query protein-coding genes overlapping a genomic window.
query_gene_window <- function(chr, pos, window, snp_id) {
  
  start <- max(1, as.integer(pos - window))
  end   <- as.integer(pos + window)
  
  url <- paste0(
    ensembl_rest,
    "/overlap/region/human/",
    chr, ":", start, "-", end
  )
  
  res <- tryCatch({
    
    response <- httr::RETRY(
      "GET",
      url,
      query = list(
        feature = "gene",
        biotype = "protein_coding"
      ),
      httr::accept_json(),
      times = 5,
      pause_base = 2,
      pause_cap = 30,
      terminate_on = c(400, 401, 403, 404),
      quiet = TRUE
    )
    
    status <- httr::status_code(response)
    
    if (status >= 400) {
      message(
        "Query failed for ", snp_id,
        " (HTTP ", status, ")"
      )
      return(NULL)
    }
    
    response_text <- httr::content(
      response,
      as = "text",
      encoding = "UTF-8"
    )
    
    jsonlite::fromJSON(
      response_text,
      simplifyDataFrame = TRUE
    )
    
  }, error = function(e) {
    message(
      "Query failed for ", snp_id, ": ",
      conditionMessage(e)
    )
    NULL
  })
  
  res <- as_tibble(res)
  
  # Handle empty responses
  if (nrow(res) == 0) {
    return(tibble(
      chromosome_name = character(),
      start_position = integer(),
      end_position = integer(),
      external_gene_name = character(),
      gene_biotype = character(),
      SNP_ID = character(),
      SNP_pos = integer(),
      distance = integer()
    ))
  }
  
  # Check that required Ensembl fields exist
  required_cols <- c(
    "seq_region_name", "start", "end",
    "external_name", "biotype"
  )
  
  missing_cols <- setdiff(required_cols, names(res))
  
  if (length(missing_cols) > 0) {
    message(
      "Unexpected Ensembl response for ", snp_id,
      ". Missing columns: ",
      paste(missing_cols, collapse = ", "),
      ". Available columns: ",
      paste(names(res), collapse = ", ")
    )
    
    return(tibble(
      chromosome_name = character(),
      start_position = integer(),
      end_position = integer(),
      external_gene_name = character(),
      gene_biotype = character(),
      SNP_ID = character(),
      SNP_pos = integer(),
      distance = integer()
    ))
  }
  
  res |>
    transmute(
      chromosome_name = as.character(seq_region_name),
      start_position = as.integer(start),
      end_position = as.integer(end),
      external_gene_name = as.character(external_name),
      gene_biotype = as.character(biotype),
      SNP_ID = snp_id,
      SNP_pos = as.integer(pos)
    ) |>
    filter(gene_biotype == "protein_coding") |>
    mutate(
      distance = if_else(
        SNP_pos >= start_position & SNP_pos <= end_position,
        0L,
        pmin(
          abs(SNP_pos - start_position),
          abs(SNP_pos - end_position)
        )
      )
    )
}

# Select the closest protein-coding gene for each SNP.
select_nearest_gene <- function(gene_data) {
  gene_data |>
    group_by(SNP_ID) |>
    slice_min(
      order_by = distance,
      n = 1,
      with_ties = FALSE
    ) |>
    ungroup() |>
    dplyr::select(
      SNP_ID,
      Gene = external_gene_name,
      distance
    )
}


positions <- summary_table |>
  separate(
    ID,
    into = c("chr", "pos", "ref", "alt"),
    sep = ":",
    remove = FALSE,
    convert = TRUE
  ) |>
  mutate(
    chr_num = as.numeric(str_remove(chr, "chr")),
    pos = as.integer(pos)
  ) |>
  filter(
    !is.na(chr_num),
    chr_num %in% 1:22,
    !is.na(pos)
  ) |>
  arrange(chr_num, pos)

# =============================================================================
# NARROW WINDOW: 10 KB
# =============================================================================

# First search within 10 kb of each selected variant to identify nearby
# protein-coding genes.
gene_annotations_per_snp <- purrr::map_dfr(
  seq_len(nrow(positions)),
  function(i) {
    query_gene_window(
      chr = positions$chr_num[i],
      pos = positions$pos[i],
      window = 10000,
      snp_id = positions$ID[i]
    )
  }
)

annotations_final <- select_nearest_gene(
  gene_annotations_per_snp
)

# =============================================================================
# WIDE WINDOW: 500 KB FOR SPECIFIED SNPS
# =============================================================================

# Re-query these prespecified variants using a wider 500 kb window.
# The wider search is an annotation exception and does not change variant
# selection or the statistical tests.
wide_snps <- positions |>
  filter(
    ID %in% c(
      "chr2:21090558:G:T",
      "chr8:125467846:G:A",
      "chr20:17864040:G:T",
      "chr19:11080521:G:A"
    )
  )

gene_annotations_wide_raw <- purrr::map_dfr(
  seq_len(nrow(wide_snps)),
  function(i) {
    query_gene_window(
      chr = wide_snps$chr_num[i],
      pos = wide_snps$pos[i],
      window = 500000,
      snp_id = wide_snps$ID[i]
    )
  }
)

gene_annotations_wide <- select_nearest_gene(
  gene_annotations_wide_raw
)


annotations_final_patched <- annotations_final |>
  filter(!SNP_ID %in% wide_snps$ID) |>
  bind_rows(gene_annotations_wide) |>
  group_by(SNP_ID) |>
  slice_min(
    order_by = distance,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    # Preserve existing manual override.
    Gene = if_else(
      SNP_ID == "chr19:11080521:G:A",
      "LDLR",
      Gene
    )
  )

cat("Final annotations:", nrow(annotations_final_patched), "\n")
print(annotations_final_patched, n = Inf)

# =============================================================================
# 7. RSID LOOKUP via NCBI dbSNP
# =============================================================================

# Query NCBI dbSNP for an rsID using chromosome and GRCh38 position.
# Variants without a returned rsID retain their chromosome-position-allele ID.
query_dbsnp <- function(chrom, pos) {
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=snp",
    "&term=", chrom, "[CHR]+AND+", pos, "[CHRPOS38]",
    "&retmode=json&retmax=5"
  )
  res <- tryCatch(GET(url, timeout(10)), error=function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(NULL)
  
  content <- fromJSON(rawToChar(res$content))
  ids <- content$esearchresult$idlist
  
  if (length(ids) == 0) return(data.frame(pos=pos, rsID=NA_character_))
  data.frame(pos=pos, rsID=paste0("rs", ids[1]))
}

rsid_results <- map_dfr(seq_len(nrow(positions)), function(i) {
  cat("Looking up rsID for", positions$ID[i], "\n")
  res <- query_dbsnp(positions$chr_num[i], positions$pos[i])
  if (!is.null(res)) res$SNP_ID <- positions$ID[i]
  Sys.sleep(0.5)  # be polite to NCBI API
  res
})

cat("rsID lookup results:\n")
print(rsid_results)

# =============================================================================
# 8. BUILD FINAL SUMMARY TABLE
# =============================================================================

# Combine GWAS statistics, gene annotations, and rsIDs into a presentation-
# ready summary table. Effect estimates and standard errors are rounded for
# display; p-values are formatted in scientific notation.
summary_table_final <- summary_table |>
  left_join(annotations_final_patched, by=c("ID"="SNP_ID")) |>
  left_join(rsid_results |> dplyr::select(SNP_ID, rsID), by=c("ID"="SNP_ID")) |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE) |>
  mutate(chr_num = as.numeric(str_remove(chr, "chr"))) |>
  arrange(chr_num, pos) |>
  mutate(
    # Use rsID where available, fall back to positional ID
    Variant = ifelse(!is.na(rsID), rsID, ID),
    P_fem = 10^(-LOG10P_fem),
    P_male = 10^(-LOG10P_male),
    P_combined = 10^(-LOG10P_combined),
    across(c(BETA_fem, BETA_male, BETA_combined), ~round(.x, 3)),
    across(c(SE_fem, SE_male, SE_combined), ~round(.x, 3)),
    across(c(P_fem, P_male, P_combined), ~formatC(.x, format="e", digits=2)),
    `Z-diff` = round(Z_diff, 3),
    P_diff_adj = formatC(P_diff_adj, format="e", digits=2),
    sex_diff_sig = ifelse(sex_diff_sig, "Yes", "No"),
    # Fix floating point underflow
    P_combined = ifelse(P_combined == "0.00e+00", "<1e-300", P_combined),
    P_fem = ifelse(P_fem == "0.00e+00", "<1e-300", P_fem),
    P_male = ifelse(P_male == "0.00e+00", "<1e-300", P_male)
  ) |>
  dplyr::select(
    Chr = chr_num,
    Position = pos,
    Variant,
    Gene,
    `Beta (Combined)` = BETA_combined,
    `SE (Combined)` = SE_combined,
    `P (Combined)` = P_combined,
    `Beta (Female)` = BETA_fem,
    `SE (Female)` = SE_fem,
    `P (Female)` = P_fem,
    `Beta (Male)` = BETA_male,
    `SE (Male)` = SE_male,
    `P (Male)` = P_male,
    `Z-diff`,
    `P-diff (adj)` = P_diff_adj,
    `Sex-differentiated` = sex_diff_sig
  )

# Write the final annotated lead-variant table to the working directory.
write_csv(summary_table_final, 'lead_snp_summary_table_final.csv')
cat("Final summary table:\n")
print(summary_table_final)

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

# Save intermediate results required by downstream analysis scripts.
saveRDS(
  list(
    female = female,
    male = male,
    gws = gws,
    lead_snp_ids = lead_snp_ids,
    sig_snps = sig_snps,
    merged_gws = merged_gws,
    annotations_final_patched = annotations_final_patched,
    rsid_results = rsid_results
  ),
  file = "ldl_sex_difference_results.rds"
)

