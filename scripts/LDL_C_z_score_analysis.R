# =============================================================================
# LDL-C Sex-Stratified Z-Score Analysis
# Tests for sex-differentiated effect sizes at independent GWS loci from
# female, male, and combined GWAS using a z-score difference test
# Input: GWAS_AoU_ADJLDL_FEMALE.rds, GWAS_AoU_ADJLDL_MALE.rds,
#        GWAS_AoU_ADJLDL_COMBINED.rds
# Output: lead_snp_summary_table_final.csv
# =============================================================================

library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)
library(ggrepel)

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

# GWS threshold (Bonferroni)
gws <- 0.05 / nrow(combined)
cat("GWS threshold:", gws, "\n")

# Convert LOG10P to p-values
female$P <- 10^(-female$LOG10P)
male$P <- 10^(-male$LOG10P)
combined$P <- 10^(-combined$LOG10P)

# =============================================================================
# 3. CLUMPING: identify independent lead SNPs
# =============================================================================

# Parse chromosome and position from ID (format: chr19:44892652:C:G)
combined_gws <- combined |>
  dplyr::filter(P < gws) |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE) |>
  mutate(chr_num = as.numeric(str_remove(chr, "chr"))) |>
  arrange(P)

# Greedy positional clumping (1MB window)
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

merged_gws <- merged[merged$ID %in% lead_snp_ids, ]
merged_gws$P_diff_adj  <- p.adjust(merged_gws$P_diff, method="bonferroni")
merged_gws$sex_diff_sig <- merged_gws$P_diff_adj < 0.05

cat("Z-diff results:\n")
print(table(merged_gws$sex_diff_sig))

# Inspect significant SNPs
sig_snps <- merged_gws[merged_gws$sex_diff_sig == TRUE, ]
sig_snps$direction <- case_when(
  10^(-sig_snps$LOG10P_fem) < gws & 10^(-sig_snps$LOG10P_male) >= gws ~ "female-specific",
  10^(-sig_snps$LOG10P_male) < gws & 10^(-sig_snps$LOG10P_fem) >= gws ~ "male-specific",
  abs(sig_snps$BETA_fem) > abs(sig_snps$BETA_male) ~ "larger effect in females",
  abs(sig_snps$BETA_male) > abs(sig_snps$BETA_fem) ~ "larger effect in males",
  sign(sig_snps$BETA_fem) != sign(sig_snps$BETA_male) ~ "opposite direction",
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
# 6. GENE ANNOTATION via biomaRt
# =============================================================================

ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

positions <- summary_table |>
  separate(ID, into=c("chr","pos","ref","alt"),
           sep=":", remove=FALSE, convert=TRUE) |>
  mutate(chr_num = as.numeric(str_remove(chr, "chr"))) |>
  arrange(chr_num, pos)

# Narrow window query (10kb)
gene_annotations_per_snp <- map_dfr(1:nrow(positions), function(i) {
  res <- tryCatch(getBM(
    attributes=c("chromosome_name","start_position","end_position",
                 "external_gene_name","gene_biotype"),
    filters=c("chromosome_name","start","end"),
    values=list(
      positions$chr_num[i],
      positions$pos[i] - 10000,
      positions$pos[i] + 10000
    ),
    mart=ensembl
  ), error=function(e) NULL)
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res |> dplyr::mutate(
    chromosome_name = as.character(chromosome_name),
    start_position = as.integer(start_position),
    end_position = as.integer(end_position),
    external_gene_name = as.character(external_gene_name),
    gene_biotype = as.character(gene_biotype),
    SNP_ID = positions$ID[i],
    SNP_pos = as.integer(positions$pos[i])
  )
})

annotations_final <- gene_annotations_per_snp |>
  dplyr::filter(gene_biotype == "protein_coding") |>
  dplyr::mutate(
    distance = pmin(abs(SNP_pos - start_position), abs(SNP_pos - end_position)),
    distance = ifelse(SNP_pos >= start_position & SNP_pos <= end_position, 0, distance)
  ) |>
  dplyr::group_by(SNP_ID) |>
  dplyr::slice_min(order_by=distance, n=1, with_ties=FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(SNP_ID, Gene=external_gene_name, distance)

# Wide window query (500kb) for unannotated SNPs + SMARCA4/LDLR check
wide_snps <- positions |>
  dplyr::filter(ID %in% c("chr2:21090558:G:T",
                          "chr8:125467846:G:A",
                          "chr20:17864040:G:T",
                          "chr19:11080521:G:A"))

gene_annotations_wide <- map_dfr(1:nrow(wide_snps), function(i) {
  res <- tryCatch(getBM(
    attributes=c("chromosome_name","start_position","end_position",
                 "external_gene_name","gene_biotype"),
    filters=c("chromosome_name","start","end"),
    values=list(
      wide_snps$chr_num[i],
      wide_snps$pos[i] - 500000,
      wide_snps$pos[i] + 500000
    ),
    mart=ensembl
  ), error=function(e) NULL)
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res |> dplyr::mutate(
    chromosome_name = as.character(chromosome_name),
    start_position = as.integer(start_position),
    end_position = as.integer(end_position),
    external_gene_name = as.character(external_gene_name),
    gene_biotype = as.character(gene_biotype),
    SNP_ID = wide_snps$ID[i],
    SNP_pos = as.integer(wide_snps$pos[i])
  )
}) |>
  dplyr::filter(gene_biotype == "protein_coding") |>
  dplyr::mutate(
    distance = pmin(abs(SNP_pos - start_position), abs(SNP_pos - end_position)),
    distance = ifelse(SNP_pos >= start_position & SNP_pos <= end_position, 0, distance)
  ) |>
  dplyr::group_by(SNP_ID) |>
  dplyr::slice_min(order_by=distance, n=1, with_ties=FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(SNP_ID, Gene=external_gene_name, distance)

# Patch: replace wide_snps entries with wide window results
annotations_final_patched <- annotations_final |>
  dplyr::filter(!SNP_ID %in% wide_snps$ID) |>
  bind_rows(gene_annotations_wide) |>
  dplyr::group_by(SNP_ID) |>
  dplyr::slice_min(order_by=distance, n=1, with_ties=FALSE) |>
  dplyr::ungroup()

# Manual override: SMARCA4 -> LDLR (SNP sits between genes, LDLR is canonical)
annotations_final_patched <- annotations_final_patched |>
  dplyr::mutate(Gene = ifelse(SNP_ID == "chr19:11080521:G:A", "LDLR", Gene))

cat("Final annotations (should be 14 rows):\n")
print(annotations_final_patched, n=14)

# =============================================================================
# 7. RSID LOOKUP via NCBI dbSNP
# =============================================================================

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

rsid_results <- map_dfr(1:nrow(positions), function(i) {
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
    P_fem = ifelse(P_fem == "0.00e+00", "<1e-300", P_fem)
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

write_csv(summary_table_final, 'lead_snp_summary_table_final.csv')
cat("Final summary table:\n")
print(summary_table_final)

