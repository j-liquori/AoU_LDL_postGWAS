# =============================================================================
# LDL-C GWS SNP Annotation
# Annotates genome-wide significant SNPs from female, male, and combined GWAS
# with rsIDs (NCBI dbSNP) and nearest protein-coding gene (Ensembl/biomaRt)
# Output: gws_female_annotated.csv, gws_male_annotated.csv,
#         gws_combined_annotated.csv
# =============================================================================

library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

female <- readRDS('GWAS_AoU_ADJLDL_FEMALE.rds')
male <- readRDS('GWAS_AoU_ADJLDL_MALE.rds')
combined <- readRDS('GWAS_AoU_ADJLDL_COMBINED.rds')

# GWS threshold (Bonferroni: 0.05 / number of variants in combined GWAS)
gws <- 0.05 / nrow(combined)
cat("GWS threshold:", gws, "\n")

# Convert LOG10P to p-values
female$P <- 10^(-female$LOG10P)
male$P <- 10^(-male$LOG10P)
combined$P <- 10^(-combined$LOG10P)

# =============================================================================
# 2. FILTER TO GWS SNPs PER ANALYSIS
# =============================================================================

female_gws <- female |> dplyr::filter(P < gws)
male_gws <- male |> dplyr::filter(P < gws)
combined_gws <- combined |> dplyr::filter(P < gws)

cat("GWS SNPs - Female:", nrow(female_gws),
    "| Male:", nrow(male_gws),
    "| Combined:", nrow(combined_gws), "\n")

# =============================================================================
# 3. HELPER: PARSE POSITION FROM ID
# =============================================================================

parse_positions <- function(df) {
  df |>
    separate(ID, into=c("chr","pos","ref","alt"),
             sep=":", remove=FALSE, convert=TRUE) |>
    mutate(
      chr_num = as.numeric(str_remove(chr, "chr")),
      pos = as.integer(pos)
    ) |>
    dplyr::filter(!is.na(chr_num), chr_num %in% 1:22) |>
    arrange(chr_num, pos)
}

female_gws <- parse_positions(female_gws)
male_gws <- parse_positions(male_gws)
combined_gws <- parse_positions(combined_gws)

# =============================================================================
# 4. HELPER: GENE ANNOTATION VIA biomaRt
# (narrow 10kb window, 500kb fallback for unannotated SNPs)
# =============================================================================

ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

query_biomart <- function(chr_num, pos, window, mart) {
  res <- tryCatch(getBM(
    attributes = c("chromosome_name", "start_position", "end_position",
                   "external_gene_name", "gene_biotype"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chr_num, pos - window, pos + window),
    mart = mart
  ), error = function(e) NULL)
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # Coerce all columns to prevent bind_rows type conflicts
  res |> mutate(
    chromosome_name = as.character(chromosome_name),
    start_position = as.integer(start_position),
    end_position = as.integer(end_position),
    external_gene_name = as.character(external_gene_name),
    gene_biotype = as.character(gene_biotype)
  )
}

get_nearest_gene <- function(results_df, snp_id, snp_pos) {
  results_df |>
    dplyr::filter(gene_biotype == "protein_coding") |>
    mutate(
      SNP_ID = snp_id,
      SNP_pos = as.integer(snp_pos),
      distance = pmin(abs(SNP_pos - start_position),
                      abs(SNP_pos - end_position)),
      distance = ifelse(SNP_pos >= start_position &
                          SNP_pos <= end_position, 0L, distance)
    ) |>
    dplyr::group_by(SNP_ID) |>
    dplyr::slice_min(order_by=distance, n=1, with_ties=FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(SNP_ID, Gene=external_gene_name, distance)
}

annotate_snps <- function(positions_df, mart) {
  cat("Querying narrow window (10kb) for", nrow(positions_df), "SNPs...\n")
  
  narrow <- map_dfr(1:nrow(positions_df), function(i) {
    res <- query_biomart(positions_df$chr_num[i],
                         positions_df$pos[i],
                         window=10000, mart=mart)
    if (is.null(res)) return(NULL)
    get_nearest_gene(res, positions_df$ID[i], positions_df$pos[i])
  })
  
  # SNPs with no annotation from narrow window
  unannotated_ids <- positions_df$ID[!positions_df$ID %in% narrow$SNP_ID]
  cat("Unannotated after narrow window:", length(unannotated_ids),
      "- running 500kb fallback...\n")
  
  if (length(unannotated_ids) > 0) {
    wide_pos <- positions_df |> dplyr::filter(ID %in% unannotated_ids)
    
    wide <- map_dfr(1:nrow(wide_pos), function(i) {
      res <- query_biomart(wide_pos$chr_num[i],
                           wide_pos$pos[i],
                           window=500000, mart=mart)
      if (is.null(res)) return(NULL)
      get_nearest_gene(res, wide_pos$ID[i], wide_pos$pos[i])
    })
    
    bind_rows(narrow, wide)
  } else {
    narrow
  }
}

# =============================================================================
# 5. HELPER: rsID LOOKUP VIA NCBI dbSNP eutils
# =============================================================================

query_dbsnp <- function(chrom, pos) {
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=snp",
    "&term=", chrom, "[CHR]+AND+", pos, "[CHRPOS38]",
    "&retmode=json&retmax=5"
  )
  res <- tryCatch(GET(url, timeout(10)), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(NA_character_)
  
  content <- fromJSON(rawToChar(res$content))
  ids <- content$esearchresult$idlist
  
  if (length(ids) == 0) return(NA_character_)
  paste0("rs", ids[1])
}

lookup_rsids <- function(positions_df) {
  cat("  Looking up rsIDs for", nrow(positions_df), "SNPs (~",
      round(nrow(positions_df) * 0.5 / 60, 1), "min)...\n")
  
  map_chr(1:nrow(positions_df), function(i) {
    if (i %% 50 == 0) cat("  ... SNP", i, "of", nrow(positions_df), "\n")
    rsid <- query_dbsnp(positions_df$chr_num[i], positions_df$pos[i])
    Sys.sleep(0.5)
    rsid
  })
}

# =============================================================================
# 6. HELPER: BUILD FINAL ANNOTATED TABLE
# =============================================================================

build_annotated_table <- function(gws_df, label) {
  cat("\n--- Annotating", label, "GWAS (", nrow(gws_df), "GWS SNPs) ---\n")
  
  cat("Step 1: Gene annotation\n")
  gene_annot <- annotate_snps(gws_df, ensembl)
  
  cat("Step 2: rsID lookup\n")
  rsids <- lookup_rsids(gws_df)
  
  cat("Step 3: Building output table\n")
  gws_df |>
    mutate(
      rsID = rsids,
      P = formatC(P, format="e", digits=2),
      P = ifelse(P == "0.00e+00", "<1e-300", P)
    ) |>
    left_join(gene_annot, by=c("ID"="SNP_ID")) |>
    dplyr::select(
      Chr = chr_num,
      Position = pos,
      ID,
      rsID,
      Gene,
      BETA,
      SE,
      LOG10P,
      P
    ) |>
    arrange(Chr, Position)
}

# =============================================================================
# 7. RUN ANNOTATION FOR EACH GWAS
# =============================================================================

female_annotated <- build_annotated_table(female_gws, "Female")
male_annotated <- build_annotated_table(male_gws, "Male")
combined_annotated <- build_annotated_table(combined_gws, "Combined")

# =============================================================================
# 8. WRITE OUTPUT
# =============================================================================

write_csv(female_annotated, "gws_female_annotated.csv")
write_csv(male_annotated, "gws_male_annotated.csv")
write_csv(combined_annotated, "gws_combined_annotated.csv")

cat("\nDone.\n")
cat("Female GWS SNPs:", nrow(female_annotated), "\n")
cat("Male GWS SNPs:", nrow(male_annotated), "\n")
cat("Combined GWS SNPs:", nrow(combined_annotated), "\n")


# =============================================================================
# 9. ADD CROSS-SEX GWS FLAG
# =============================================================================

female_annotated <- read_csv("gws_female_annotated.csv")
male_annotated <- read_csv("gws_male_annotated.csv")

female_annotated <- female_annotated |>
  mutate(Also_GWS_in_males = ID %in% male_annotated$ID)

male_annotated <- male_annotated |>
  mutate(Also_GWS_in_females = ID %in% female_annotated$ID)

write_csv(female_annotated, "gws_female_annotated.csv")
write_csv(male_annotated,   "gws_male_annotated.csv")

cat("Cross-sex GWS flags added.\n")
cat("Female SNPs also GWS in males:", sum(female_annotated$Also_GWS_in_males), "\n")
cat("Male SNPs also GWS in females:", sum(male_annotated$Also_GWS_in_females), "\n")