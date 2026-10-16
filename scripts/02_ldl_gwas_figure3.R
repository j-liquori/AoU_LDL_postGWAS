# =============================================================================
# LDL-C GWAS Figure 3: Sex-Stratified Miami Plot and Effect-Size Comparison
#
# Purpose:
# Visualize sex-stratified GWAS results and compare female and male effect
# estimates at selected loci.
#
# Panels:
#   (A) Miami plot of female and male GWAS association results.
#   (B) Scatterplot comparing male and female effect estimates at selected loci.
#
# Inputs:
#   GWAS summary statistics and selected-locus results from the preceding
#   analysis, including gene annotations and sex-difference test results.
#
# Outputs:
#   miami.png: Sex-stratified Miami plot.
#   effect_size_comparison.png: Scatterplot of male versus female 
#   effect estimates.
#
# Notes:
#   Male -log10(p) values are plotted below zero to create the Miami layout.
#   Highlighted variants are those classified as sex-differentiated by the
#   preceding z-score difference analysis.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

# Install required packages if they are not already available.
required_packages <- c("tidyverse", "ggrepel", "ggrastr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(ggrepel)
library(ggrastr)

# Load results from the preceding sex-difference analysis.
analysis_results <- readRDS("ldl_sex_difference_results.rds")

female <- analysis_results$female
male <- analysis_results$male
gws <- analysis_results$gws
lead_snp_ids <- analysis_results$lead_snp_ids
sig_snps <- analysis_results$sig_snps
merged_gws <- analysis_results$merged_gws
annotations_final_patched <- analysis_results$annotations_final_patched
rsid_results <- analysis_results$rsid_results

gws_log10p <- -log10(gws)

# =============================================================================
# 1. DATA PREP
# =============================================================================

# Function to prepare GWAS data for Miami plot
prep_manhattan <- function(df, sex_label) {
  df |>
    dplyr::filter(!is.na(LOG10P)) |>
    separate(ID, into=c("chr","pos","ref","alt"),
             sep=":", remove=FALSE, convert=TRUE) |>
    mutate(
      chr_num = as.numeric(str_remove(chr, "chr")),
      pos = as.numeric(pos),
      sex = sex_label
    ) |>
    dplyr::filter(!is.na(chr_num), chr_num %in% 1:22) |>
    arrange(chr_num, pos)
}

female_plot <- prep_manhattan(female, "Female")
male_plot <- prep_manhattan(male, "Male")
combined_plot <- bind_rows(female_plot, male_plot)

# Calculate cumulative chromosome offsets from the maximum observed position
# in the plotted variants, then use chromosome midpoints for axis labels.
chr_lengths <- combined_plot |>
  dplyr::group_by(chr_num) |>
  dplyr::summarise(max_pos = max(pos), .groups="drop") |>
  arrange(chr_num) |>
  mutate(cum_start = lag(cumsum(as.numeric(max_pos)), default=0))

combined_plot <- combined_plot |>
  left_join(chr_lengths |> dplyr::select(chr_num, cum_start), by="chr_num") |>
  mutate(cum_pos = pos + cum_start)

# Chromosome midpoints for x-axis labels
chr_mids <- chr_lengths |>
  mutate(mid = cum_start + max_pos / 2)

# =============================================================================
# 2. LABEL PREP
# =============================================================================

# Prepare labels for selected loci that reach genome-wide significance
# in at least one sex-stratified GWAS. Retain separate rows for female
# and male results so labels can be placed on the appropriate side
# of the Miami plot.
# Alternating colors by chromosome and sex
chr_colors <- c(
  "Female_light" = "#88CCEE",
  "Female_dark"  = "#332288",
  "Male_light"   = "#44AA99",
  "Male_dark"    = "#117733"
)

combined_plot <- combined_plot |>
  mutate(
    chr_color = case_when(
      sex == "Female" & chr_num %% 2 == 1 ~ "Female_light",
      sex == "Female" & chr_num %% 2 == 0 ~ "Female_dark",
      sex == "Male"   & chr_num %% 2 == 1 ~ "Male_light",
      sex == "Male"   & chr_num %% 2 == 0 ~ "Male_dark"
    )
  )

# Flip male LOG10P for Miami plot
plot_data <- combined_plot |>
  mutate(LOG10P_plot = ifelse(sex == "Male", -LOG10P, LOG10P))


# Y axis limits
y_max <- max(plot_data$LOG10P_plot, na.rm=TRUE) * 1.15
y_min <- min(plot_data$LOG10P_plot, na.rm=TRUE) * 1.15

cat("Rows in lead_snp_positions:", nrow(lead_snp_positions), "\n")
cat("Rows in lead_snp_female:", nrow(lead_snp_female), "\n")
cat("Rows in lead_snp_male:", nrow(lead_snp_male), "\n")

# Check if IDs are matching between combined_plot and lead_snp_ids
cat("Lead SNP IDs in combined_plot:\n")
print(sum(combined_plot$ID %in% lead_snp_ids))

# Check LOG10P range in combined_plot
cat("LOG10P range:\n")
print(summary(combined_plot$LOG10P))

# Check what lead_snp_positions looks like before filtering
lead_snp_positions_debug <- combined_plot |>
  dplyr::filter(ID %in% lead_snp_ids) |>
  left_join(annotations_final_patched, by=c("ID"="SNP_ID")) |>
  dplyr::select(ID, sex, LOG10P, Gene)

cat("Lead SNPs found in plot data:\n")
print(lead_snp_positions_debug)


# Fix label preparation - check GWS in either sex separately
lead_snp_positions <- lead_snp_positions_debug |>
  left_join(
    lead_snp_positions_debug |>
      dplyr::group_by(ID) |>
      dplyr::summarise(max_log10p = max(LOG10P), .groups="drop"),
    by="ID"
  ) |>
  dplyr::filter(max_log10p >= gws_log10p) |>
  left_join(combined_plot |>
              dplyr::select(ID, sex, cum_pos) |>
              distinct(),
            by=c("ID","sex")) |>
  mutate(is_sex_diff = ID %in% sig_snps$ID)

cat("Lead SNP positions after fix:\n")
print(lead_snp_positions |> dplyr::select(ID, sex, LOG10P, Gene, max_log10p))

# Re-split into female and male
lead_snp_female <- lead_snp_positions |> dplyr::filter(sex == "Female")
lead_snp_male   <- lead_snp_positions |> dplyr::filter(sex == "Male") |>
  mutate(LOG10P = -LOG10P)

cat("Female label rows:", nrow(lead_snp_female), "\n")
cat("Male label rows:", nrow(lead_snp_male), "\n")

# =============================================================================
# 3. MIAMI PLOT 
# =============================================================================


miami <- ggplot(plot_data, aes(x=cum_pos, y=LOG10P_plot)) +
  
  # Rasterized points - much faster to render
  rasterise(geom_point(aes(color=chr_color), size=0.4, alpha=0.6), dpi=300) +
  scale_color_manual(values = chr_colors, guide = "none") +
  
  # GWS lines
  geom_hline(yintercept= gws_log10p, linetype="dashed",
             color="red", linewidth=0.4) +
  geom_hline(yintercept=-gws_log10p, linetype="dashed",
             color="red", linewidth=0.4) +
  
  # Dividing line
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  
  # Labels - female panel (deduplicated)
  geom_text_repel(
    data=lead_snp_female |>
      dplyr::arrange(desc(LOG10P)) |>
      dplyr::distinct(Gene, .keep_all=TRUE),  # one label per gene
    aes(x=cum_pos, y=LOG10P,
        label=Gene,
        fontface=ifelse(is_sex_diff, "bold", "plain")),
    color=ifelse(
      (lead_snp_female |>
         dplyr::arrange(desc(LOG10P)) |>
         dplyr::distinct(Gene, .keep_all=TRUE))$is_sex_diff,
      "#AA4499", "black"),
    size=2.5,
    max.overlaps=20,
    segment.size=0.3,
    segment.color="grey50",
    box.padding=0.4,
    point.padding=0.2,
    min.segment.length=0,
    nudge_y=3
  ) +
  
  # Labels - male panel (deduplicated)
  geom_text_repel(
    data=lead_snp_male |>
      dplyr::arrange(LOG10P) |>
      dplyr::distinct(Gene, .keep_all=TRUE),  # one label per gene
    aes(x=cum_pos, y=LOG10P,
        label=Gene,
        fontface=ifelse(is_sex_diff, "bold", "plain")),
    color=ifelse(
      (lead_snp_male |>
         dplyr::arrange(LOG10P) |>
         dplyr::distinct(Gene, .keep_all=TRUE))$is_sex_diff,
      "#AA4499", "black"),
    size=2.5,
    max.overlaps=20,
    segment.size=0.3,
    segment.color="grey50",
    box.padding=0.4,
    point.padding=0.2,
    min.segment.length=0,
    nudge_y=-3
  ) +
  
  # X axis
  scale_x_continuous(
    breaks=chr_mids$mid,
    labels=chr_mids$chr_num,
    expand=c(0.01, 0)
  ) +
  
  # Y axis
  scale_y_continuous(
    breaks=seq(ceiling(y_min), floor(y_max), by=20),
    labels=abs(seq(ceiling(y_min), floor(y_max), by=20)),
    limits=c(y_min, y_max)
  ) +
  
  # Female/Male labels centered above and below
  annotate("text",
           x=max(plot_data$cum_pos) / 2,  # centered on x axis
           y=y_max * 0.97,                # just below top
           label="Female",
           fontface="bold",
           hjust=0.5, size=4) +
  annotate("text",
           x=max(plot_data$cum_pos) / 2,  # centered on x axis
           y=y_min * 0.97,                # just above bottom
           label="Male",
           fontface="bold",
           hjust=0.5, size=4) +
  
  labs(x="Chromosome",
       y=expression(-log[10](p))) +
  theme_classic() +
  theme(
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10),
    panel.grid=element_blank()
  )


# =============================================================================
# 4. EFFECT SCATTERPLOT
# =============================================================================

# Build label dataframe by joining merged_gws with rsid and gene info

plot_labels <- merged_gws |>
  left_join(annotations_final_patched, by = c("ID" = "SNP_ID")) |>
  mutate(
    label = case_when(
      Gene == "APOE" ~ "APOE/APOC1",
      TRUE ~ Gene
    )
  ) |>
  filter(!is.na(label), label != "")

effect <- ggplot(plot_labels, aes(x=BETA_male, y=BETA_fem)) +
  geom_point(aes(color=sex_diff_sig), size=2, alpha=0.8) +
  scale_color_manual(values=c("FALSE"="grey60", "TRUE"="#AA4499"),
                     labels=c("Not significant", "Sex-differentiated"),
                     name="Z-diff test") +
  geom_text_repel(
    aes(label=label),
    size=3,
    max.overlaps=20,
    segment.color="grey50",
    segment.size=0.3,
    # Make sex-diff label bold and red
    fontface=ifelse(plot_labels$sex_diff_sig, "bold", "plain"),
    color=ifelse(plot_labels$sex_diff_sig, "#AA4499", "grey30")
  ) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x="Male Effect Size",
       y="Female Effect Size") +
  theme_classic() +
  theme(legend.position="bottom")

# =============================================================================
# 5. SAVE PLOTS
# =============================================================================

ggsave(
  "miami.png",
  plot = miami,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  "effect_size_comparison.png",
  plot = effect,
  width = 7,
  height = 6,
  dpi = 300
)






