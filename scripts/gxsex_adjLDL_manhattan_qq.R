# =============================================================================
# Script: gxsex_adjLDL_manhattan_qq.R
#
# Purpose:
#   Visualize SNP-by-sex interaction results from the GxSex GWAS of
#   statin-adjusted LDL-C using a Manhattan plot and a quantile-quantile
#   (QQ) plot, and calculate the genomic inflation factor (lambda).
#
# Input:
#   - GxSex_AoU_ADJLDL.rds
#
# Outputs:
#   - gxsex_adjLDL_man.png
#   - gxsex_adjLDL_qq.png
#
# Methods:
#   - Uses the REGENIE SNP-by-sex interaction test
#     (ADD-INT_SNPxsex_at_birth=Male), with female sex as the reference.
#   - Applies a Bonferroni significance threshold based on the number of
#     autosomal interaction tests retained after filtering.
#   - Calculates lambda from the median chi-squared statistic, assuming
#     1 degree of freedom.
#
# Notes:
#   - Input data must contain the REGENIE columns used in this script,
#     including TEST, ID, LOG10P, BETA, and SE.
#   - The analysis is restricted to autosomes 1–22.
#   - Chromosome positions are cumulative offsets based on the maximum
#     observed variant position in each chromosome.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

lapply(c('tidyverse', 'ggrepel'),
       function(pkg_name) { if(! pkg_name %in% installed.packages()) 
         { install.packages(pkg_name)} } )

library(tidyverse)
library(ggrepel)

# =============================================================================
# 1. DATA PREP - GxSex statin-adjusted LDL-C
# =============================================================================

gxsex <- readRDS('GxSex_AoU_ADJLDL.rds')

gxsex_plot <- gxsex |>
  filter(
    # female sex is the reference
    TEST == "ADD-INT_SNPxsex_at_birth=Male",
    !is.na(LOG10P)
  ) |>
  separate(
    ID,
    into = c("chr", "pos", "ref", "alt"),
    sep = ":",
    remove = FALSE,
    convert = TRUE
  ) |>
  mutate(
    chr_num = as.numeric(str_remove(chr, "chr")),
    pos = as.numeric(pos),
    P = 10^(-LOG10P)
  ) |>
  filter(
    !is.na(chr_num),
    chr_num %in% 1:22,
    !is.na(pos)
  ) |>
  arrange(chr_num, pos)

# =============================================================================
# 2. CUMULATIVE GENOMIC POSITIONS - GxSex statin-adjusted LDL-C
# =============================================================================

chr_lengths <- gxsex_plot |>
  group_by(chr_num) |>
  summarise(
    max_pos = max(pos),
    .groups = "drop"
  ) |>
  arrange(chr_num) |>
  mutate(
    cum_start = lag(cumsum(as.numeric(max_pos)), default = 0)
  )

gxsex_plot <- gxsex_plot |>
  left_join(
    chr_lengths |>
      select(chr_num, cum_start),
    by = "chr_num"
  ) |>
  mutate(
    cum_pos = pos + cum_start
  )

# Chromosome midpoints for x-axis labels
chr_mids <- chr_lengths |>
  mutate(
    mid = cum_start + max_pos / 2
  )

# =============================================================================
# 3. GENOME-WIDE SIGNIFICANCE THRESHOLD - GxSex statin-adjusted LDL-C
# =============================================================================

gws <- 0.05 / nrow(gxsex_plot)

gws_log10p <- -log10(gws)

cat("Number of interaction tests:", nrow(gxsex_plot), "\n")
cat("GWS threshold:", gws, "\n")
cat("-log10(GWS):", gws_log10p, "\n")

# =============================================================================
# 4. LABEL PREP - GxSex statin-adjusted LDL-C
# =============================================================================

# Label all GWS SNPs with their variant ID.
label_data <- gxsex_plot |>
  filter(LOG10P >= gws_log10p) |>
  mutate(
    label = ID
  )

cat("Number of GWS SNPs:", nrow(label_data), "\n")

cat("\nGWS SNPs:\n")
print(
  label_data |>
    select(
      Chr = chr_num,
      Position = pos,
      ID,
      BETA,
      SE,
      LOG10P,
      P
    )
)

# =============================================================================
# 5. CHROMOSOME COLORS - GxSex statin-adjusted LDL-C
# =============================================================================

chr_colors <- c(
  "light" = "#44AA99",
  "dark"  = "#332288"
)

gxsex_plot <- gxsex_plot |>
  mutate(
    chr_color = ifelse(
      chr_num %% 2 == 0,
      "dark",
      "light"
    )
  )

# =============================================================================
# 6. MANHATTAN PLOT - GxSex statin-adjusted LDL-C
# =============================================================================
options(repr.plot.width = 15, repr.plot.height = 10)

manhattan_gxsex <- ggplot(
  gxsex_plot,
  aes(
    x = cum_pos,
    y = LOG10P
  )
) +
  
  geom_point(
    aes(color = chr_color),
    size = 0.4,
    alpha = 0.6
  ) +
  
  scale_color_manual(
    values = chr_colors,
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = gws_log10p,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  
  geom_point(
    data = label_data,
    aes(
      x = cum_pos,
      y = LOG10P
    ),
    size = 1.5,
    color = "#AA4499"
  ) +
  
  geom_text_repel(
    data = label_data |>
      arrange(desc(LOG10P)),
    aes(
      x = cum_pos,
      y = LOG10P,
      label = label
    ),
    size = 3.5,
    color = "#AA4499",
    max.overlaps = 20,
    segment.size = 0.25,
    segment.color = "grey50",
    box.padding = 0.15,
    point.padding = 0.1,
    min.segment.length = 0,
    nudge_y = 0.5,
    force = 0.5,
    force_pull = 1
  ) +
  
  scale_x_continuous(
    breaks = chr_mids$mid,
    labels = chr_mids$chr_num,
    expand = c(0.01, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(max(gxsex_plot$LOG10P, na.rm = TRUE)),
      by = 2
    ),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid = element_blank()
  )

manhattan_gxsex

ggsave(
  filename = "gxsex_adjLDL_man.png",
  plot = manhattan_gxsex,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600
)

# =============================================================================
# 7. LAMBDA - GxSex statin-adjusted LDL-C
# =============================================================================

# Get SNP × sex interaction p-values
qq_data <- gxsex |>
  filter(TEST == "ADD-INT_SNPxsex_at_birth=Male") |>
  mutate(
    P = 10^(-LOG10P)
  ) |>
  filter(
    is.finite(P),
    P > 0,
    P <= 1
  ) |>
  arrange(P) |>
  mutate(
    expected = -log10(ppoints(n())),
    observed = -log10(P)
  )

# Calculate lambda
chisq <- qchisq(1 - qq_data$P, df = 1)

lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

lambda

# =============================================================================
# 8. QQ PLOT - GxSex statin-adjusted LDL-C
# =============================================================================

options(repr.plot.width = 10, repr.plot.height = 10)

qqsex <- ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 0.3, alpha = 0.3, color = "#AA4499") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("Lambda = %.2f", lambda),
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  labs(
    x = expression(Expected~~-log[10](p)),
    y = expression(Observed~~-log[10](p))
  ) +
  theme_classic(base_size = 14)

qqsex

ggsave(
  filename = "gxsex_adjLDL_qq.png",
  plot = qqsex,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)