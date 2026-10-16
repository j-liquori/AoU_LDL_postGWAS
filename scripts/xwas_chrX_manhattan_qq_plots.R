# =============================================================================
# Script: xwas_chrX_manhattan_qq_plots.R
#
# Purpose:
#   Visualize chromosome X GWAS results for LDL-C using Manhattan plots
#   and quantile-quantile (QQ) plots, and calculate genomic inflation
#   factors (lambda) for pooled, female, and male analyses.
#
# Inputs:
#   Three REGENIE chromosome X GWAS results saved as RDS files:
#   - XWAS_AoU_ADJLDL_COMBINED.rds
#   - XWAS_AoU_ADJLDL_FEMALE.rds
#   - XWAS_AoU_ADJLDL_MALE.rds
#
# Outputs:
#   XWAS_combined_manhattan.png
#   XWAS_combined_qq.png
#
# Methods:
#   - Genome-wide significance threshold: 0.05 / 365,774 variants.
#   - Manhattan plots display -log10(p-values) by cumulative genomic position.
#   - QQ plots compare observed and expected -log10(p-values).
#   - Genomic inflation factor (lambda) is calculated from the median
#     observed chi-squared statistic divided by its expected median under
#     the null distribution, assuming 1 degree of freedom.
#
# Notes:
#   - Input files must contain the REGENIE columns used in this script,
#     including CHROM, GENPOS, ID, and LOG10P.
#   - The Bonferroni threshold assumes 365,774 tested variants.
#   - Annotation of rs73461969 is performed for the male Manhattan plot.
#   - This script expects the input RDS files to be available locally.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

lapply(c('tidyverse', 'patchwork', 'ggrepel'),
       function(pkg_name) { if(! pkg_name %in% installed.packages()) 
         { install.packages(pkg_name)} } )

library(tidyverse)
library(patchwork)
library(ggrepel)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

results <- list(
  Pooled = readRDS('XWAS_AoU_ADJLDL_COMBINED.rds'),
  Female = readRDS('XWAS_AoU_ADJLDL_FEMALE.rds'),
  Male = readRDS('XWAS_AoU_ADJLDL_MALE.rds')
)

# =============================================================================
# 2. PREP DATA
# =============================================================================

results <- lapply(results, function(x) {
  
  x |>
    mutate(
      p_val = 10^(-LOG10P),
      RSID = paste0("chr", CHROM, ":", GENPOS)
    )
})

# greatest amount of tested SNPS?

nrow(results$Pooled)
nrow(results$Female)
nrow(results$Male)

# -----------------------------
# Bonferroni-significant SNPs
# -----------------------------

bonferroni_snps <- lapply(results, function(x) {
  
  threshold <- 0.05 / 365774
  
  x[x$p_val < threshold, ]
})

# how many significant in each grp?

sapply(bonferroni_snps, nrow)

# what was the sig male SNP?

bonferroni_snps$Male

# what are the pooled and female p-values of this SNP?

results$Pooled |> filter(ID == 'chrX:22957855:G:A')
results$Female |> filter(ID == 'chrX:22957855:G:A')
results$Male |> filter(ID == 'chrX:22957855:G:A')

# =============================================================================
# 3. MANHATTAN PLOT FUNCTION
# =============================================================================

plot_gg_manhattan <- function(regenie_results,
                              significant_snps,
                              analysis_name) {
  
  bonferroni_threshold <- 0.05 / 365774
  bonferroni_line <- -log10(bonferroni_threshold)
  
  # Set color based on analysis
  analysis_color <- switch(
    tolower(analysis_name),
    "pooled" = "#AA4499",
    "female" = "#332288",
    "male" = "#44AA99",
    stop("analysis_name must be 'pooled', 'female', or 'male'")
  )
  
  # Prepare data
  plot_data <- regenie_results %>%
    mutate(
      CHROM = as.numeric(CHROM),
      GENPOS = as.numeric(GENPOS),
      logp = -log10(p_val)
    ) %>%
    arrange(CHROM, GENPOS)
  
  # Calculate cumulative genomic positions
  chr_offsets <- plot_data %>%
    group_by(CHROM) %>%
    summarise(
      chr_max = max(GENPOS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      offset = cumsum(lag(chr_max, default = 0))
    )
  
  plot_data <- plot_data %>%
    left_join(
      chr_offsets %>% select(CHROM, offset),
      by = "CHROM"
    ) %>%
    mutate(
      BPcum = GENPOS + offset
    )
  
  # Chromosome centers for x-axis
  axis_data <- plot_data %>%
    group_by(CHROM) %>%
    summarise(
      center = mean(range(BPcum)),
      .groups = "drop"
    )
  
  # Significant SNPs
  sig_data <- plot_data %>%
    filter(RSID %in% significant_snps$RSID)
  
  # Base plot
  p <- ggplot(plot_data, aes(x = BPcum, y = logp)) +
    
    geom_point(
      color = analysis_color,
      size = 1.3,
      alpha = 0.8
    ) +
    
    geom_hline(
      yintercept = bonferroni_line,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    scale_x_continuous(
      breaks = seq(
        0,
        max(plot_data$BPcum, na.rm = TRUE),
        by = 20e6
      ),
      labels = function(x) paste0(x / 1e6, " Mb"),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    
    labs(
      x = "Chromosome X",
      y = expression(-log[10](italic(P)))
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  # Highlight significant SNPs
  if (nrow(sig_data) > 0) {
    p <- p +
      geom_point(
        data = sig_data,
        color = analysis_color,
        size = 1.3,
        alpha = 0.8
      )
  }
  
  # Annotate rs73461969 in males only
  if (tolower(analysis_name) == "male") {
    
    annotation_data <- plot_data %>%
      filter(RSID == "chr23:22957855")
    
    if (nrow(annotation_data) > 0) {
      p <- p +
        geom_point(
          data = annotation_data,
          color = analysis_color,
          size = 1.3,
          alpha = 0.8
        ) +
        geom_text_repel(
          data = annotation_data,
          aes(label = "rs73461969"),
          size = 4,
          box.padding = 0.25,
          point.padding = 0.3,
          min.segment.length = 0,
          max.overlaps = Inf
        )
    }
  }
  
  return(p)
}

# =============================================================================
# 4. MAKE PLOTS
# =============================================================================

plots <- list()

for (analysis_name in names(results)) {
  
  plots[[analysis_name]] <- plot_gg_manhattan(
    results[[analysis_name]],
    bonferroni_snps[[analysis_name]],
    analysis_name
  )
}

# =============================================================================
# 5. COMBINE PLOTS
# =============================================================================
options(repr.plot.width = 15, repr.plot.height = 10)

combined_plot <- wrap_plots(
  plots,
  ncol = 1
) +
  plot_annotation(tag_levels = "a",          
                  tag_prefix  = "(",
                  tag_suffix  = ")",
                  theme = theme(
                    plot.tag = element_text(size = 18, face = "bold")
                  )) 

combined_plot

ggsave(
  filename = "XWAS_combined_manhattan.png",
  plot = combined_plot,
  width = 15,
  height = 10,
  units = "in",
  dpi = 600,
  device = "png"
)

# =============================================================================
# 6. LAMBDA
# =============================================================================

calculate_lambda <- function(gwas_results) {
  
  observed_pvalues <- gwas_results$p_val
  
  observed_pvalues <- observed_pvalues[
    is.finite(observed_pvalues) &
      observed_pvalues > 0 &
      observed_pvalues <= 1
  ]
  
  observed_chisq <- qchisq(
    1 - observed_pvalues,
    df = 1
  )
  
  median(observed_chisq) / qchisq(0.5, df = 1)
}

lambda <- sapply(
  results,
  calculate_lambda
)

lambda

# =============================================================================
# 7. QQ PLOT FUNCTION
# =============================================================================

plot_gg_qq <- function(regenie_results, analysis_name) {
  
  # Set color based on analysis
  analysis_color <- switch(
    tolower(analysis_name),
    "pooled" = "#AA4499",
    "female" = "#332288",
    "male" = "#44AA99",
    stop("analysis_name must be 'pooled', 'female', or 'male'")
  )
  
  # Calculate lambda
  lambda <- calculate_lambda(regenie_results)
  
  # Prepare p-values
  observed_pvalues <- regenie_results$p_val
  
  observed_pvalues <- observed_pvalues[
    is.finite(observed_pvalues) &
      observed_pvalues > 0 &
      observed_pvalues <= 1
  ]
  
  # Expected and observed -log10 p-values
  n <- length(observed_pvalues)
  
  expected <- -log10(
    ppoints(n)
  )
  
  observed <- -log10(
    sort(observed_pvalues)
  )
  
  qq_data <- data.frame(
    expected = expected,
    observed = observed
  )
  
  # Create QQ plot
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    
    geom_point(
      color = analysis_color,
      size = 1.3,
      alpha = 0.8
    ) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = paste0(
        "lambda == ",
        round(lambda, 2)
      ),
      hjust = -0.1,
      vjust = 1.5,
      size = 5,
      fontface = "bold",
      parse = TRUE
    ) +
    
    labs(
      x = expression(Expected~~-log[10](P)),
      y = expression(Observed~~-log[10](P))
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  return(p)
}

# =============================================================================
# 8. PLOT
# =============================================================================

qq_plots <- list()

for (analysis_name in names(results)) {
  
  qq_plots[[analysis_name]] <- plot_gg_qq(
    results[[analysis_name]],
    analysis_name
  )
}

# =============================================================================
# 9. COMBINE
# =============================================================================

options(repr.plot.width = 10, repr.plot.height = 10)

combined_qq <- wrap_plots(
  qq_plots,
  ncol = 1
) +
  plot_annotation(tag_levels = "a",          # lowercase a, b, c tags
                  tag_prefix  = "(",
                  tag_suffix  = ")",
                  theme = theme(
                    plot.tag = element_text(size = 18, face = "bold")
                  )) 

combined_qq

ggsave(
  filename = "XWAS_combined_qq.png",
  plot = combined_qq,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600,
  device = "png"
)