# AoU_LDL_postGWAS
Post-GWAS analysis scripts for [*Phenotypic and genetic sex differences in LDL-C and statin treatment in All of Us*](https://doi.org/10.64898/2026.01.26.26344867). Includes annotation, sex-stratified effect size comparisons, and visualization of GWAS, XWAS, and GxSex summary statistics of statin-adjusted LDL-C generated using REGENIE v2.2.4 and v3.1.4 from the All of Us Research Program (v7).

## Analysis and execution order
### Post GWAS analysis
The numbered scripts form a sequential workflow and should be run in ascending numerical order.

| Script | Purpose | 
|--------|---------|
| `00_annotate_gws_snps.R` | Identifies genome-wide significant variants in the female, male, and combined GWAS results and annotates variants with rsIDs and nearby genes.|
| `01_ldl_sex_difference_zscore.R` | Compares female and male effect estimates at selected loci using a z-score difference test with Bonferroni correction.|
| `02_ldl_gwas_figure3.R` | Generates a sex-stratified Miami plot and an effect-size comparison plot using GWAS results and the preceding analysis outputs.|

### XWAS
`xwas_chrX_manhattan_qq_plots.R`

Generates Manhattan and quantile-quantile (QQ) plots for combined, female, and male chromosome X association results and calculates genomic inflation factors.

This script can be run independently of the other repository scripts, provided the required XWAS input files are available.

### GxSex

`gxsex_adjLDL_manhattan_qq.R`

Generates Manhattan and QQ plots for SNP-by-sex interaction results and calculates the genomic inflation factor. The analysis visualizes REGENIE interaction results using female sex as the reference and is restricted to autosomes 1–22.

This script can be run independently of the other repository scripts, provided the required GxSex input file is available.

## Data Availability
The summary statistics used in this analysis are available on Zenodo:  
[https://doi.org/10.5281/zenodo.19333986](https://doi.org/10.5281/zenodo.19333986)

Individual-level data, including genomic and phenotype data, from the All of Us Research Program are not included in this repository. Access to restricted data requires authorization through the All of Us Researcher Workbench.

This repository currently contains post-GWAS, XWAS, and GxSex analyses. Additional analyses, including heritability summaries, may be incorporated in future updates.

## References

**Manuscript**

Liquori, J. L., Tian, X., & Musharoff, S. A. (2026). Phenotypic and genetic sex differences in LDL-C and statin treatment in All of Us, medRxiv 2026.01.26.26344867; doi: https://doi.org/10.64898/2026.01.26.26344867

**Summary statistics**

Liquori, J. L., Tian, X., & Musharoff, S. A. (2026). Phenotypic and genetic sex differences in LDL-C and statin treatment in All of Us [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.23415946
