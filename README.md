# TB Whole-Blood Transcriptomics: GSE79362 & GSE94438

Analysis of two prospective whole-blood RNA-seq cohorts from `curatedTBData`, comparing active/progressive TB against latent infection or healthy household contacts.

## Datasets

| Dataset | Cohort | Contrast | N |
|---------|--------|----------|---|
| GSE79362 | South African adults (Zak et al. 2016) | PTB vs LTBI | 355 samples (110 PTB / 245 LTBI) |
| GSE94438 | Pan-African household contacts (Suliman et al. 2018) | PTB/progressors vs controls | 428 labelled samples (101 PTB / 327 Control) |

## Key Results

- **GSE79362 baseline DEGs:** 30 (FDR < 0.05, |logFC| > 1); 29 up in PTB, 1 down
- **GSE79362 adjusted DEGs:** 9 after adjusting for timepoint and patient-level repeated measures using limma-voom duplicateCorrelation
- **GSE94438 baseline DEGs:** 43 (FDR < 0.05, |logFC| > 1); all 43 up in PTB/progressors
- **GSE94438 adjusted DEGs:** 41 after adjustment for site/country, sex and age

| Signature | GSE79362 AUC | GSE94438 AUC |
|-----------|-------------|-------------|
| Zak16 | 0.765 | 0.686 |
| RISK4 | 0.759 | 0.687 |
| Eleven_gene | 0.771 | 0.694 |

Cross-dataset logFC Spearman r = 0.641 · Hallmark NES Spearman r = 0.815  
Shared tested genes = 14,128 · Shared strict DEGs = 15 · Top-100 DEG overlap = 40

## Biological Conclusion

Across both whole-blood RNA-seq cohorts, active or progressive TB is marked by a reproducible **interferon–inflammatory–complement–myeloid programme**. This signal distinguishes active disease or progression from LTBI/exposed non-progression, but it likely reflects both immune activation and shifts in circulating leukocyte composition.

## Repository Structure

```
code/              R analysis scripts (English)
code_original/     Original R scripts
results/           CSV output tables (DEG, GSEA, signature AUC, cross-dataset)
figures/           All output plots (PNG)
RESULTS_SUMMARY.md Detailed results narrative
```

## Running the Analysis

Requires R ≥ 4.2 with Bioconductor packages: `curatedTBData`, `edgeR`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `msigdbr`.

```r
# Install
BiocManager::install(c("curatedTBData","edgeR","limma","clusterProfiler",
                       "org.Hs.eg.db","ReactomePA"))
install.packages(c("msigdbr","pROC","ggplot2","dplyr","pheatmap"))

# Run main analysis
source("code/tb_analysis.R")
```

## Limitations

- Whole-blood RNA-seq cannot distinguish cell-intrinsic regulation from leukocyte composition changes without deconvolution.
- GSE79362 contains repeated longitudinal samples; repeated-measures correction substantially reduces the strict DEG count.
- GSE94438 spans three African sites; site is a required covariate.
- The two datasets use different label definitions and should not be treated as identical biological contrasts.
