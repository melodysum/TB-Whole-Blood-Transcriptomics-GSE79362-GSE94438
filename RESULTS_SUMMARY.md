# TB whole-blood transcriptomics analysis: GSE79362 and GSE94438

This folder contains the English deliverable package for the curatedTBData analysis of two TB whole-blood RNA-seq datasets.

## Dataset interpretation

- **GSE79362** compares **PTB versus LTBI** in a South African cohort. It is best interpreted as an active-disease signal on a latent-infection background.
- **GSE94438** compares **PTB/progressors versus household-contact controls** across Ethiopia, South Africa, and The Gambia. It is best interpreted as a progression or disease-emergence signal among exposed contacts.
- The GSE94438 control group should not be described as strictly uninfected without additional TST/QFT-based stratification.

## Key sample counts

- **GSE79362:** 355 labelled samples: 110 PTB and 245 LTBI.
- **GSE94438:** 434 hg38 count samples in curatedTBData; 428 samples have TBStatus labels: 101 PTB and 327 Control.

## Main differential expression results

- **GSE79362 baseline edgeR model:** 30 strict DEGs at FDR < 0.05 and |logFC| > 1; 29 up-regulated and 1 down-regulated in PTB.
- **GSE79362 adjusted model:** 9 strict DEGs after adjusting for timepoint and patient-level repeated measures with duplicateCorrelation.
- **GSE94438 baseline edgeR model:** 43 strict DEGs; all 43 are up-regulated in PTB/progressors.
- **GSE94438 adjusted model:** 41 strict DEGs after adjustment for site/country, sex, and age.

## Signature validation

- **Zak16 AUC:** 0.765 in GSE79362 and 0.686 in GSE94438.
- **RISK4 AUC:** 0.759 in GSE79362 and 0.687 in GSE94438.
- **11-gene AUC:** 0.771 in GSE79362 and 0.694 in GSE94438.

The signature scores perform better in the PTB/LTBI contrast than in the household-contact progression contrast, which is biologically expected because the latter is a more heterogeneous and clinically harder prediction task.

## Cross-dataset consistency

- Shared tested genes: 14,128.
- Spearman correlation of gene-level logFC values: 0.641.
- Shared strict DEGs: 15.
- Top-100 DEG overlap: 40 genes.
- Spearman correlation of Hallmark pathway NES values: 0.815.

Pathway-level agreement is stronger than individual-gene overlap, indicating that the two datasets converge on shared TB biology even when the exact DEG lists differ.

## Biological conclusion

Across both prospective whole-blood RNA-seq cohorts, active or progressive TB is marked by a reproducible interferon-inflammatory-complement-myeloid program. This signal distinguishes active disease or progression from LTBI/exposed non-progression, but it likely reflects both immune activation and changes in circulating leukocyte composition.

## Important limitations

- Whole-blood RNA-seq cannot distinguish cell-intrinsic transcriptional regulation from changes in cell-type abundance without deconvolution or single-cell/sorted-cell validation.
- GSE79362 contains repeated longitudinal samples, and adjustment for repeated measures substantially reduces the strict DEG count.
- GSE94438 spans three African sites, so site/country must remain a priority covariate in final models.
- The label definitions differ between datasets, so AUC and DEG results should not be interpreted as identical biological contrasts.
