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

## Gene Set Databases Used in GSEA

| Role | Gene Set Database | Notes |
|------|-------------------|-------|
| **Primary database** | **MSigDB Hallmark gene sets** | Main results; used for cross-dataset pathway-level comparison. Only 50 gene sets with low redundancy — well suited for summarising the core TB immune programmes: IFN, inflammation, complement, TNF/NF-κB, and IL6/JAK/STAT3. |
| **Supporting database 1** | **Reactome pathways** | Auxiliary pathway-level validation; confirms whether the immune/inflammatory themes seen in Hallmark are reproducible in a more granular pathway database. |
| **Supporting database 2** | **GO Biological Process (GO BP)** | Auxiliary biological-process-level validation; checks whether broader immune and inflammatory biological processes are consistent across datasets. |

## Dataset Comparison: GSE79362 vs GSE94438

| Dimension | GSE79362 | GSE94438 | Notes |
|-----------|----------|----------|-------|
| Dataset role | Active TB vs latent TB infection | TB patients / progressors vs exposed household contacts / controls | Same topic, different control group definitions |
| Core biological question | In an already-infected background, what signals distinguish **latent TB** from **active pulmonary TB**? | In an exposed population, what signals associate with **TB progression / active disease**? | GSE79362 = disease-state comparison; GSE94438 = progression/exposure cohort |
| Primary contrast | **PTB vs LTBI** | **PTB / progressors vs household contacts / controls** | Cannot treat both as simply "patients vs healthy" |
| Geographic origin | South Africa | Ethiopia, South Africa, The Gambia | GSE94438 has greater geographic heterogeneity |
| Tissue | Whole blood | Whole blood | Identical — key reason cross-dataset comparison is valid |
| Sequencing platform | Illumina RNA-seq | Illumina RNA-seq | Consistent platform — another reason for comparability |
| Data type | Bulk whole-blood RNA-seq | Bulk whole-blood RNA-seq | Neither is single-cell; cell-type attribution requires deconvolution |
| Sample size | 355 labelled samples | 428 labelled samples (434 total; 6 excluded for missing TBStatus) | — |
| Control group nature | LTBI: infected but not progressed to active TB | Household contacts: TB-exposed but not necessarily truly naive/uninfected | GSE94438 controls are not strict "healthy naive" individuals |
| Cohort structure | Single country, relatively homogeneous | Multi-country, multi-site, higher heterogeneity | Site/country effect must be handled carefully in GSE94438 |
| Repeated measures | Marked repeated longitudinal measurements | Some longitudinal structure, but main confounder is site/country | GSE79362 most needs PatientID blocking |
| Main confounders | Timepoint, PatientID repeated measures | Site/country, sex, age | The two datasets cannot share an identical adjusted model |
| Adjustment method | limma-voom + duplicateCorrelation for PatientID block / repeated measures | edgeR model adjusted for site/country, sex, age | — |
| Baseline DEGs | 30 strict DEGs | 43 strict DEGs | — |
| Adjusted DEGs | 9 strict DEGs | 41 strict DEGs | GSE79362 drops sharply; GSE94438 retains most signal after adjustment |
| Effect of adjustment | 30 → 9 after adjustment | 43 → 41 after adjustment | GSE79362 is highly sensitive to repeated measures; GSE94438 disease signal is not driven by site/sex/age |
| DEG direction | Baseline: 29 up, 1 down; fewer after adjustment | Baseline: 43 up; adjusted: 41 up | Both datasets are dominated by upregulated signal in active/progressive TB |
| Single-gene stability | Strongly affected by repeated measures | Affected by geographic heterogeneity, but largely retained after adjustment | Individual DEGs are not the most stable layer |
| GSEA / pathway results | IFN, inflammation, complement, TNF/NF-κB, IL6/JAK/STAT3 enriched | Same: IFN, inflammation, complement, TNF/NF-κB, IL6/JAK/STAT3 enriched | Pathway-level results are more consistent than DEG lists |
| Pathway consistency | Highly consistent with GSE94438 | Highly consistent with GSE79362 | Core rationale for analysing both datasets together |
| Signature AUC | Zak16, RISK4, Eleven_gene: ~0.76–0.77 | Zak16, RISK4, Eleven_gene: ~0.68–0.69 | Signatures perform better in GSE79362 |
| Why higher AUC in GSE79362 | PTB vs LTBI is a more direct biological contrast | Exposed household contacts may already have TB-related immune priming | Lower AUC in GSE94438 is biologically expected |
| Best suited for | Disease-specific whole-blood signal: active TB vs latent infection | TB risk / disease-associated signal in an exposure/progression context | Complementary, not interchangeable |
| Main strength | Clean contrast, single-country cohort, direct PTB vs LTBI biological interpretation | Larger sample size, multi-country, closer to real-world exposure/progression setting | One is cleaner; the other is more realistic |
| Main limitation | Repeated measures have large impact; ignoring them inflates DEG counts | Strong site/country/population stratification; controls are not naive uninfected individuals | Both have limitations that must be stated |
| Reporting role | Primary analysis: "active disease vs latent infection" | Cross-cohort validation or secondary analysis: "progression-exposure comparison" | GSE79362 as primary; GSE94438 as validation/complement |
| Interpretation risk | Without PatientID / timepoint adjustment, DEGs may be inflated | Without site/country adjustment, geographic differences may be mistaken for TB biology | Neither dataset should be run with a simple `~ TBStatus` model only |
| Biological conclusion | Active PTB vs LTBI shows a stronger interferon–inflammatory–myeloid signal | Progressive/active TB vs exposed controls shows a similar immune programme | Together they support a conserved whole-blood TB immune programme |
| Value for downstream analysis | Defines active TB disease-associated markers | Tests whether those markers replicate in an exposure/progression cohort | Both together are more convincing than either alone |
