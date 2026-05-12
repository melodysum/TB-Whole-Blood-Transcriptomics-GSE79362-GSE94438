suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(pheatmap)
  library(matrixStats)
  library(pROC)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

outdir <- "/private/tmp/tb_curated_results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))

load_dataset <- function(study) {
  obj <- readRDS(file.path("/private/tmp", paste0(study, "_curatedTBData.rds")))[[study]]
  counts <- as.matrix(experiments(obj)[["assay_reprocess_hg38"]])
  meta <- as.data.frame(colData(obj))
  meta$sample_id <- rownames(meta)
  common <- intersect(colnames(counts), meta$sample_id)
  counts <- counts[, common, drop = FALSE]
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]
  keep <- !is.na(meta$TBStatus)
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  meta$TBStatus <- factor(meta$TBStatus)
  meta$Gender <- factor(meta$Gender)
  meta$GeographicalRegion <- factor(meta$GeographicalRegion)
  list(study = study, counts = counts, meta = meta, raw_n = ncol(experiments(obj)[["assay_reprocess_hg38"]]))
}

collapse_duplicate_genes <- function(counts) {
  rowsum(counts, group = rownames(counts), reorder = FALSE)
}

filter_counts <- function(counts, group) {
  counts <- collapse_duplicate_genes(counts)
  min_n <- min(table(group))
  keep <- rowSums(cpm(counts) > 1) >= min_n
  list(counts = counts[keep, , drop = FALSE], keep = keep, min_group_n = min_n)
}

audit_dataset <- function(x, geo_total) {
  m <- x$meta
  audit <- list(
    study = x$study,
    geo_total = geo_total,
    curated_hg38_total = x$raw_n,
    retained_for_labelled_analysis = nrow(m),
    tbstatus = as.data.frame(table(m$TBStatus)) %>% setNames(c("TBStatus", "n")),
    site_status = as.data.frame(table(m$GeographicalRegion, m$TBStatus)) %>% setNames(c("site", "TBStatus", "n")),
    unique_patients = length(unique(m$PatientID)),
    repeated_patients = sum(table(m$PatientID) > 1),
    max_samples_per_patient = max(table(m$PatientID)),
    measurement_time = as.data.frame(table(m$MeasurementTime, useNA = "ifany")) %>% setNames(c("MeasurementTime", "n")),
    progression_status = as.data.frame(table(m$Progression, m$TBStatus, useNA = "ifany")) %>% setNames(c("Progression", "TBStatus", "n"))
  )
  write.csv(audit$tbstatus, file.path(outdir, paste0(x$study, "_TBStatus_counts.csv")), row.names = FALSE)
  write.csv(audit$site_status, file.path(outdir, paste0(x$study, "_site_by_TBStatus.csv")), row.names = FALSE)
  write.csv(as.data.frame(colnames(m)), file.path(outdir, paste0(x$study, "_metadata_columns.csv")), row.names = FALSE)
  audit
}

plot_qc <- function(x, filt) {
  counts <- collapse_duplicate_genes(x$counts)
  meta <- x$meta
  meta$library_size <- colSums(counts)[meta$sample_id]
  meta$low_lib <- meta$library_size < 1e6
  write.csv(meta %>% dplyr::select(sample_id, TBStatus, GeographicalRegion, PatientID, Age, Gender, library_size, low_lib),
            file.path(outdir, paste0(x$study, "_sample_qc.csv")), row.names = FALSE)

  p_lib <- ggplot(meta, aes(library_size, fill = TBStatus)) +
    geom_histogram(bins = 35, alpha = 0.75, colour = "white") +
    geom_vline(xintercept = 1e6, colour = "red", linewidth = 0.7) +
    scale_x_continuous(labels = scales::comma) +
    labs(title = paste0(x$study, ": library size"), subtitle = "Red line flags <1e6 reads", x = "Library size", y = "Samples") +
    theme(legend.position = "top")
  ggsave(file.path(figdir, paste0(x$study, "_library_size.png")), p_lib, width = 7, height = 4, dpi = 220)

  y <- DGEList(filt$counts)
  y <- calcNormFactors(y)
  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  pca <- prcomp(t(logcpm), scale. = TRUE)
  var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pc <- as.data.frame(pca$x[, 1:4]) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta, by = "sample_id")
  outlier <- abs(scale(pc$PC1)) > 2.5 | abs(scale(pc$PC2)) > 2.5

  p_status <- ggplot(pc, aes(PC1, PC2, colour = TBStatus)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(data = pc[outlier, ], aes(label = sample_id), size = 2.7, max.overlaps = 20) +
    labs(title = paste0(x$study, ": PCA by TBStatus"), x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
    theme(legend.position = "top")
  ggsave(file.path(figdir, paste0(x$study, "_PCA_TBStatus.png")), p_status, width = 6, height = 5, dpi = 220)

  p_site <- ggplot(pc, aes(PC1, PC2, colour = GeographicalRegion)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(data = pc[outlier, ], aes(label = sample_id), size = 2.7, max.overlaps = 20) +
    labs(title = paste0(x$study, ": PCA by site/country"), x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
    theme(legend.position = "top")
  ggsave(file.path(figdir, paste0(x$study, "_PCA_site.png")), p_site, width = 6, height = 5, dpi = 220)

  vars <- rowVars(logcpm)
  top <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
  png(file.path(figdir, paste0(x$study, "_sample_correlation_heatmap.png")), width = 1600, height = 1400, res = 180)
  pheatmap(cor(logcpm[top, ]), show_colnames = FALSE, show_rownames = FALSE,
           annotation_col = meta %>% dplyr::select(TBStatus, GeographicalRegion) %>% as.data.frame() %>% `rownames<-`(meta$sample_id),
           main = paste0(x$study, ": sample correlation, top 500 variable genes"))
  dev.off()

  list(logcpm = logcpm, pca = pc, pca_var = var)
}

run_baseline <- function(x, filt) {
  meta <- x$meta
  counts <- filt$counts
  # Reference = Control/LTBI, coefficient = PTB higher vs reference.
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$TBStatus <- relevel(factor(meta$TBStatus), ref = ref)
  y <- DGEList(counts = counts, samples = meta)
  y <- calcNormFactors(y)
  design <- model.matrix(~ TBStatus, data = meta)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = grep("^TBStatus", colnames(design)))
  res <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(FDR = p.adjust(PValue, "BH"), sig = FDR < 0.05 & abs(logFC) > 1)
  write.csv(res, file.path(outdir, paste0(x$study, "_baseline_edgeR_DEG.csv")), row.names = FALSE)

  lab <- res %>% arrange(PValue) %>% slice_head(n = 20)
  p_volc <- ggplot(res, aes(logFC, -log10(FDR), colour = sig)) +
    geom_point(alpha = 0.55, size = 1) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.8, max.overlaps = 30) +
    scale_colour_manual(values = c("grey65", "#B2182B")) +
    labs(title = paste0(x$study, ": baseline DE volcano"), x = "logFC: PTB/progressor vs reference", y = "-log10(FDR)") +
    theme(legend.position = "none")
  ggsave(file.path(figdir, paste0(x$study, "_volcano.png")), p_volc, width = 6.5, height = 5.2, dpi = 220)

  p_ma <- ggplot(res, aes(logCPM, logFC, colour = sig)) +
    geom_point(alpha = 0.55, size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_colour_manual(values = c("grey65", "#2166AC")) +
    labs(title = paste0(x$study, ": MA plot"), x = "logCPM", y = "logFC") +
    theme(legend.position = "none")
  ggsave(file.path(figdir, paste0(x$study, "_MA.png")), p_ma, width = 6.5, height = 5.2, dpi = 220)

  res
}

run_adjusted <- function(x, filt, baseline) {
  meta <- x$meta
  counts <- filt$counts
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$TBStatus <- relevel(factor(meta$TBStatus), ref = ref)
  y <- DGEList(counts = counts, samples = meta)
  y <- calcNormFactors(y)

  if (x$study == "GSE94438") {
    design <- model.matrix(~ GeographicalRegion + Gender + Age + TBStatus, data = meta)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef = grep("^TBStatus", colnames(design)))
    method <- "edgeR adjusted for site + sex + age"
    res <- topTags(test, n = Inf)$table %>% rownames_to_column("gene") %>% mutate(FDR = p.adjust(PValue, "BH"))
  } else {
    design <- model.matrix(~ MeasurementTime + TBStatus, data = meta)
    v <- voom(y, design, plot = FALSE)
    dup <- duplicateCorrelation(v, design, block = meta$PatientID)
    fit <- lmFit(v, design, block = meta$PatientID, correlation = dup$consensus)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = grep("^TBStatus", colnames(design)), number = Inf, sort.by = "none") %>%
      rownames_to_column("gene") %>%
      rename(PValue = P.Value, FDR = adj.P.Val, logCPM = AveExpr)
    method <- paste0("limma-voom adjusted for timepoint with PatientID block, correlation=", round(dup$consensus, 3))
  }

  res <- res %>% mutate(sig = FDR < 0.05 & abs(logFC) > 1)
  write.csv(res, file.path(outdir, paste0(x$study, "_adjusted_DEG.csv")), row.names = FALSE)
  comp <- baseline %>% dplyr::select(gene, logFC_baseline = logFC, FDR_baseline = FDR) %>%
    inner_join(res %>% dplyr::select(gene, logFC_adjusted = logFC, FDR_adjusted = FDR), by = "gene") %>%
    mutate(sig_baseline = FDR_baseline < 0.05 & abs(logFC_baseline) > 1,
           sig_adjusted = FDR_adjusted < 0.05 & abs(logFC_adjusted) > 1,
           change = case_when(sig_baseline & !sig_adjusted ~ "Lost",
                              !sig_baseline & sig_adjusted ~ "Gained",
                              sig_baseline & sig_adjusted ~ "Stable significant",
                              TRUE ~ "Not significant"))
  write.csv(comp, file.path(outdir, paste0(x$study, "_baseline_vs_adjusted.csv")), row.names = FALSE)
  p <- ggplot(comp, aes(logFC_baseline, logFC_adjusted, colour = change)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(title = paste0(x$study, ": baseline vs adjusted logFC"), subtitle = method,
         x = "Baseline logFC", y = "Adjusted logFC") +
    theme(legend.position = "top")
  ggsave(file.path(figdir, paste0(x$study, "_baseline_vs_adjusted.png")), p, width = 6.5, height = 5.2, dpi = 220)
  attr(res, "method") <- method
  res
}

run_signature <- function(x, logcpm, genes, name) {
  meta <- x$meta
  present <- intersect(genes, rownames(logcpm))
  missing <- setdiff(genes, rownames(logcpm))
  score <- colMeans(logcpm[present, , drop = FALSE])
  meta$score <- score[meta$sample_id]
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$outcome <- factor(ifelse(meta$TBStatus == ref, ref, "PTB"), levels = c(ref, "PTB"))
  rocobj <- roc(meta$outcome, meta$score, levels = c(ref, "PTB"), direction = "<", quiet = TRUE)
  coords_best <- as.data.frame(t(coords(rocobj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")))
  p <- ggplot(meta, aes(TBStatus, score, fill = TBStatus)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1) +
    labs(title = paste0(x$study, ": ", name, " score"), y = "Mean log2 CPM score", x = NULL) +
    theme(legend.position = "none")
  ggsave(file.path(figdir, paste0(x$study, "_", name, "_boxplot.png")), p, width = 5.5, height = 4.2, dpi = 220)

  if (length(levels(meta$GeographicalRegion)) > 1) {
    ps <- p + facet_wrap(~ GeographicalRegion) + labs(title = paste0(x$study, ": ", name, " score by site"))
    ggsave(file.path(figdir, paste0(x$study, "_", name, "_boxplot_site.png")), ps, width = 8, height = 4.6, dpi = 220)
  }

  tibble(study = x$study, Signature = name, genes_total = length(genes),
         genes_present = length(present), missing_genes = paste(missing, collapse = ";"),
         AUC = as.numeric(auc(rocobj)),
         Sensitivity = coords_best$sensitivity,
         Specificity = coords_best$specificity,
         WilcoxP = wilcox.test(score ~ meta$outcome)$p.value)
}

risk4_ratio <- function(x, logcpm) {
  genes <- c("GBP5", "SEPTIN4", "CDO1", "TRAV27")
  meta <- x$meta
  if (!all(genes %in% rownames(logcpm))) return(NULL)
  meta$score <- (logcpm["GBP5", ] + logcpm["SEPTIN4", ]) / (logcpm["CDO1", ] + logcpm["TRAV27", ])
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$outcome <- factor(ifelse(meta$TBStatus == ref, ref, "PTB"), levels = c(ref, "PTB"))
  rocobj <- roc(meta$outcome, meta$score, levels = c(ref, "PTB"), direction = "<", quiet = TRUE)
  coords_best <- as.data.frame(t(coords(rocobj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")))
  p <- ggplot(meta, aes(TBStatus, score, fill = TBStatus)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1) +
    labs(title = paste0(x$study, ": RISK4 ratio"), y = "(GBP5 + SEPTIN4) / (CDO1 + TRAV27)", x = NULL) +
    theme(legend.position = "none")
  ggsave(file.path(figdir, paste0(x$study, "_RISK4_ratio_boxplot.png")), p, width = 5.5, height = 4.2, dpi = 220)
  tibble(study = x$study, Signature = "RISK4 ratio", genes_total = 4, genes_present = 4,
         missing_genes = "", AUC = as.numeric(auc(rocobj)),
         Sensitivity = coords_best$sensitivity,
         Specificity = coords_best$specificity,
         WilcoxP = wilcox.test(meta$score ~ meta$outcome)$p.value)
}

simple_gsea_hallmark <- function(deg, study) {
  ranks <- deg %>%
    filter(!is.na(PValue), PValue > 0, !is.na(logFC)) %>%
    mutate(score = sign(logFC) * -log10(PValue)) %>%
    group_by(gene) %>% slice_max(abs(score), n = 1, with_ties = FALSE) %>% ungroup()
  gl <- ranks$score
  names(gl) <- ranks$gene
  gl <- sort(gl, decreasing = TRUE)
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
  g <- GSEA(gl, TERM2GENE = hallmark, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE)
  df <- as.data.frame(g) %>% transmute(study = study, pathway = ID, NES, pvalue, p.adjust, leading_edge_genes = core_enrichment)
  write.csv(df, file.path(outdir, paste0(study, "_GSEA_Hallmark.csv")), row.names = FALSE)
  top <- df %>% arrange(p.adjust) %>% slice_head(n = 20)
  p <- ggplot(top, aes(reorder(pathway, NES), NES, colour = p.adjust)) +
    geom_point(size = 3) +
    coord_flip() +
    scale_colour_viridis_c(direction = -1) +
    labs(title = paste0(study, ": Hallmark GSEA top pathways"), x = NULL, y = "NES")
  ggsave(file.path(figdir, paste0(study, "_GSEA_Hallmark_top20.png")), p, width = 8, height = 5.5, dpi = 220)
  df
}

compare_de <- function(deg1, deg2) {
  shared <- inner_join(deg1, deg2, by = "gene", suffix = c("_79362", "_94438"))
  shared <- shared %>%
    mutate(sig_79362 = FDR_79362 < 0.05 & abs(logFC_79362) > 1,
           sig_94438 = FDR_94438 < 0.05 & abs(logFC_94438) > 1,
           sig_status = case_when(sig_79362 & sig_94438 ~ "Both",
                                  sig_79362 ~ "GSE79362 only",
                                  sig_94438 ~ "GSE94438 only",
                                  TRUE ~ "Neither"))
  p <- ggplot(shared, aes(logFC_79362, logFC_94438, colour = sig_status)) +
    geom_point(alpha = 0.55, size = 1) +
    geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(title = "Cross-dataset DE consistency", x = "GSE79362 logFC", y = "GSE94438 logFC") +
    theme(legend.position = "top")
  ggsave(file.path(figdir, "cross_dataset_logFC_scatter.png"), p, width = 6.5, height = 5.2, dpi = 220)
  write.csv(shared, file.path(outdir, "cross_dataset_DEG_comparison.csv"), row.names = FALSE)
  shared
}

compare_gsea <- function(g1, g2) {
  shared <- inner_join(g1, g2, by = "pathway", suffix = c("_79362", "_94438")) %>%
    mutate(sig = case_when(p.adjust_79362 < 0.25 & p.adjust_94438 < 0.25 ~ "Both FDR<0.25",
                           p.adjust_79362 < 0.25 | p.adjust_94438 < 0.25 ~ "One dataset",
                           TRUE ~ "Neither"))
  key <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING")
  p <- ggplot(shared, aes(NES_79362, NES_94438, colour = sig)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(data = shared %>% filter(pathway %in% key), aes(label = pathway), size = 2.7, max.overlaps = 20) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(title = "Hallmark pathway consistency", x = "GSE79362 NES", y = "GSE94438 NES") +
    theme(legend.position = "top")
  ggsave(file.path(figdir, "cross_dataset_Hallmark_NES_scatter.png"), p, width = 7, height = 5.5, dpi = 220)
  write.csv(shared, file.path(outdir, "cross_dataset_Hallmark_GSEA_comparison.csv"), row.names = FALSE)
  shared
}

studies <- list(GSE79362 = load_dataset("GSE79362"), GSE94438 = load_dataset("GSE94438"))
audits <- list(GSE79362 = audit_dataset(studies$GSE79362, 355), GSE94438 = audit_dataset(studies$GSE94438, 434))

filtered <- lapply(studies, function(x) filter_counts(x$counts, x$meta$TBStatus))
qc <- Map(plot_qc, studies, filtered)
baseline <- Map(run_baseline, studies, filtered)
adjusted <- Map(run_adjusted, studies, filtered, baseline)
gsea <- Map(simple_gsea_hallmark, baseline, names(baseline))
decomp <- compare_de(baseline$GSE79362, baseline$GSE94438)
gseacomp <- compare_gsea(gsea$GSE79362, gsea$GSE94438)

genes_zak16 <- c("GBP5", "BATF2", "FCGR1B", "SCARF1", "TRAV27", "ISG15", "ANKRD22",
                 "ETV7", "SERPING1", "SAMD9L", "IFIT2", "IFIT3", "IFI44L", "CXCL10", "HERC5", "OAS1")
genes_risk4 <- c("GBP5", "SEPTIN4", "CDO1", "TRAV27")
genes_11 <- c("GBP5", "BATF2", "FCGR1B", "ANKRD22", "ETV7", "SERPING1", "SAMD9L",
              "IFI44L", "CXCL10", "HERC5", "OAS1")

sigsum <- bind_rows(
  run_signature(studies$GSE79362, qc$GSE79362$logcpm, genes_zak16, "Zak16"),
  run_signature(studies$GSE94438, qc$GSE94438$logcpm, genes_zak16, "Zak16"),
  run_signature(studies$GSE79362, qc$GSE79362$logcpm, genes_risk4, "RISK4_mean"),
  run_signature(studies$GSE94438, qc$GSE94438$logcpm, genes_risk4, "RISK4_mean"),
  risk4_ratio(studies$GSE79362, qc$GSE79362$logcpm),
  risk4_ratio(studies$GSE94438, qc$GSE94438$logcpm),
  run_signature(studies$GSE79362, qc$GSE79362$logcpm, genes_11, "Eleven_gene"),
  run_signature(studies$GSE94438, qc$GSE94438$logcpm, genes_11, "Eleven_gene")
)
write.csv(sigsum, file.path(outdir, "signature_AUC_summary.csv"), row.names = FALSE)

summary_rows <- bind_rows(lapply(names(studies), function(s) {
  m <- studies[[s]]$meta
  f <- filtered[[s]]
  b <- baseline[[s]]
  a <- adjusted[[s]]
  tibble(
    study = s,
    geo_total = ifelse(s == "GSE79362", 355, 434),
    curated_hg38_total = studies[[s]]$raw_n,
    retained_labelled = nrow(m),
    groups = paste(capture.output(print(table(m$TBStatus))), collapse = " "),
    sites = paste(capture.output(print(table(m$GeographicalRegion, useNA = "ifany"))), collapse = " "),
    unique_patients = length(unique(m$PatientID)),
    repeated_patients = sum(table(m$PatientID) > 1),
    genes_before = nrow(collapse_duplicate_genes(studies[[s]]$counts)),
    genes_after_filter = nrow(f$counts),
    baseline_deg = sum(b$FDR < 0.05 & abs(b$logFC) > 1),
    baseline_up = sum(b$FDR < 0.05 & b$logFC > 1),
    baseline_down = sum(b$FDR < 0.05 & b$logFC < -1),
    adjusted_deg = sum(a$FDR < 0.05 & abs(a$logFC) > 1),
    adjusted_method = attr(a, "method")
  )
}))
write.csv(summary_rows, file.path(outdir, "analysis_summary.csv"), row.names = FALSE)

cross_summary <- tibble(
  metric = c("shared_genes_tested", "spearman_logFC", "shared_sig_DEGs",
             "shared_up_DEGs", "shared_down_DEGs", "top100_overlap",
             "hallmark_NES_spearman"),
  value = c(nrow(decomp),
            cor(decomp$logFC_79362, decomp$logFC_94438, method = "spearman"),
            sum(decomp$sig_79362 & decomp$sig_94438),
            length(intersect(baseline$GSE79362 %>% filter(FDR < 0.05, logFC > 1) %>% pull(gene),
                             baseline$GSE94438 %>% filter(FDR < 0.05, logFC > 1) %>% pull(gene))),
            length(intersect(baseline$GSE79362 %>% filter(FDR < 0.05, logFC < -1) %>% pull(gene),
                             baseline$GSE94438 %>% filter(FDR < 0.05, logFC < -1) %>% pull(gene))),
            length(intersect(baseline$GSE79362 %>% arrange(PValue) %>% slice_head(n = 100) %>% pull(gene),
                             baseline$GSE94438 %>% arrange(PValue) %>% slice_head(n = 100) %>% pull(gene))),
            cor(gseacomp$NES_79362, gseacomp$NES_94438, method = "spearman"))
)
write.csv(cross_summary, file.path(outdir, "cross_dataset_summary.csv"), row.names = FALSE)

sink(file.path(outdir, "plain_language_interpretation.txt"))
cat("TB curated analysis interpretation\n\n")
print(summary_rows)
cat("\nSignature AUC summary:\n")
print(sigsum)
cat("\nCross-dataset summary:\n")
print(cross_summary)
cat("\nKey interpretation:\n")
cat("GSE79362 compares PTB vs LTBI in a South African cohort with repeated longitudinal samples; it is best interpreted as disease-associated signal on an infected/latent background.\n")
cat("GSE94438 compares progressors/PTB vs household-contact controls across Ethiopia, South Africa, and The Gambia; it is best interpreted as progression-risk/disease-emergence signal among exposed contacts, with site as a major confounder.\n")
cat("Shared whole-blood signals should be interpreted as TB disease/progression biology: interferon response, inflammatory/myeloid activation, complement and innate immune pathways. These may partly reflect neutrophil/monocyte abundance shifts rather than purely within-cell regulation.\n")
sink()

message("Done. Outputs in: ", outdir)
