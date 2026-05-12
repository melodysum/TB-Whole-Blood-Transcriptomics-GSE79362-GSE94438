suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(pROC)
  library(msigdbr)
  library(clusterProfiler)
})

outdir <- "/private/tmp/tb_curated_results"
figdir <- file.path(outdir, "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

load_x <- function(study) {
  obj <- readRDS(file.path("/private/tmp", paste0(study, "_curatedTBData.rds")))[[study]]
  counts <- as.matrix(experiments(obj)[["assay_reprocess_hg38"]])
  meta <- as.data.frame(colData(obj)); meta$sample_id <- rownames(meta)
  common <- intersect(colnames(counts), meta$sample_id)
  counts <- counts[, common]; meta <- meta[match(common, meta$sample_id),]
  keep <- !is.na(meta$TBStatus)
  counts <- counts[, keep]; meta <- meta[keep,]
  rows <- rowSums(counts)
  counts <- rowsum(counts[rows > 0,], rownames(counts)[rows > 0], reorder = FALSE)
  meta$TBStatus <- factor(meta$TBStatus)
  meta$Gender <- factor(meta$Gender)
  meta$GeographicalRegion <- factor(meta$GeographicalRegion)
  list(study = study, counts = counts, meta = meta, raw_n = ncol(experiments(obj)[["assay_reprocess_hg38"]]))
}

filter_x <- function(x) {
  min_n <- min(table(x$meta$TBStatus))
  keep <- rowSums(cpm(x$counts) > 1) >= min_n
  x$counts[keep,]
}

logcpm_x <- function(counts) {
  y <- DGEList(counts); y <- calcNormFactors(y); cpm(y, log = TRUE, prior.count = 1)
}

run_adj_fast <- function(x, counts, base_file) {
  meta <- x$meta
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$TBStatus <- relevel(factor(meta$TBStatus), ref)
  y <- DGEList(counts, samples = meta); y <- calcNormFactors(y)
  if (x$study == "GSE94438") {
    design <- model.matrix(~ GeographicalRegion + Gender + Age + TBStatus, meta)
    method <- "edgeR adjusted for country/site, sex, age"
  } else {
    design <- model.matrix(~ MeasurementTime + Gender + Age + TBStatus, meta)
    method <- "edgeR adjusted for timepoint, sex, age; repeated PatientID not modelled in this fast run"
  }
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = grep("^TBStatus", colnames(design)))
  res <- topTags(qlf, n = Inf)$table %>% rownames_to_column("gene") %>% mutate(FDR = p.adjust(PValue, "BH"))
  write.csv(res, file.path(outdir, paste0(x$study, "_adjusted_fast_DEG.csv")), row.names = FALSE)
  base <- read.csv(base_file)
  comp <- base %>% select(gene, logFC_baseline = logFC, FDR_baseline = FDR) %>%
    inner_join(res %>% select(gene, logFC_adjusted = logFC, FDR_adjusted = FDR), by = "gene") %>%
    mutate(sig_baseline = FDR_baseline < 0.05 & abs(logFC_baseline) > 1,
           sig_adjusted = FDR_adjusted < 0.05 & abs(logFC_adjusted) > 1)
  write.csv(comp, file.path(outdir, paste0(x$study, "_baseline_vs_adjusted_fast.csv")), row.names = FALSE)
  p <- ggplot(comp, aes(logFC_baseline, logFC_adjusted, colour = sig_baseline != sig_adjusted)) +
    geom_point(alpha = 0.55, size = 1) + geom_abline(linetype = 2) +
    labs(title = paste0(x$study, ": baseline vs adjusted"), subtitle = method,
         x = "Baseline logFC", y = "Adjusted logFC") + theme_bw()
  ggsave(file.path(figdir, paste0(x$study, "_baseline_vs_adjusted_fast.png")), p, width = 6.5, height = 5, dpi = 220)
  attr(res, "method") <- method
  res
}

sig_score <- function(x, logcpm, genes, name) {
  present <- intersect(genes, rownames(logcpm))
  meta <- x$meta
  meta$score <- colMeans(logcpm[present,,drop=FALSE])
  ref <- if ("Control" %in% levels(meta$TBStatus)) "Control" else "LTBI"
  meta$outcome <- factor(ifelse(meta$TBStatus == ref, ref, "PTB"), levels = c(ref, "PTB"))
  rocobj <- roc(meta$outcome, meta$score, levels = c(ref, "PTB"), direction = "<", quiet = TRUE)
  cc <- as.data.frame(t(coords(rocobj, "best", ret = c("threshold","sensitivity","specificity"), best.method = "youden")))
  p <- ggplot(meta, aes(TBStatus, score, fill = TBStatus)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width=.15, alpha=.45, size=.9) +
    labs(title = paste0(x$study, ": ", name), y = "Mean log2 CPM score", x = NULL) +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(figdir, paste0(x$study, "_", name, "_score.png")), p, width = 5.5, height = 4.2, dpi = 220)
  if (length(levels(meta$GeographicalRegion)) > 1) {
    ggsave(file.path(figdir, paste0(x$study, "_", name, "_score_by_site.png")),
           p + facet_wrap(~GeographicalRegion), width = 8, height = 4.5, dpi = 220)
  }
  tibble(study=x$study, Signature=name, genes_present=length(present),
         missing_genes=paste(setdiff(genes, rownames(logcpm)), collapse=";"),
         AUC=as.numeric(auc(rocobj)), Sensitivity=cc$sensitivity,
         Specificity=cc$specificity, WilcoxP=wilcox.test(meta$score ~ meta$outcome)$p.value)
}

gsea_h <- function(study) {
  deg <- read.csv(file.path(outdir, paste0(study, "_baseline_edgeR_DEG.csv")))
  ranks <- deg %>% filter(PValue > 0) %>% mutate(score = sign(logFC) * -log10(PValue)) %>%
    group_by(gene) %>% slice_max(abs(score), n=1, with_ties=FALSE) %>% ungroup()
  gl <- ranks$score; names(gl) <- ranks$gene; gl <- sort(gl, decreasing=TRUE)
  h <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, gene_symbol)
  gs <- GSEA(gl, TERM2GENE=h, pvalueCutoff=1, verbose=FALSE)
  df <- as.data.frame(gs) %>% transmute(study=study, pathway=ID, NES, pvalue, p.adjust, leading_edge_genes=core_enrichment)
  write.csv(df, file.path(outdir, paste0(study, "_GSEA_Hallmark.csv")), row.names=FALSE)
  top <- df %>% arrange(p.adjust) %>% slice_head(n=18)
  p <- ggplot(top, aes(reorder(pathway, NES), NES, colour=p.adjust)) +
    geom_point(size=3) + coord_flip() + theme_bw() + scale_colour_viridis_c(direction=-1) +
    labs(title=paste0(study, ": Hallmark GSEA"), x=NULL, y="NES")
  ggsave(file.path(figdir, paste0(study, "_GSEA_Hallmark.png")), p, width=8, height=5.2, dpi=220)
  df
}

x1 <- load_x("GSE79362"); x2 <- load_x("GSE94438")
c1 <- filter_x(x1); c2 <- filter_x(x2)
l1 <- logcpm_x(c1); l2 <- logcpm_x(c2)
a1 <- run_adj_fast(x1, c1, file.path(outdir, "GSE79362_baseline_edgeR_DEG.csv"))
a2 <- run_adj_fast(x2, c2, file.path(outdir, "GSE94438_baseline_edgeR_DEG.csv"))

zak16 <- c("GBP5","BATF2","FCGR1B","SCARF1","TRAV27","ISG15","ANKRD22","ETV7","SERPING1","SAMD9L","IFIT2","IFIT3","IFI44L","CXCL10","HERC5","OAS1")
risk4 <- c("GBP5","SEPTIN4","CDO1","TRAV27")
gene11 <- c("GBP5","BATF2","FCGR1B","ANKRD22","ETV7","SERPING1","SAMD9L","IFI44L","CXCL10","HERC5","OAS1")
sigsum <- bind_rows(sig_score(x1,l1,zak16,"Zak16"), sig_score(x2,l2,zak16,"Zak16"),
                    sig_score(x1,l1,risk4,"RISK4"), sig_score(x2,l2,risk4,"RISK4"),
                    sig_score(x1,l1,gene11,"Eleven_gene"), sig_score(x2,l2,gene11,"Eleven_gene"))
write.csv(sigsum, file.path(outdir, "signature_AUC_summary.csv"), row.names=FALSE)

g1 <- gsea_h("GSE79362"); g2 <- gsea_h("GSE94438")
gs <- inner_join(g1,g2,by="pathway",suffix=c("_79362","_94438"))
key <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_COMPLEMENT")
p <- ggplot(gs, aes(NES_79362,NES_94438,colour=(p.adjust_79362<.25 & p.adjust_94438<.25))) +
  geom_point(size=2.2) + geom_text_repel(data=gs %>% filter(pathway %in% key), aes(label=pathway), size=2.7) +
  geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=2)+theme_bw()+
  labs(title="Hallmark pathway NES consistency", x="GSE79362 NES", y="GSE94438 NES", colour="Both FDR<0.25")
ggsave(file.path(figdir, "cross_dataset_Hallmark_NES_scatter.png"), p, width=7, height=5.5, dpi=220)
write.csv(gs, file.path(outdir, "cross_dataset_Hallmark_GSEA_comparison.csv"), row.names=FALSE)

d1 <- read.csv(file.path(outdir, "GSE79362_baseline_edgeR_DEG.csv"))
d2 <- read.csv(file.path(outdir, "GSE94438_baseline_edgeR_DEG.csv"))
dc <- inner_join(d1,d2,by="gene",suffix=c("_79362","_94438")) %>%
  mutate(sig_79362=FDR_79362<.05 & abs(logFC_79362)>1,
         sig_94438=FDR_94438<.05 & abs(logFC_94438)>1,
         status=case_when(sig_79362 & sig_94438 ~ "Both", sig_79362 ~ "GSE79362 only", sig_94438 ~ "GSE94438 only", TRUE ~ "Neither"))
write.csv(dc, file.path(outdir, "cross_dataset_DEG_comparison.csv"), row.names=FALSE)
p <- ggplot(dc, aes(logFC_79362,logFC_94438,colour=status)) + geom_point(alpha=.55,size=1) +
  geom_smooth(method="loess",se=FALSE,colour="black") + theme_bw() +
  labs(title="Cross-dataset gene-level logFC comparison", x="GSE79362 logFC", y="GSE94438 logFC")
ggsave(file.path(figdir, "cross_dataset_logFC_scatter.png"), p, width=6.5, height=5, dpi=220)

summary <- bind_rows(lapply(list(x1,x2), function(x) {
  deg <- read.csv(file.path(outdir, paste0(x$study, "_baseline_edgeR_DEG.csv")))
  adj <- if (x$study=="GSE79362") a1 else a2
  tibble(study=x$study, geo_total=ifelse(x$study=="GSE79362",355,434), curated_hg38_total=x$raw_n,
         retained_labelled=nrow(x$meta), TBStatus=paste(names(table(x$meta$TBStatus)), as.integer(table(x$meta$TBStatus)), collapse="; "),
         sites=paste(names(table(x$meta$GeographicalRegion)), as.integer(table(x$meta$GeographicalRegion)), collapse="; "),
         unique_patients=length(unique(x$meta$PatientID)), repeated_patients=sum(table(x$meta$PatientID)>1),
         genes_after_filter=ifelse(x$study=="GSE79362", nrow(c1), nrow(c2)),
         baseline_deg=sum(deg$FDR<.05 & abs(deg$logFC)>1),
         baseline_up=sum(deg$FDR<.05 & deg$logFC>1), baseline_down=sum(deg$FDR<.05 & deg$logFC< -1),
         adjusted_fast_deg=sum(adj$FDR<.05 & abs(adj$logFC)>1))
}))
write.csv(summary, file.path(outdir, "analysis_summary.csv"), row.names=FALSE)

cross <- tibble(metric=c("shared_genes","spearman_logFC","shared_sig_DEGs","top100_overlap","hallmark_NES_spearman"),
                value=c(nrow(dc), cor(dc$logFC_79362,dc$logFC_94438,method="spearman"),
                        sum(dc$sig_79362 & dc$sig_94438),
                        length(intersect(d1 %>% arrange(PValue) %>% slice_head(n=100) %>% pull(gene),
                                         d2 %>% arrange(PValue) %>% slice_head(n=100) %>% pull(gene))),
                        cor(gs$NES_79362,gs$NES_94438,method="spearman")))
write.csv(cross, file.path(outdir, "cross_dataset_summary.csv"), row.names=FALSE)
message("Fast finish complete: ", outdir)
