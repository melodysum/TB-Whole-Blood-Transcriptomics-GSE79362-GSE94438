suppressPackageStartupMessages({
  library(MultiAssayExperiment); library(edgeR); library(dplyr); library(tibble)
  library(ggplot2); library(pROC); library(msigdbr); library(clusterProfiler); library(ggrepel)
})
outdir <- "/private/tmp/tb_curated_results"; figdir <- file.path(outdir, "figures")
load_lcpm <- function(study) {
  obj <- readRDS(file.path("/private/tmp", paste0(study, "_curatedTBData.rds")))[[study]]
  counts <- as.matrix(experiments(obj)[["assay_reprocess_hg38"]])
  meta <- as.data.frame(colData(obj)); meta$sample_id <- rownames(meta)
  common <- intersect(colnames(counts), meta$sample_id)
  counts <- counts[,common]; meta <- meta[match(common, meta$sample_id),]
  keep <- !is.na(meta$TBStatus); counts <- counts[,keep]; meta <- meta[keep,]
  counts <- rowsum(counts[rowSums(counts)>0,], rownames(counts)[rowSums(counts)>0], reorder=FALSE)
  min_n <- min(table(meta$TBStatus)); counts <- counts[rowSums(cpm(counts)>1) >= min_n,]
  y <- DGEList(counts); y <- calcNormFactors(y)
  list(meta=meta, logcpm=cpm(y, log=TRUE, prior.count=1), raw_n=ncol(as.matrix(experiments(obj)[["assay_reprocess_hg38"]])))
}
sig <- function(study, dat, genes, name) {
  meta <- dat$meta; present <- intersect(genes, rownames(dat$logcpm))
  meta$score <- colMeans(dat$logcpm[present,,drop=FALSE])
  ref <- if ("Control" %in% unique(meta$TBStatus)) "Control" else "LTBI"
  outcome <- factor(ifelse(meta$TBStatus == ref, ref, "PTB"), levels=c(ref,"PTB"))
  rocobj <- roc(outcome, meta$score, levels=c(ref,"PTB"), direction="<", quiet=TRUE)
  cc <- as.data.frame(t(coords(rocobj, "best", ret=c("threshold","sensitivity","specificity"), best.method="youden")))
  p <- ggplot(meta, aes(TBStatus, score, fill=TBStatus)) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=.15, alpha=.45, size=.9) + theme_bw() + theme(legend.position="none") +
    labs(title=paste0(study, ": ", name), y="Mean log2 CPM score", x=NULL)
  ggsave(file.path(figdir, paste0(study, "_", name, "_score.png")), p, width=5.5, height=4.2, dpi=220)
  if (length(unique(na.omit(meta$GeographicalRegion))) > 1) {
    ggsave(file.path(figdir, paste0(study, "_", name, "_score_by_site.png")),
           p + facet_wrap(~GeographicalRegion), width=8, height=4.5, dpi=220)
  }
  tibble(study=study, Signature=name, genes_present=length(present),
         missing_genes=paste(setdiff(genes, rownames(dat$logcpm)), collapse=";"),
         AUC=as.numeric(auc(rocobj)), Sensitivity=cc$sensitivity, Specificity=cc$specificity,
         WilcoxP=wilcox.test(meta$score ~ outcome)$p.value)
}
gsea_h <- function(study) {
  deg <- read.csv(file.path(outdir, paste0(study, "_baseline_edgeR_DEG.csv")))
  ranks <- deg %>% filter(PValue > 0) %>% mutate(score=sign(logFC)*-log10(PValue)) %>%
    group_by(gene) %>% slice_max(abs(score), n=1, with_ties=FALSE) %>% ungroup()
  gl <- ranks$score; names(gl) <- ranks$gene; gl <- sort(gl, decreasing=TRUE)
  h <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, gene_symbol)
  res <- GSEA(gl, TERM2GENE=h, pvalueCutoff=1, verbose=FALSE)
  df <- as.data.frame(res) %>% transmute(study=study, pathway=ID, NES, pvalue, p.adjust, leading_edge_genes=core_enrichment)
  write.csv(df, file.path(outdir, paste0(study, "_GSEA_Hallmark.csv")), row.names=FALSE)
  top <- df %>% arrange(p.adjust) %>% slice_head(n=18)
  p <- ggplot(top, aes(reorder(pathway,NES), NES, colour=p.adjust)) + geom_point(size=3) +
    coord_flip() + theme_bw() + scale_colour_viridis_c(direction=-1) + labs(title=paste0(study, ": Hallmark GSEA"), x=NULL, y="NES")
  ggsave(file.path(figdir, paste0(study, "_GSEA_Hallmark.png")), p, width=8, height=5.2, dpi=220)
  df
}
d1 <- load_lcpm("GSE79362"); d2 <- load_lcpm("GSE94438")
zak16 <- c("GBP5","BATF2","FCGR1B","SCARF1","TRAV27","ISG15","ANKRD22","ETV7","SERPING1","SAMD9L","IFIT2","IFIT3","IFI44L","CXCL10","HERC5","OAS1")
risk4 <- c("GBP5","SEPTIN4","CDO1","TRAV27")
gene11 <- c("GBP5","BATF2","FCGR1B","ANKRD22","ETV7","SERPING1","SAMD9L","IFI44L","CXCL10","HERC5","OAS1")
sigsum <- bind_rows(sig("GSE79362",d1,zak16,"Zak16"), sig("GSE94438",d2,zak16,"Zak16"),
                    sig("GSE79362",d1,risk4,"RISK4"), sig("GSE94438",d2,risk4,"RISK4"),
                    sig("GSE79362",d1,gene11,"Eleven_gene"), sig("GSE94438",d2,gene11,"Eleven_gene"))
write.csv(sigsum, file.path(outdir, "signature_AUC_summary.csv"), row.names=FALSE)
g1 <- gsea_h("GSE79362"); g2 <- gsea_h("GSE94438")
gs <- inner_join(g1,g2,by="pathway",suffix=c("_79362","_94438"))
write.csv(gs, file.path(outdir, "cross_dataset_Hallmark_GSEA_comparison.csv"), row.names=FALSE)
key <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_COMPLEMENT")
p <- ggplot(gs, aes(NES_79362,NES_94438,colour=(p.adjust_79362<.25 & p.adjust_94438<.25))) +
  geom_point(size=2.2) + geom_text_repel(data=gs %>% filter(pathway %in% key), aes(label=pathway), size=2.7) +
  geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=2)+theme_bw()+
  labs(title="Hallmark pathway NES consistency", x="GSE79362 NES", y="GSE94438 NES", colour="Both FDR<0.25")
ggsave(file.path(figdir, "cross_dataset_Hallmark_NES_scatter.png"), p, width=7, height=5.5, dpi=220)
message("signature/gsea only done")
