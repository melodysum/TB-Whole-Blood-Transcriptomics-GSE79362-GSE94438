suppressPackageStartupMessages(library(dplyr))
out <- "/private/tmp/tb_curated_results"
for (s in c("GSE79362", "GSE94438")) {
  cat("\n==", s, "==\n")
  print(read.csv(file.path(out, paste0(s, "_TBStatus_counts.csv"))))
  print(read.csv(file.path(out, paste0(s, "_site_by_TBStatus.csv"))))
  d <- read.csv(file.path(out, paste0(s, "_baseline_edgeR_DEG.csv")))
  cat("baseline DEG", sum(d$FDR < .05 & abs(d$logFC) > 1),
      "up", sum(d$FDR < .05 & d$logFC > 1),
      "down", sum(d$FDR < .05 & d$logFC < -1), "\n")
  a <- read.csv(file.path(out, paste0(s, "_adjusted_fast_DEG.csv")))
  cat("adjusted fast DEG", sum(a$FDR < .05 & abs(a$logFC) > 1), "\n")
  print(head(d[order(d$PValue), c("gene", "logFC", "PValue", "FDR")], 10))
}
cat("\nSignature\n")
print(read.csv(file.path(out, "signature_AUC_summary.csv")))
cat("\nGSEA key\n")
for (s in c("GSE79362", "GSE94438")) {
  g <- read.csv(file.path(out, paste0(s, "_GSEA_Hallmark.csv")))
  print(g %>%
          filter(grepl("INTERFERON|INFLAMMATORY|TNFA|COMPLEMENT|IL6", pathway)) %>%
          dplyr::select(study, pathway, NES, p.adjust) %>%
          arrange(p.adjust))
}
d1 <- read.csv(file.path(out, "GSE79362_baseline_edgeR_DEG.csv"))
d2 <- read.csv(file.path(out, "GSE94438_baseline_edgeR_DEG.csv"))
dc <- inner_join(d1, d2, by = "gene", suffix = c("_79362", "_94438"))
cat("\nshared genes", nrow(dc),
    "spearman", cor(dc$logFC_79362, dc$logFC_94438, method = "spearman"),
    "shared sig", sum(dc$FDR_79362 < .05 & abs(dc$logFC_79362) > 1 &
                        dc$FDR_94438 < .05 & abs(dc$logFC_94438) > 1),
    "top100 overlap", length(intersect(head(d1$gene[order(d1$PValue)], 100),
                                        head(d2$gene[order(d2$PValue)], 100))), "\n")
g1 <- read.csv(file.path(out, "GSE79362_GSEA_Hallmark.csv"))
g2 <- read.csv(file.path(out, "GSE94438_GSEA_Hallmark.csv"))
gs <- inner_join(g1, g2, by = "pathway", suffix = c("_79362", "_94438"))
cat("NES spearman", cor(gs$NES_79362, gs$NES_94438, method = "spearman"), "\n")
