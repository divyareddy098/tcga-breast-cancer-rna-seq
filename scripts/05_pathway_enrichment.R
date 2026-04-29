suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/enrichment", recursive = TRUE, showWarnings = FALSE)

comparisons <- c(
  "Basal_vs_LumA",
  "Her2_vs_LumA",
  "LumB_vs_LumA",
  "Normal_vs_LumA"
)

for (comp in comparisons) {
  message("Running GO enrichment for: ", comp)

  df <- read.csv(paste0("results/biological/top_genes_", comp, ".csv"))

  genes <- df %>%
    filter(significant == "Significant") %>%
    pull(SYMBOL) %>%
    unique()

  entrez <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  ego <- enrichGO(
    gene = entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  out_csv <- paste0("results/enrichment/GO_BP_", comp, ".csv")
  write.csv(as.data.frame(ego), out_csv, row.names = FALSE)

  if (nrow(as.data.frame(ego)) > 0) {
    p <- dotplot(ego, showCategory = 15) +
      ggtitle(paste("GO Biological Process:", comp))

    ggsave(
      paste0("figures/enrichment/GO_BP_dotplot_", comp, ".png"),
      p,
      width = 9,
      height = 6,
      dpi = 300
    )
  }
}

message("✅ Pathway enrichment complete")
message("Saved tables in results/enrichment/")
message("Saved plots in figures/enrichment/")
