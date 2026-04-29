suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

PADJ_CUTOFF <- 0.01
LOG2FC_CUTOFF <- 1.5

message("Loading processed data...")

counts <- read.csv(
  "data/processed/counts_filtered.csv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  "data/processed/metadata_filtered.csv",
  check.names = FALSE
)

common_samples <- intersect(colnames(counts), meta$sample_id)

if (length(common_samples) == 0) {
  stop("No matching samples found between counts and metadata.")
}

meta <- meta[match(common_samples, meta$sample_id), ]
counts <- counts[, common_samples, drop = FALSE]

stopifnot(all(colnames(counts) == meta$sample_id))

meta$Subtype <- factor(meta$Subtype)
meta$Subtype <- relevel(meta$Subtype, ref = "LumA")

message("Subtype counts:")
print(table(meta$Subtype))

message("Creating DESeq2 dataset...")

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ Subtype
)

dds <- dds[rowSums(counts(dds)) > 10, ]

message("Running DESeq2...")
dds <- DESeq(dds)

saveRDS(dds, "results/deseq2/dds_object.rds")

message("Running VST for PCA...")
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "results/deseq2/vst_object.rds")

pca <- plotPCA(vsd, intgroup = "Subtype", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

p <- ggplot(pca, aes(PC1, PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of TCGA-BRCA RNA-seq Samples by PAM50 Subtype") +
  theme_minimal()

ggsave("figures/pca_plot.png", p, width = 7, height = 5, dpi = 300)

comparisons <- list(
  c("Subtype", "Basal", "LumA"),
  c("Subtype", "Her2", "LumA"),
  c("Subtype", "LumB", "LumA"),
  c("Subtype", "Normal", "LumA")
)

summary_table <- data.frame()

for (comp in comparisons) {
  group1 <- comp[2]
  group2 <- comp[3]

  message("Running comparison: ", group1, " vs ", group2)

  res <- results(
    dds,
    contrast = comp,
    alpha = PADJ_CUTOFF
  )

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  res_df <- res_df %>%
    arrange(padj)

  res_df$significant <- ifelse(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_CUTOFF &
      abs(res_df$log2FoldChange) > LOG2FC_CUTOFF,
    "Significant",
    "Not significant"
  )

  outfile <- paste0("results/deseq2/", group1, "_vs_", group2, "_strict.csv")
  write.csv(res_df, outfile, row.names = FALSE)

  top_degs <- res_df %>%
    filter(significant == "Significant") %>%
    arrange(padj) %>%
    head(100)

  top_file <- paste0("results/deseq2/top100_", group1, "_vs_", group2, ".csv")
  write.csv(top_degs, top_file, row.names = FALSE)

  n_sig <- sum(res_df$significant == "Significant", na.rm = TRUE)
  n_up <- sum(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_CUTOFF &
      res_df$log2FoldChange > LOG2FC_CUTOFF,
    na.rm = TRUE
  )
  n_down <- sum(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_CUTOFF &
      res_df$log2FoldChange < -LOG2FC_CUTOFF,
    na.rm = TRUE
  )

  summary_table <- rbind(
    summary_table,
    data.frame(
      Comparison = paste0(group1, "_vs_", group2),
      Padj_Cutoff = PADJ_CUTOFF,
      Log2FC_Cutoff = LOG2FC_CUTOFF,
      Significant_DEGs = n_sig,
      Upregulated = n_up,
      Downregulated = n_down
    )
  )

  v <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(
      values = c("Not significant" = "grey70", "Significant" = "red")
    ) +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste(group1, "vs", group2),
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "DEG status"
    )

  plot_file <- paste0("figures/volcano_", group1, "_vs_", group2, "_strict.png")
  ggsave(plot_file, v, width = 7, height = 5, dpi = 300)
}

write.csv(summary_table, "results/deseq2/deseq2_summary_strict.csv", row.names = FALSE)

message("Advanced DESeq2 analysis complete")
print(summary_table)
message("Saved strict DEG results in results/deseq2/")
message("Saved PCA and volcano plots in figures/")
