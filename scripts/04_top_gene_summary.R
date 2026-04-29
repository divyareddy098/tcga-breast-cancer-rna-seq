suppressPackageStartupMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

dir.create("results/biological", recursive = TRUE, showWarnings = FALSE)

convert_genes <- function(df) {
  df$ENSEMBL <- gsub("\\..*", "", df$gene)

  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(df$ENSEMBL),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME")
  )

  # Remove duplicate mappings
  gene_map <- gene_map %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)

  df <- dplyr::left_join(df, gene_map, by = "ENSEMBL")

  return(df)
}

comparisons <- c(
  "Basal_vs_LumA",
  "Her2_vs_LumA",
  "LumB_vs_LumA",
  "Normal_vs_LumA"
)

all_top <- data.frame()

for (comp in comparisons) {
  message("Processing: ", comp)

  input_file <- paste0("results/deseq2/top100_", comp, ".csv")

  df <- read.csv(input_file)
  df <- convert_genes(df)

  df_clean <- df %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::select(
      SYMBOL,
      GENENAME,
      ENSEMBL,
      log2FoldChange,
      padj,
      significant
    ) %>%
    dplyr::arrange(padj)

  write.csv(
    df_clean,
    paste0("results/biological/top_genes_", comp, ".csv"),
    row.names = FALSE
  )

  top20 <- df_clean %>% dplyr::slice_head(n = 20)
  top20$Comparison <- comp

  all_top <- rbind(all_top, top20)
}

write.csv(
  all_top,
  "results/biological/top20_all_subtypes.csv",
  row.names = FALSE
)

message("Top gene summary complete")
message("Saved outputs in results/biological/")
