suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

message("Querying TCGA-BRCA raw count data from GDC...")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

message("Downloading data. This may take time...")
GDCdownload(query)

message("Preparing SummarizedExperiment object...")
se <- GDCprepare(query)

saveRDS(se, file = "data/raw/tcga_brca_star_counts_se.rds")

message("Saved file: data/raw/tcga_brca_star_counts_se.rds")
message("Genes: ", nrow(se))
message("Samples: ", ncol(se))
