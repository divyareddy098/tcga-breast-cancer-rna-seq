suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

RAW_RDS <- "data/raw/tcga_brca_star_counts_se.rds"

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

message("Loading data...")
se <- readRDS(RAW_RDS)

counts <- assay(se, "unstranded")
meta_raw <- as.data.frame(colData(se))

meta_raw$sample_id <- colnames(counts)
meta_raw$patient_id <- substr(meta_raw$sample_id, 1, 12)
meta_raw$sample_type_code <- substr(meta_raw$sample_id, 14, 15)

# Keep primary tumor only
meta_raw <- meta_raw %>% filter(sample_type_code == "01")
counts <- counts[, meta_raw$sample_id, drop = FALSE]

# Correct PAM50 subtype column
meta_raw$Subtype <- meta_raw$paper_BRCA_Subtype_PAM50

# Keep samples with subtype
meta_raw <- meta_raw %>% filter(!is.na(Subtype), Subtype != "")
counts <- counts[, meta_raw$sample_id, drop = FALSE]

# Clean subtype names
meta_raw$Subtype <- dplyr::case_when(
  meta_raw$Subtype == "Luminal A" ~ "LumA",
  meta_raw$Subtype == "Luminal B" ~ "LumB",
  meta_raw$Subtype == "Basal-like" ~ "Basal",
  meta_raw$Subtype == "HER2-enriched" ~ "Her2",
  meta_raw$Subtype == "Normal-like" ~ "Normal",
  TRUE ~ meta_raw$Subtype
)

# Remove duplicate patients
meta_raw <- meta_raw %>% arrange(patient_id, sample_id)
meta_raw <- meta_raw[!duplicated(meta_raw$patient_id), ]
counts <- counts[, meta_raw$sample_id, drop = FALSE]

# Keep only simple metadata columns
meta <- data.frame(
  sample_id = meta_raw$sample_id,
  patient_id = meta_raw$patient_id,
  Subtype = meta_raw$Subtype,
  sample_type = meta_raw$sample_type,
  gender = meta_raw$gender,
  age_at_diagnosis = meta_raw$age_at_diagnosis,
  ajcc_pathologic_stage = meta_raw$ajcc_pathologic_stage,
  stringsAsFactors = FALSE
)

# Integer counts
counts <- round(counts)
storage.mode(counts) <- "integer"

# Filter low-expression genes
keep_genes <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep_genes, ]

# QC summaries
sample_qc <- data.frame(
  sample_id = colnames(counts_filtered),
  total_counts = colSums(counts_filtered),
  detected_genes = colSums(counts_filtered > 0)
)

gene_qc <- data.frame(
  gene_id = rownames(counts_filtered),
  total_counts = rowSums(counts_filtered),
  samples_detected = rowSums(counts_filtered > 0)
)

subtype_counts <- as.data.frame(table(meta$Subtype))
colnames(subtype_counts) <- c("Subtype", "N")

write.csv(counts_filtered, "data/processed/counts_filtered.csv")
write.csv(meta, "data/processed/metadata_filtered.csv", row.names = FALSE)
write.csv(sample_qc, "results/sample_qc_summary.csv", row.names = FALSE)
write.csv(gene_qc, "results/gene_qc_summary.csv", row.names = FALSE)
write.csv(subtype_counts, "results/subtype_counts.csv", row.names = FALSE)

message("DONE")
message("Genes retained: ", nrow(counts_filtered))
message("Samples retained: ", ncol(counts_filtered))
print(subtype_counts)
