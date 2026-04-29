library(dplyr)

dir.create("results/final", recursive = TRUE, showWarnings = FALSE)

files <- list.files("results/enrichment", pattern = "GO_BP_.*\\.csv$", full.names = TRUE)

summary_list <- lapply(files, function(f) {
  df <- read.csv(f)

  if (nrow(df) == 0) return(NULL)

  top_terms <- df %>%
    arrange(p.adjust) %>%
    head(5) %>%
    pull(Description)

  data.frame(
    Comparison = gsub("GO_BP_|\\.csv", "", basename(f)),
    Top_GO_Terms = paste(top_terms, collapse = "; ")
  )
})

summary_df <- bind_rows(summary_list)

write.csv(summary_df, "results/final/top_GO_terms_summary.csv", row.names = FALSE)

print(summary_df)
