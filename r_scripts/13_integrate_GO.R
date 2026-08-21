library(tidyverse)


# Load ranking

ranking <- read.csv(
  "results/final_candidate_ranking.csv"
)


# Load GO enrichment

drought_bp <- read.csv(
  "results/Drought_down_BP.csv"
)


# Extract all drought GO genes

go_genes <- drought_bp$geneID %>%
  str_split("/") %>%
  unlist() %>%
  unique()


# Add GO flag

ranking_GO <- ranking %>%
  mutate(
    GO_drought_hit = Gene %in% go_genes
  )


# Save

write.csv(
  ranking_GO,
  "results/final_ranking_GO.csv",
  row.names=FALSE
)


cat(
  "GO genes:",
  length(go_genes),
  "\n"
)

cat(
  "Hub genes overlapping GO:",
  sum(ranking_GO$GO_drought_hit),
  "\n"
)