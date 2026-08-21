library(dplyr)


hubs <- read.csv(
  "results/robust_hubs.csv"
)


proteins <- read.csv(
  "results/IFM89_protein_sequences.csv"
)


hub_proteins <- proteins %>%
  filter(
    Gene %in% hubs$Gene
  )


write.csv(
  hub_proteins,
  "results/robust_hub_proteins.csv",
  row.names=FALSE
)


cat(
  "Hub proteins:",
  nrow(hub_proteins)
)