library(dplyr)
library(igraph)


candidates <- read.csv(
  "results/top5_candidate_genes.csv"
)

network <- read.csv(
  "results/network_rewiring.csv"
)


top_genes <- candidates$Gene


# find candidate connections
candidate_network <- network %>%
  filter(
    Gene %in% top_genes
  )


write.csv(
  candidate_network,
  "results/top5_candidate_network.csv",
  row.names=FALSE
)


cat(
  "Candidate network rows:",
  nrow(candidate_network),
  "\n"
)

print(head(candidate_network))