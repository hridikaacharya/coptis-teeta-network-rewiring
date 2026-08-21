library(WGCNA)
library(igraph)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# expression matrix from previous script
expr <- logCPM[candidate_genes, ]

# one robust correlation network
cor_mat <- bicor(
  t(expr),
  use = "pairwise.complete.obs",
  maxPOutliers = 0.1
)

cor_mat[is.na(cor_mat)] <- 0

cutoff <- 0.7
adj <- abs(cor_mat)
adj[adj < cutoff] <- 0
diag(adj) <- 0

g <- graph_from_adjacency_matrix(
  adj,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

deg <- degree(g)

network_results <- data.frame(
  Gene = names(deg),
  Degree = deg
)

network_results <- network_results[
  order(network_results$Degree, decreasing = TRUE),
]

dir.create("results", showWarnings = FALSE)

write.csv(
  network_results,
  "results/network_centrality.csv",
  row.names = FALSE
)

head(network_results)


pathway <- read.csv("results/pathway_gene_level_candidates.csv")

network <- read.csv("results/network_centrality.csv")

merged <- merge(
  pathway,
  network,
  by.x = "Gene",
  by.y = "Gene"
)

merged <- merged[order(-merged$Degree), ]

head(merged, 20)

write.csv(
  merged,
  "results/pathway_network_hubs.csv",
  row.names = FALSE
)

top500 <- network_results[1:500, ]

write.csv(
  top500,
  "results/top500_hubs.csv",
  row.names = FALSE
)

list.files("results", full.names = TRUE)

