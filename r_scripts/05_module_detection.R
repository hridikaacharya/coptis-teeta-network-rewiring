library(igraph)
library(WGCNA)

options(stringsAsFactors = FALSE)

library(edgeR)
library(WGCNA)
library(igraph)

# recreate expression matrix from saved counts

counts <- read.delim(
  "data/counts/gene_counts.txt",
  comment.char="#",
  check.names=FALSE
)

count_matrix <- counts[, -(1:6)]
rownames(count_matrix) <- counts$Geneid
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

colnames(count_matrix) <- c(
  "Control_1","Control_2","Control_3",
  "Drought_1","Drought_2","Drought_3",
  "Heat_1","Heat_2","Heat_3"
)

group <- factor(c(
  rep("Control",3),
  rep("Drought",3),
  rep("Heat",3)
))

dge <- DGEList(counts=count_matrix, group=group)

keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE)

# candidate genes come from final ranking

ranking <- read.csv("results/final_candidate_ranking.csv")

candidate_genes <- ranking$Gene

expr <- logCPM[candidate_genes, ]

# correlation network
cor_mat <- bicor(
  t(expr),
  use = "pairwise.complete.obs",
  maxPOutliers = 0.1
)

cor_mat[is.na(cor_mat)] <- 0

adj <- abs(cor_mat)
adj[adj < 0.95] <- 0
diag(adj) <- 0

cat("Non-zero entries in adjacency:",
    sum(adj > 0),
    "\n")

cat("Dimension:",
    dim(adj),
    "\n")

cat("Average degree estimate:",
    rowSums(adj > 0) |> mean(),
    "\n")

g <- graph_from_adjacency_matrix(
  adj,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

cat("Vertices:", vcount(g), "\n")
cat("Edges:", ecount(g), "\n")

cat("Connected components:",
    components(g)$no,
    "\n")

summary(E(g)$weight)

cat("Density:",
    edge_density(g),
    "\n")

# Leiden communities
modules <- cluster_leiden(
  g,
  weights = E(g)$weight
)

module_table <- data.frame(
  Gene = names(membership(modules)),
  Module = membership(modules)
)

dir.create("results", showWarnings = FALSE)

write.csv(
  module_table,
  "results/network_modules.csv",
  row.names = FALSE
)

sizes <- sort(sizes(modules), decreasing = TRUE)

print(head(sizes,20))

head(module_table)

table(module_table$Module)