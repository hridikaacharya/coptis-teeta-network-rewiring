# ============================================
# Tests stability of hub-gene rankings across
# multiple correlation thresholds
# ============================================

library(WGCNA)
library(edgeR)
library(igraph)
library(dplyr)

options(stringsAsFactors = FALSE)

setwd("D:/Lumiere/Project")

dir.create("results/robustness", showWarnings = FALSE)

counts <- read.delim(
  "data/counts/gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
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

ranking <- read.csv(
  "results/final_candidate_ranking.csv"
)

candidate_genes <- ranking$Gene

expr <- logCPM[candidate_genes, ]

cor_mat <- bicor(
  t(expr),
  use="pairwise.complete.obs",
  maxPOutliers=0.1
)

cor_mat[is.na(cor_mat)] <- 0

diag(cor_mat) <- 0

thresholds <- c(
  0.60,
  0.70,
  0.80,
  0.90
)

robustness <- data.frame()

for(th in thresholds){
  
  adj <- abs(cor_mat)
  
  adj[adj < th] <- 0
  
  g <- graph_from_adjacency_matrix(
    adj,
    mode="undirected",
    weighted=TRUE,
    diag=FALSE
  )
  
  deg <- degree(g)
  
  top <- names(
    sort(
      deg,
      decreasing=TRUE
    )[1:50]
  )
  
  robustness <- rbind(
    robustness,
    data.frame(
      Threshold=th,
      Gene=top,
      Degree=deg[top]
    )
  )
  
  cat(
    "Finished threshold:",
    th,
    "\n"
  )
}

write.csv(
  robustness,
  "results/robustness/network_robustness.csv",
  row.names=FALSE
)

# ============================================
# Robust hub frequency
# ============================================

hub_frequency <- robustness |>
  count(Gene, name = "Appearances") |>
  arrange(desc(Appearances))

write.csv(
  hub_frequency,
  "results/robustness/hub_frequency.csv",
  row.names = FALSE
)

cat("\nMost stable hub genes:\n")
print(head(hub_frequency, 20))

library(ggplot2)

stable <- hub_frequency |>
  filter(Appearances >= 3)

p <- ggplot(
  stable,
  aes(
    x = reorder(Gene, Appearances),
    y = Appearances
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Robust Hub Genes Across Network Thresholds",
    x = "Gene",
    y = "Number of Thresholds"
  )

ggsave(
  filename = "04_Figures/Figure4_Robust_Hubs.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

print(p)
