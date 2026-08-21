library(edgeR)
library(WGCNA)
library(igraph)
library(dplyr)

options(stringsAsFactors = FALSE)

#########################################################
# Load counts
#########################################################

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

#########################################################
# Normalize
#########################################################

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

#########################################################
# Candidate genes
#########################################################

ranking <- read.csv(
  "results/final_candidate_ranking.csv"
)

genes <- ranking$Gene

expr <- logCPM[genes,]

#########################################################
# Helper function
#########################################################

build_network <- function(mat){
  
  cor_mat <- bicor(
    t(mat),
    use="pairwise.complete.obs",
    maxPOutliers=0.1
  )
  
  cor_mat[is.na(cor_mat)] <- 0
  
  adj <- abs(cor_mat)
  
  adj[adj < 0.6] <- 0
  
  diag(adj) <- 0
  
  g <- graph_from_adjacency_matrix(
    adj,
    mode="undirected",
    weighted=TRUE,
    diag=FALSE
  )
  
  degree_tbl <- data.frame(
    
    Gene=names(degree(g)),
    Degree=degree(g)
    
  )
  
  degree_tbl
  
}

#########################################################
# Drought network
#########################################################

expr_drought <- expr[,1:6]

drought_degree <- build_network(expr_drought)

#########################################################
# Heat network
#########################################################

expr_heat <- expr[,c(1:3,7:9)]

heat_degree <- build_network(expr_heat)

#########################################################
# Merge
#########################################################

comparison <-
  
  drought_degree %>%
  
  rename(
    Degree_Drought=Degree
  ) %>%
  
  left_join(
    
    heat_degree %>%
      
      rename(
        Degree_Heat=Degree
      ),
    
    by="Gene"
    
  ) %>%
  
  mutate(
    
    Degree_Difference=
      Degree_Drought-Degree_Heat
    
  )

#########################################################
# Save
#########################################################

write.csv(
  comparison,
  "results/network_rewiring.csv",
  row.names=FALSE
)

#########################################################
# Top rewired genes
#########################################################

comparison %>%
  
  arrange(desc(abs(Degree_Difference))) %>%
  
  head(25) %>%
  
  print()