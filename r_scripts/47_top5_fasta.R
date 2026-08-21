library(dplyr)

top5 <- read.csv(
  "results/top5_candidate_genes.csv"
)

proteins <- read.csv(
  "results/robust_hub_representative_proteins.csv"
)


top5_protein <- proteins %>%
  filter(
    Gene %in% top5$Gene
  )


write.csv(
  top5_protein,
  "results/top5_candidate_proteins.csv",
  row.names=FALSE
)


fasta <- c()

for(i in 1:nrow(top5_protein)){
  
  fasta <- c(
    fasta,
    paste0(">",top5_protein$Gene[i]),
    top5_protein$Protein_clean[i]
  )
}


writeLines(
  fasta,
  "results/top5_candidate_proteins.fasta"
)