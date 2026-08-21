library(dplyr)


proteins <- read.csv(
  "results/robust_hub_representative_proteins.csv",
  stringsAsFactors = FALSE
)


fasta <- c()


for(i in 1:nrow(proteins)){
  
  gene <- proteins$Gene[i]
  
  seq <- proteins$Protein_clean[i]
  
  
  fasta <- c(
    fasta,
    paste0(">",gene),
    seq
  )
  
}


writeLines(
  fasta,
  "results/robust_hub_proteins.fasta"
)


cat(
  "FASTA proteins:",
  nrow(proteins),
  "\n"
)