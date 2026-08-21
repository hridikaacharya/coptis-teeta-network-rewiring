library(dplyr)

master <- read.csv("results/final_master_table.csv")

robust <- read.csv("results/robust_hubs.csv")

###########################################################
# Build publication table
###########################################################

table1 <-
  
  robust %>%
  
  select(
    
    Gene,
    
    Integrated_Score,
    
    Degree,
    
    Appearances,
    
    logFC_Drought,
    
    FDR_Drought,
    
    logFC_Heat,
    
    FDR_Heat,
    
    logFC_HeatVsDrought,
    
    FDR_HeatVsDrought,
    
    GO_drought_hit,
    
    Pathway_hit,
    
    TF_annotation,
    
    Product,
    
    Module
    
  ) %>%
  
  arrange(
    
    desc(Integrated_Score)
    
  )

###########################################################
# Save
###########################################################

write.csv(
  
  table1,
  
  "results/Table1_Robust_Hub_Genes.csv",
  
  row.names=FALSE
  
)

cat("Rows:",nrow(table1),"\n")

print(table1)