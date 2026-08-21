library(dplyr)

############################################################
# Load data
############################################################

hubs <- read.csv("results/Table1_Robust_Hub_Genes.csv")

anno <- read.csv(
  "results/IFM89_function_annotation.csv",
  stringsAsFactors = FALSE
)

############################################################
# Check annotation columns
############################################################

cat("Annotation columns:\n")
print(names(anno))

############################################################
# Try automatic merge
############################################################

gene_col <- names(anno)[
  grepl("gene", names(anno), ignore.case = TRUE)
][1]

if (is.na(gene_col)) {
  
  stop("No gene column found in annotation table.")
  
}

annotation_col <- names(anno)[
  grepl("product|function|annotation|description",
        names(anno),
        ignore.case = TRUE)
][1]

if (is.na(annotation_col)) {
  
  annotation_col <- names(anno)[2]
  
}

############################################################
# Merge
############################################################

anno_small <-
  
  anno %>%
  
  select(
    
    Gene = all_of(gene_col),
    
    Annotation = all_of(annotation_col)
    
  )

merged <-
  
  hubs %>%
  
  left_join(
    
    anno_small,
    
    by = "Gene"
    
  )

############################################################
# Save
############################################################

write.csv(
  
  merged,
  
  "results/Table1_Robust_Hub_Genes_Annotated.csv",
  
  row.names = FALSE
  
)

############################################################
# Statistics
############################################################

cat("\nAnnotated genes:\n")

print(
  
  table(
    
    is.na(merged$Annotation)
    
  )
  
)

############################################################
# Show table
############################################################

print(
  
  merged %>%
    
    select(
      
      Gene,
      
      Annotation,
      
      Integrated_Score,
      
      Degree,
      
      logFC_Drought,
      
      logFC_Heat,
      
      logFC_HeatVsDrought
      
    )
  
)