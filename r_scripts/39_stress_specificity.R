library(dplyr)

master <- read.csv(
  "results/final_master_table.csv"
)

rewiring <- read.csv(
  "results/network_rewiring.csv"
)


specificity <- master %>%
  mutate(
    
    Stress_Response =
      
      case_when(
        
        abs(logFC_Drought) > 2 &
          abs(logFC_Heat) < 1
        ~ "Drought_specific",
        
        abs(logFC_Heat) > 2 &
          abs(logFC_Drought) < 1
        ~ "Heat_specific",
        
        abs(logFC_Drought) > 2 &
          abs(logFC_Heat) > 2
        ~ "Shared_stress",
        
        TRUE
        ~ "Weak_response"
        
      )
    
  )


write.csv(
  specificity,
  "results/stress_specificity_table.csv",
  row.names=FALSE
)


print(
  table(
    specificity$Stress_Response
  )
)