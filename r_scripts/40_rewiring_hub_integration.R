library(dplyr)


master <- read.csv(
  "results/final_master_table.csv"
)


rewiring <- read.csv(
  "results/network_rewiring.csv"
)


rewired_master <- master %>%
  left_join(
    rewiring,
    by="Gene"
  )


rewired_master <- rewired_master %>%
  mutate(
    
    Rewiring_Category =
      
      case_when(
        
        Degree_Difference > 1000 &
          Integrated_Score > 2
        ~ "High_confidence_rewired_hub",
        
        Degree_Difference > 1000
        ~ "Strong_rewiring",
        
        Integrated_Score > 2
        ~ "Robust_hub",
        
        TRUE
        ~ "Other"
        
      )
    
  )


write.csv(
  rewired_master,
  "results/hub_rewiring_integrated.csv",
  row.names=FALSE
)


print(
  table(
    rewired_master$Rewiring_Category
  )
)