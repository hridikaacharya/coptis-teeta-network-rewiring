library(dplyr)

master <- read.csv(
  "results/final_master_table.csv"
)

rewiring <- read.csv(
  "results/network_rewiring.csv"
)


final <- master %>%
  left_join(
    rewiring,
    by="Gene"
  )


candidates <- final %>%
  mutate(
    Priority_Score =
      Integrated_Score +
      abs(logFC_Drought) +
      abs(logFC_HeatVsDrought)/2 +
      abs(Degree_Difference)/1000
  ) %>%
  arrange(desc(Priority_Score)) %>%
  slice(1:5)


write.csv(
  candidates,
  "results/top5_candidate_genes.csv",
  row.names=FALSE
)


print(
  candidates[,c(
    "Gene",
    "Integrated_Score",
    "logFC_Drought",
    "logFC_HeatVsDrought",
    "Degree_Difference",
    "Priority_Score"
  )]
)