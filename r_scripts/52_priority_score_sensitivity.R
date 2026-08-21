library(dplyr)

###########################################################
# Load
###########################################################

master <- read.csv(
  "results/final_master_table.csv"
)

rewiring <- read.csv(
  "results/network_rewiring.csv"
)

final <-
  master %>%
  left_join(
    rewiring,
    by="Gene"
  )

###########################################################
# Three scoring models
###########################################################

model1 <-
  final %>%
  mutate(
    Score =
      Integrated_Score +
      abs(logFC_Drought) +
      abs(logFC_HeatVsDrought)/2 +
      abs(Degree_Difference)/1000
  ) %>%
  arrange(desc(Score))

model2 <-
  final %>%
  mutate(
    Score =
      Integrated_Score +
      1.5*abs(logFC_Drought) +
      abs(logFC_HeatVsDrought) +
      abs(Degree_Difference)/2000
  ) %>%
  arrange(desc(Score))

model3 <-
  final %>%
  mutate(
    Score =
      Integrated_Score +
      abs(logFC_Drought) +
      abs(logFC_HeatVsDrought)/2 +
      abs(Degree_Difference)/500
  ) %>%
  arrange(desc(Score))

###########################################################
# Save rankings
###########################################################

write.csv(model1,
          "results/ranking_model1.csv",
          row.names=FALSE)

write.csv(model2,
          "results/ranking_model2.csv",
          row.names=FALSE)

write.csv(model3,
          "results/ranking_model3.csv",
          row.names=FALSE)

###########################################################
# Top20 overlap
###########################################################

top1 <- model1$Gene[1:20]
top2 <- model2$Gene[1:20]
top3 <- model3$Gene[1:20]

cat("\nOverlap M1 vs M2:",
    length(intersect(top1,top2)),
    "/20\n")

cat("Overlap M1 vs M3:",
    length(intersect(top1,top3)),
    "/20\n")

cat("Overlap M2 vs M3:",
    length(intersect(top2,top3)),
    "/20\n")

###########################################################
# Spearman correlations
###########################################################

merge12 <-
  inner_join(
    model1[,c("Gene","Score")],
    model2[,c("Gene","Score")],
    by="Gene",
    suffix=c("_1","_2")
  )

merge13 <-
  inner_join(
    model1[,c("Gene","Score")],
    model3[,c("Gene","Score")],
    by="Gene",
    suffix=c("_1","_3")
  )

cat("\nSpearman M1 vs M2:\n")
print(cor.test(
  merge12$Score_1,
  merge12$Score_2,
  method="spearman"
))

cat("\nSpearman M1 vs M3:\n")
print(cor.test(
  merge13$Score_1,
  merge13$Score_3,
  method="spearman"
))

###########################################################
# Top 5 comparison
###########################################################

cat("\nTop5 Model1\n")
print(model1$Gene[1:5])

cat("\nTop5 Model2\n")
print(model2$Gene[1:5])

cat("\nTop5 Model3\n")
print(model3$Gene[1:5])