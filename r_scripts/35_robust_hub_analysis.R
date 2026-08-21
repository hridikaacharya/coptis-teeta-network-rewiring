library(dplyr)
library(ggplot2)

master <- read.csv("results/final_master_table.csv")

#######################################################
# Robust hubs
#######################################################

robust <- master %>%
  filter(Appearances >= 3) %>%
  arrange(desc(Integrated_Score))

cat("Robust hubs:", nrow(robust), "\n")

write.csv(
  robust,
  "results/robust_hubs.csv",
  row.names = FALSE
)

#######################################################
# Summary
#######################################################

summary_table <- robust %>%
  summarise(
    Mean_Degree = mean(Degree),
    Mean_logFC = mean(abs(logFC_HeatVsDrought)),
    Mean_Score = mean(Integrated_Score)
  )

print(summary_table)

#######################################################
# Figure
#######################################################

p <- robust %>%
  slice(1:20) %>%
  ggplot(
    aes(
      x=reorder(Gene, Integrated_Score),
      y=Integrated_Score,
      fill=Appearances
    )
  )+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(
    title="Robust Stress-Regulatory Hub Genes",
    x="Gene",
    y="Integrated Score"
  )

ggsave(
  "04_Figures/Figure5_Robust_Hubs.pdf",
  p,
  width=9,
  height=7
)

ggsave(
  "04_Figures/Figure5_Robust_Hubs.png",
  p,
  width=9,
  height=7,
  dpi=300
)

print(p)