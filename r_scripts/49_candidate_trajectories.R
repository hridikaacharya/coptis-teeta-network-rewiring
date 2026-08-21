library(dplyr)
library(ggplot2)
library(tidyr)


candidates <- read.csv(
  "results/top5_candidate_genes.csv"
)


plot_data <- candidates %>%
  select(
    Gene,
    logFC_Drought,
    logFC_Heat
  ) %>%
  pivot_longer(
    cols=c(
      logFC_Drought,
      logFC_Heat
    ),
    names_to="Condition",
    values_to="logFC"
  )


plot_data$Condition <- factor(
  plot_data$Condition,
  levels=c(
    "logFC_Drought",
    "logFC_Heat"
  )
)


p <- ggplot(
  plot_data,
  aes(
    x=Condition,
    y=logFC,
    group=Gene,
    color=Gene
  )
)+
  geom_line(size=1.2)+
  geom_point(size=3)+
  theme_minimal()+
  labs(
    title="Stress-response trajectories of top candidate regulators",
    x="Stress condition",
    y="log2 fold change"
  )


ggsave(
  "04_Figures/Figure7_Candidate_Trajectories.png",
  p,
  width=8,
  height=6,
  dpi=300
)


ggsave(
  "04_Figures/Figure7_Candidate_Trajectories.pdf",
  p,
  width=8,
  height=6
)

print(p)