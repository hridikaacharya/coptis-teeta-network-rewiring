library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)


###################################################
# Load data
###################################################

specificity <- read.csv(
  "results/stress_specificity_table.csv"
)

rewiring <- read.csv(
  "results/hub_rewiring_integrated.csv"
)


###################################################
# Panel A
# Stress categories
###################################################

category <- specificity %>%
  count(Stress_Response)


p1 <- ggplot(
  category,
  aes(
    x=reorder(Stress_Response,n),
    y=n
  )
)+
  geom_col()+
  coord_flip()+
  labs(
    title="Stress response categories",
    x="",
    y="Number of genes"
  )+
  theme_bw()



###################################################
# Panel B
# Drought vs Heat expression
###################################################

p2 <- ggplot(
  specificity,
  aes(
    x=logFC_Drought,
    y=logFC_Heat,
    color=Stress_Response
  )
)+
  geom_point(
    alpha=0.4
  )+
  geom_hline(
    yintercept=0
  )+
  geom_vline(
    xintercept=0
  )+
  labs(
    title="Stress transcriptional landscape",
    x="Drought log2FC",
    y="Heat log2FC"
  )+
  theme_bw()



###################################################
# Panel C
# Network rewiring
###################################################

top_rewire <- rewiring %>%
  arrange(
    desc(abs(Degree_Difference))
  )%>%
  slice(1:20)



p3 <- ggplot(
  top_rewire,
  aes(
    x=reorder(Gene,
              Degree_Difference),
    y=Degree_Difference
  )
)+
  geom_col()+
  coord_flip()+
  labs(
    title="Top network rewiring genes",
    x="",
    y="Degree difference"
  )+
  theme_bw()



###################################################
# Combine
###################################################

figure <- p1 + p2 + p3 +
  plot_layout(
    ncol=3
  )


ggsave(
  "04_Figures/Figure6_Stress_Landscape.png",
  figure,
  width=16,
  height=6,
  dpi=300
)


ggsave(
  "04_Figures/Figure6_Stress_Landscape.pdf",
  figure,
  width=16,
  height=6
)


print(figure)