library(dplyr)

robust <- read.csv(
  "results/robust_hubs.csv"
)

rewiring <- read.csv(
  "results/network_rewiring.csv"
)


# top rewiring genes
top_rewire <- rewiring %>%
  arrange(desc(abs(Degree_Difference))) %>%
  slice(1:50)


overlap <- robust %>%
  filter(
    Gene %in% top_rewire$Gene
  )


cat(
  "Robust hubs:",
  nrow(robust),
  "\n"
)

cat(
  "Top rewiring genes:",
  nrow(top_rewire),
  "\n"
)

cat(
  "Overlap:",
  nrow(overlap),
  "\n"
)


write.csv(
  overlap,
  "results/hub_rewiring_overlap.csv",
  row.names=FALSE
)