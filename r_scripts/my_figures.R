# ============================================================
# Coptis teeta Preliminary Findings
# Figure Generation Script
# ============================================================


# ------------------------------
# Working directory
# ------------------------------

setwd("D:/Lumiere/Project")

cat("Working directory:\n")
print(getwd())


# ------------------------------
# Packages
# ------------------------------

library(ggplot2)
library(dplyr)
library(readr)
library(igraph)


# ------------------------------
# Create figure folder
# ------------------------------

if (!dir.exists("04_Figures")) {
  dir.create("04_Figures")
}


# ------------------------------
# Check files
# ------------------------------

required_files <- c(
  "results/final_integrated_candidates.csv",
  "results/network_centrality.csv",
  "results/candidate_TF_pathway_network.csv"
)


cat("\nChecking files:\n")

for (f in required_files) {
  cat(f, " : ",
      file.exists(f),
      "\n")
}


if (!all(file.exists(required_files))) {
  stop("Missing required files. Check folder structure.")
}



# ============================================================
# Load data
# ============================================================


ranking <- read.csv(
  "results/final_integrated_candidates.csv"
)


centrality <- read.csv(
  "results/network_centrality.csv"
)


tf <- read.csv(
  "results/candidate_TF_pathway_network.csv"
)


cat("\nData dimensions:\n")

cat("Ranking:")
print(dim(ranking))

cat("Centrality:")
print(dim(centrality))

cat("TF network:")
print(dim(tf))



# ============================================================
# Figure 1
# Hub Stress Response
# ============================================================


top_hubs <- ranking %>%
  arrange(desc(Degree)) %>%
  head(30)


p1 <- ggplot(
  top_hubs,
  aes(
    x = reorder(Gene, Degree),
    y = Degree,
    fill = logFC_HeatVsDrought
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top Hub Genes Under Stress Transition",
    x = "Gene",
    y = "Network Degree",
    fill = "Heat vs Drought log2FC"
  )


ggsave(
  "04_Figures/Figure1_Hub_Stress_Response.png",
  p1,
  width = 10,
  height = 8,
  dpi = 300
)


ggsave(
  "04_Figures/Figure1_Hub_Stress_Response.pdf",
  p1,
  width = 10,
  height = 8
)


cat("\nFigure 1 saved\n")



# ============================================================
# Figure 2
# Differential Stress Landscape
# ============================================================


p2 <- ggplot(
  ranking,
  aes(
    x = logFC_Drought,
    y = logFC_Heat,
    color = Degree
  )
) +
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  scale_color_gradient(
    low = "lightblue",
    high = "darkblue"
  ) +
  theme_bw() +
  labs(
    title = "Differential Stress Response Landscape",
    x = "Drought log2 Fold Change",
    y = "Heat log2 Fold Change",
    color = "Network Degree"
  )


ggsave(
  "04_Figures/Figure2_Stress_Landscape.png",
  p2,
  width = 8,
  height = 7,
  dpi = 300
)


cat("Figure 2 saved\n")

# ============================================================
# Figure 3
# Candidate TF–Pathway Network
# ============================================================


tf_plot <- tf %>%
  arrange(desc(edge_score)) %>%
  head(25)


edges <- tf_plot %>%
  select(
    from = Gene.Names_x,
    to = Gene.Names_y,
    edge_score
  )


net <- graph_from_data_frame(
  edges,
  directed = TRUE
)


png(
  "04_Figures/Figure3_TF_Pathway_Network.png",
  width = 2500,
  height = 2000,
  res = 300
)


plot(
  net,
  layout = layout_with_fr(net),
  vertex.size = 20,
  vertex.label.cex = 0.8,
  vertex.color = "skyblue",
  edge.width = E(net)$edge_score / max(E(net)$edge_score) * 8,
  edge.arrow.size = 0.4,
  main = "Candidate TF–Alkaloid Biosynthesis Regulatory Network"
)


dev.off()


cat("Figure 3 saved\n")



# ============================================================
# Final check
# ============================================================


cat("\nGenerated figures:\n")

print(
  list.files(
    "04_Figures"
  )
)