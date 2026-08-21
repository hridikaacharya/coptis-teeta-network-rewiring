library(tidyverse)
library(edgeR)
library(WGCNA)
library(igraph)
library(pheatmap)

options(stringsAsFactors = FALSE)

counts <- read.delim(
  "data/counts/gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)
dim(counts)
head(counts)

gene_info <- counts[, 1:6]

count_matrix <- counts[, -(1:6)]

rownames(count_matrix) <- counts$Geneid

count_matrix <- as.matrix(count_matrix)

storage.mode(count_matrix) <- "integer"

dim(count_matrix)
colnames(count_matrix)

colnames(count_matrix) <- c(
  "Control_1",
  "Control_2",
  "Control_3",
  "Drought_1",
  "Drought_2",
  "Drought_3",
  "Heat_1",
  "Heat_2",
  "Heat_3"
)

colnames(count_matrix)

group <- factor(c(
  "Control","Control","Control",
  "Drought","Drought","Drought",
  "Heat","Heat","Heat"
))

library(edgeR)

dge <- DGEList(
  counts = count_matrix,
  group = group
)

keep <- filterByExpr(dge)

table(keep)

dge <- dge[keep, , keep.lib.sizes = FALSE]


dge <- calcNormFactors(dge)

dge$samples

table(keep) 

dim(dge)

logCPM <- cpm(dge, log = TRUE)

boxplot(
  logCPM,
  las = 2,
  col = c("skyblue","skyblue","skyblue",
          "tomato","tomato","tomato",
          "seagreen","seagreen","seagreen"),
  main = "Normalized logCPM distributions"
)

plotMDS(
  dge,
  col = c(rep("blue",3),
          rep("red",3),
          rep("darkgreen",3)),
  pch = 19
)

plotMDS(
  dge,
  labels = colnames(dge),
  col = c(rep("blue",3),
          rep("red",3),
          rep("darkgreen",3)),
  pch = 19
)

logCPM <- cpm(dge, log = TRUE)

gene_var <- apply(logCPM, 1, var)

summary(gene_var)

cutoff <- median(gene_var)

expr <- logCPM[gene_var > cutoff, ]

dim(expr)

# -----------------------------
# Differential expression setup
# -----------------------------

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

design

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)

contrast_drought <- makeContrasts(
  Drought - Control,
  levels = design
)

contrast_heat <- makeContrasts(
  Heat - Control,
  levels = design
)

contrast_heat_vs_drought <- makeContrasts(
  Heat - Drought,
  levels = design
)

de_drought <- glmQLFTest(fit, contrast = contrast_drought)

de_heat <- glmQLFTest(fit, contrast = contrast_heat)

de_hd <- glmQLFTest(fit, contrast = contrast_heat_vs_drought)

tab_drought <- topTags(
  de_drought,
  n = Inf
)$table

tab_heat <- topTags(
  de_heat,
  n = Inf
)$table

tab_hd <- topTags(
  de_hd,
  n = Inf
)$table

head(tab_drought)

head(tab_heat)

sig_drought <- rownames(tab_drought[tab_drought$FDR < 0.05, ])
sig_heat <- rownames(tab_heat[tab_heat$FDR < 0.05, ])
sig_hd <- rownames(tab_hd[tab_hd$FDR < 0.05, ])

candidate_genes <- unique(c(
  sig_drought,
  sig_heat,
  sig_hd
))

length(candidate_genes)

sig_drought <- rownames(tab_drought[tab_drought$FDR < 0.05, ])
sig_heat <- rownames(tab_heat[tab_heat$FDR < 0.05, ])
sig_hd <- rownames(tab_hd[tab_hd$FDR < 0.05, ])

candidate_genes <- unique(c(
  sig_drought,
  sig_heat,
  sig_hd
))

length(candidate_genes)

expr_network <- expr[rownames(expr) %in% candidate_genes, ]

dim(expr_network)

length(candidate_genes)

dim(expr_network)

head(rownames(expr_network))

dev.off()

pathway <- read.csv("results/pathway_genes_clean.csv")

dim(pathway)

head(pathway)

colnames(pathway)

pathway <- read.csv("results/pathway_genes_clean.csv")

dim(pathway)

colnames(pathway)

head(pathway)

graphics.off()

pathway <- read.csv("results/pathway_genes_clean.csv")

dim(pathway)
colnames(pathway)
head(pathway)

print(9)
