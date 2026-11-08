library(dplyr)

robust <- read.csv(
  "results/robust_hubs.csv"
)

rewiring <- read.csv(
  "results/network_rewiring.csv"
)

top_rewire <- rewiring %>%
  arrange(desc(abs(Degree_Difference))) %>%
  slice(1:50)


# Universe
all_genes <- unique(rewiring$Gene)


a <- length(intersect(
  robust$Gene,
  top_rewire$Gene
))

b <- length(setdiff(
  robust$Gene,
  top_rewire$Gene
))

c <- length(setdiff(
  top_rewire$Gene,
  robust$Gene
))

d <- length(setdiff(
  all_genes,
  union(robust$Gene, top_rewire$Gene)
))


matrix <- matrix(
  c(a,b,c,d),
  nrow=2
)


print(matrix)

test <- fisher.test(matrix)

print(test)

###Most robust network hubs were distinct from stress-induced rewiring
###genes. However, their overlap was significantly enriched compared with random 
###expectation (Fisher's exact test, odds ratio = 7.17, p = 0.038), suggesting 
###that a subset of highly connected genes may act as structural regulators of 
###stress-driven network remodeling.
