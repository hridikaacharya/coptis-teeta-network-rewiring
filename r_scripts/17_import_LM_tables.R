library(readxl)
library(tidyverse)

# Change this if needed
table_path <- "D:/Lumiere/Project/data/lmtables/"

P1 <- read_excel(file.path(table_path, "P1.xlsx"))
P2 <- read_excel(file.path(table_path, "P2.xlsx"))
P3 <- read_excel(file.path(table_path, "P3.xlsx"))
P4 <- read_excel(file.path(table_path, "P4.xlsx"))

cat("P1:", dim(P1), "\n")
cat("P2:", dim(P2), "\n")
cat("P3:", dim(P3), "\n")
cat("P4:", dim(P4), "\n")

write.csv(P1, "results/P1.csv", row.names = FALSE)
write.csv(P2, "results/P2.csv", row.names = FALSE)
write.csv(P3, "results/P3.csv", row.names = FALSE)
write.csv(P4, "results/P4.csv", row.names = FALSE)