library(readxl)
library(janitor)

table_path <- "D:/Lumiere/Project/data/lmtables/"

P1 <- read_excel(file.path(table_path,"P1.xlsx"))
P2 <- read_excel(file.path(table_path,"P2.xlsx"))
P3 <- read_excel(file.path(table_path,"P3.xlsx"))
P4 <- read_excel(file.path(table_path,"P4.xlsx"))

fix_table <- function(df){
  
  names(df) <- as.character(df[1,])
  
  df <- df[-1,]
  
  df <- clean_names(df)
  
  rownames(df) <- NULL
  
  df
  
}

P1 <- fix_table(P1)
P2 <- fix_table(P2)
P3 <- fix_table(P3)
P4 <- fix_table(P4)

write.csv(P1,"results/P1_clean.csv",row.names=FALSE)
write.csv(P2,"results/P2_clean.csv",row.names=FALSE)
write.csv(P3,"results/P3_clean.csv",row.names=FALSE)
write.csv(P4,"results/P4_clean.csv",row.names=FALSE)

cat("Done\n")