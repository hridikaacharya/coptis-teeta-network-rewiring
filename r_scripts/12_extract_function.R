library(tidyverse)


gb <- readLines(
  "D:/Lumiere/ncbi_dataset/data/GCA_015680905.1/genomic.gbff"
)


# Combine lines into one text block
gb_text <- paste(gb, collapse="\n")


# Split each gene feature
features <- str_split(
  gb_text,
  "\\n     gene"
)[[1]]


annotation <- map_dfr(
  features,
  function(x){
    
    gene <- str_extract(
      x,
      "IFM89_[0-9]+"
    )
    
    product <- str_extract(
      x,
      "/product=\"[^\"]+\""
    )
    
    tibble(
      Gene = gene,
      Product = product
    )
  }
)


annotation <- annotation %>%
  filter(!is.na(Gene))


annotation$Product <- gsub(
  "/product=\"|\"",
  "",
  annotation$Product
)


dim(annotation)

head(annotation)


write.csv(
  annotation,
  "results/IFM89_function_annotation.csv",
  row.names=FALSE
)

annotation %>%
  filter(
    grepl(
      "berberine|norcoclaurine|bridge|stylopine|methyltransferase|OMT|coclaurine",
      Product,
      ignore.case=TRUE
    )
  )