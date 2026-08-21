library(dplyr)
library(stringr)


proteins <- read.csv(
  "results/robust_hub_proteins.csv",
  stringsAsFactors = FALSE
)


# remove spaces/newlines from sequence

proteins <- proteins %>%
  mutate(
    Protein_clean = gsub(
      "\\s+",
      "",
      Protein
    ),
    Length = nchar(Protein_clean)
  )


# keep longest protein per gene

representative <- proteins %>%
  group_by(Gene) %>%
  slice_max(
    Length,
    n=1,
    with_ties = FALSE
  ) %>%
  ungroup()


write.csv(
  representative,
  "results/robust_hub_representative_proteins.csv",
  row.names=FALSE
)


cat(
  "Representative proteins:",
  nrow(representative),
  "\n"
)