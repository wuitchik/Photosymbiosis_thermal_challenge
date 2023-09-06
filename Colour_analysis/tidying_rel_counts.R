# tidy up tsv

library(tidyverse)
data = read_tsv("rel_counts.txt") %>%
  mutate(sample = gsub(".* ","", Sym)) %>%
  mutate(sample = gsub(".fastq_relative_Counts.tsv", "", sample)) %>%
  mutate(Sym = gsub(" .*", "", Sym))

write.csv(data, "astrangia_relative_counts.csv")
