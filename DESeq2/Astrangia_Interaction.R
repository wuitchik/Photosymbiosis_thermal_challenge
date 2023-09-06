## Differential expression in Astrangia poculata
## Interaction between temperature and symbiont state
## Likelihood ratio test

library(DESeq2)
library(tidyverse)

counts = read.csv("Astrangia/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Astrangia/outlier_removed_expDesign.csv") %>%
  select(-X)

# Interaction between treatment and symbiotic state

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status + Treatment:Sym.Status)
dds = DESeq(dds)

dds = DESeq(dds,
            test="LRT",
            reduced = ~ Treatment + Sym.Status) 

tempXsym = results(dds) 

summary(tempXsym)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.003%
# outliers [1]       : 13, 0.039%
# low counts [2]     : 2, 0.0059%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

tempXsym = as.data.frame(tempXsym) %>%
  arrange(padj)

write.csv(tempXsym, "Astrangia/tempXsym_interaction.csv")

save.image("Astrangia/interaction_DDS.RData")
