## Differential expression in Astrangia poculata
## Cold results
## Wald test


library(DESeq2)
library(tidyverse)

counts = read.csv("Astrangia/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Astrangia/outlier_removed_expDesign.csv") %>%
  select(-X)

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status)
dds = DESeq(dds)

## Cold 
cold_results = results(dds, alpha = 0.05, contrast = c("Treatment", "Cold", "Control"))
cold_summary = summary(cold_results)

# out of 33650 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2549, 7.6%
# LFC < 0 (down)     : 4141, 12%
# outliers [1]       : 20, 0.059%
# low counts [2]     : 14103, 42%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

cold_results = as.data.frame(cold_results) %>%
  arrange(padj)

write.csv(cold_results, "Astrangia/cold_results.csv")

## Heat 
heat_results = results(dds, alpha = 0.05, contrast = c("Treatment", "Heat", "Control"))
heat_summary = summary(heat_results)

# out of 33650 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 552, 1.6%
# LFC < 0 (down)     : 459, 1.4%
# outliers [1]       : 20, 0.059%
# low counts [2]     : 19865, 59%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

heat_results = as.data.frame(heat_results) %>%
  arrange(padj)

write.csv(heat_results, "Astrangia/heat_results.csv")

save.image("Astrangia/temperature_DDS.RData")