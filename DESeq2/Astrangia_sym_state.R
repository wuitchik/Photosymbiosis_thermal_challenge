## Differential expression in Astrangia poculata
## Symbiont state
library(DESeq2)
library(tidyverse)

counts = read.csv("Astrangia/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Astrangia/outlier_removed_expDesign.csv") %>%
  select(-X)

# Symbiont state across all treatments

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status)
dds = DESeq(dds)

# Comparing how the white coral respond relative to brown corals

sym_results = results(dds, alpha = 0.05, contrast = c("Sym.Status", "White", "Brown"))
sym_summary = summary(sym_results)

# out of 33650 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 63, 0.19%
# LFC < 0 (down)     : 79, 0.23%
# outliers [1]       : 20, 0.059%
# low counts [2]     : 22426, 67%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


sym_results = as.data.frame(sym_results) %>%
  arrange(padj)

write.csv(sym_results, "Astrangia/sym_allTreatments_results.csv")

## Symbiont state in control

expDesign = expDesign %>%
  mutate(Treatment_by_Sym.Status = paste(Treatment,Sym.Status, sep = "_"))

dds.sym = DESeqDataSetFromMatrix(countData = counts,
                                 colData = expDesign,
                                 design = ~ Treatment_by_Sym.Status)
dds.sym = DESeq(dds.sym)

# save for other analyses
save(counts, expDesign, dds.sym, file = "./Astrangia/sym_dds.RData")

# Calculate results

control_sym_results = results(dds.sym,
                              alpha = 0.05,
                              contrast = c("Treatment_by_Sym.Status", "Control_White", "Control_Brown"))
control_sym_summary = summary(control_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 2, 0.0059%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

control_sym_results = as.data.frame(control_sym_results) %>%
  arrange(padj)

write.csv(control_sym_results, "Astrangia/sym_control_results.csv")

# Symbiont state in heat

heat_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Heat_White", "Heat_Brown"))
heat_sym_summary = summary(heat_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 8, 0.024%
# LFC < 0 (down)     : 91, 0.27%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 23711, 70%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

heat_sym_results = as.data.frame(heat_sym_results) %>%
  arrange(padj)

write.csv(heat_sym_results, "Astrangia/heat_sym_results.csv")

# Symbiont state in cold

cold_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Cold_White", "Cold_Brown"))
cold_sym_summary = summary(cold_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.003%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 2, 0.0059%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

cold_sym_results = as.data.frame(cold_sym_results) %>%
  arrange(padj)

write.csv(cold_sym_results, "Astrangia/cold_sym_results.csv")


# Brown control vs brown hot, controls go second

sym_control_vs_hot_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Heat_Brown", "Control_Brown"))
sym_control_vs_hot_summary = summary(sym_control_vs_hot_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 351, 1%
# LFC < 0 (down)     : 207, 0.62%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 24352, 72%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sym_control_vs_hot_results = as.data.frame(sym_control_vs_hot_results) %>%
  arrange(padj)

write.csv(sym_control_vs_hot_results, "Astrangia/sym_control_vs_hot_results.csv")


# White control vs white hot

apo_control_vs_hot_results = results(dds.sym,
                                     alpha = 0.05,
                                     contrast = c("Treatment_by_Sym.Status", "Heat_White", "Control_White"))
apo_control_vs_hot_summary = summary(apo_control_vs_hot_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 61, 0.18%
# LFC < 0 (down)     : 111, 0.33%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 26273, 78%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

apo_control_vs_hot_results = as.data.frame(apo_control_vs_hot_results) %>%
  arrange(padj)

write.csv(apo_control_vs_hot_results, "Astrangia/apo_control_vs_hot_results.csv")


# Brown control vs brown cold

sym_control_vs_cold_results = results(dds.sym,
                                     alpha = 0.05,
                                     contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Control_Brown"))
sym_control_vs_cold_summary = summary(sym_control_vs_cold_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1477, 4.4%
# LFC < 0 (down)     : 2023, 6%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 17947, 53%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sym_control_vs_cold_results = as.data.frame(sym_control_vs_cold_results) %>%
  arrange(padj)

write.csv(sym_control_vs_cold_results, "Astrangia/sym_control_vs_cold_results.csv")


# White control vs white hot

apo_control_vs_cold_results = results(dds.sym,
                                     alpha = 0.05,
                                     contrast = c("Treatment_by_Sym.Status", "Cold_White", "Control_White"))
apo_control_vs_cold_summary = summary(apo_control_vs_cold_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1814, 5.4%
# LFC < 0 (down)     : 2375, 7.1%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 16026, 48%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

apo_control_vs_cold_results = as.data.frame(apo_control_vs_cold_results) %>%
  arrange(padj)

write.csv(apo_control_vs_cold_results, "Astrangia/apo_control_vs_cold_results.csv")

