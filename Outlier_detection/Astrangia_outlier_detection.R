### Here, we explore for potential outlier samples in Astrangia poculata that result from poor sequencing or other such errors

library(DESeq)
library(arrayQualityMetrics)
library(tidyverse)

counts = read.csv("Astrangia/host_counts.csv", row.names = 1) %>%
  rename_all(funs(str_replace_all(.,".counts.txt", "")))

expDesign = read.csv("Astrangia/expDesign.csv")

## using OG DESeq (not DESEq2) and arrayQualityMetrics to isolate potential outliers

real = newCountDataSet(counts, expDesign) 
real = estimateSizeFactors(real)

cds = estimateDispersions(real,method="blind")
vsdBlind = varianceStabilizingTransformation(cds)
arrayQualityMetrics(vsdBlind,intgroup=c("Treatment"), force=TRUE, outdir = "Astrangia/arrayQualityMetrics") 

## Visit the .html output in Astrangia/arrayQualityMetrics/index.html
# from this, two outliers are identified: AF2, AS3 and removed in subsequent analyses

counts = counts %>%
  select(-AF2,-AS3,-AP4,-AC1)

expDesign = expDesign %>%
  filter(Sample %in% colnames(counts)) 

write.csv(counts, "Astrangia/outlier_removed_host_counts.csv")
write.csv(expDesign, "Astrangia/outlier_removed_expDesign.csv")



