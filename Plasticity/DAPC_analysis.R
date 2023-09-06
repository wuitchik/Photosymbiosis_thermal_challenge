# Discriminant analysis of principal components

library(tidyverse)
library(adegenet)
library(DESeq2)

# load dds data
load("../DESeq2/Astrangia/sym_dds.RData")
expDesign = expDesign %>%
  dplyr::rename(Phenotype = Sym.Status) %>%
  dplyr::rename(Phenotype_within_treatment = Treatment_by_Sym.Status)

#Remove low count genes
counts$low = apply(counts,1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 

#35 is approximately 90% of 39 samples - get rid of count less than 10 in more than 90% of samples
counts = counts[-which(counts$low>37),]
nrow(counts) #10627 genes pass filter for host

counts = counts %>%
  select(-low)

# re-run DESEqobject
dds = DESeqDataSetFromMatrix(countData = counts, 
                                colData = expDesign,
                                design = ~ Treatment + Phenotype + Treatment:Phenotype)

# Transform the data
vst = vst(dds, blind = TRUE)
dat = as.data.frame(assay(vst))
colnames(dat) = names(counts)


# find the number of clusters 

clus = find.clusters(t(dat), n.pca = 50, n.clust = 6) #found out that 15 pcs, saw a relatively small BIC around 6 clusters


expDesign = expDesign %>%
  mutate(Levels = paste(Phenotype, Treatment))
clus$grp = expDesign$Levels
dapcall = dapc(t(dat),clus$grp, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall) 

# optimal number of PCs: 8 based on a-score optimisation
dapc1 <- dapc(t(dat), clus$grp, n.pca = 8, n.da = 6 )

myCol = c("Brown Cold" = "dodgerblue4",
          "Brown Control" = "grey40",
          "Brown Heat" = "firebrick",
          "White Cold" = "dodgerblue1",
          "White Control" = "grey80",
          "White Heat" = "firebrick1")

scatter(dapc1,2,1, col=myCol, legend = TRUE)
