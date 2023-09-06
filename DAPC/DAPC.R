# DAPC 

# libraries
library(DESeq2) 
library(adegenet)
library(tidyverse)

# load DESEq2 data
load("../DESeq2/Astrangia/sym_dds.RData")

# remove low counts
counts$low = apply(counts,1,
                   function(x)
                     {sum(x<=10)})

#35 is 90% of 39 samples - get rid of count less than 10 in more than 90% of samples
counts = counts[-which(counts$low>35),]
nrow(counts) # 8857 genes pass filter for host

# transform counts data using DESeq2
counts = counts %>%
  select(-low)

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status + Treatment:Sym.Status)

# rlog transform
rlog_counts = rlog(dds, blind = TRUE)
dat = as.data.frame(assay(rlog_counts))

# make pca
pcp=prcomp(t(dat), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores = pcp$x
screeplot(pcp, bstick = T)

# adegenet: finding clusters (even though we know what clusters we want) 
clus=find.clusters(t(dat),max.n.clus=30) #20PC, 3cluster
clus$grp = expDesign$Treatment_by_Sym.Status

# using a-score to find number of PCs
dapc_1 = dapc(t(dat), n.da = 100, n.pca= 50, clus$grp)
temp = optim.a.score(dapc_1) #optimal number of PCs = 10

# Shape pallet
shapes = c("White" = 1, "Brown" = 19)


# colours 

myCol = c("#68a2ff","#68a2ff","grey","grey","#ea6227","#ea6227")

# now lets build a discriminant function for these six groups:
dp=dapc(t(dat),clus$grp) # 10pcs, 6 discriminant functions

pdf("dapc.pdf")
scatter(dp, bg = "white", scree.da = TRUE, legend = TRUE, posi.leg = "topleft", pch = shapes,
        solid = 0.8, col = myCol, scree.pca = TRUE,
        clabel = 0)
dev.off()

save.image("DAPC.Rdata")


