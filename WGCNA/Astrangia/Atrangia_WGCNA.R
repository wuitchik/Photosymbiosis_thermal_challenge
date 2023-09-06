### WGCNA on Astrangia poculata ###

library(DESeq2)
library(WGCNA)
library(tidyverse)


# Load dds object from previous analyses
load("../../DESeq2/Astrangia/temperature_DDS.RData")

# Variance stabilizing transformation, it's important to make sure blind = TRUE 
vst = vst(dds, blind = TRUE)

# Let WGCNA use multi threads (do only once)
allowWGCNAThreads()  

# transpose data frame
input = t(as.data.frame(assay(vst)))

# filter genes that are useful in WGCNA
gsg = goodSamplesGenes(input)
gsg$allOK #FALSE

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(input)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(input)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  input = input[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(input)
gsg$allOK #TRUE

# read in trait spreadsheet
traits = read.delim("astrangia_expDesign.csv", sep = ",", row.names = 1)

# cluster samples
sampleTree = hclust(dist(input), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
pdf("figures/Astrangia_dendro_and_traits.pdf")

plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = names(traits),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5, networkType = "signed")


#    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#      1 1.00e-01 10.7000          0.937 17000.0  17100.00  18200
#      2 5.78e-06  0.0441          0.962  8900.0   8880.00  10400
#      3 8.88e-02 -3.5500          0.857  4800.0   4760.00   6360
#      4 2.77e-01 -4.6200          0.811  2670.0   2630.00   4090
#      5 4.65e-01 -4.8100          0.827  1530.0   1490.00   2750
#      6 6.43e-01 -4.8300          0.884   907.0    873.00   1920
#      7 7.52e-01 -4.6300          0.928   554.0    521.00   1380
#      8 8.11e-01 -4.3700          0.957   348.0    321.00   1030
#      9 8.47e-01 -4.1200          0.974   227.0    202.00    789
#     10 8.53e-01 -3.8400          0.980   152.0    129.00    618
#     12 8.16e-01 -3.0600          0.933    76.6     56.20    403
#     14 8.91e-01 -2.0100          0.897    44.4     25.90    282
#     16 9.58e-01 -1.5700          0.951    29.5     12.60    241
#     18 9.49e-01 -1.3900          0.936    22.0      6.41    226
#     20 6.17e-01 -1.4700          0.538    18.0      3.37    215

# Plot the soft threshold scale independence
pdf("figures/scale_independence.pdf")
cex1 = 0.9;
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");
abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h

dev.off()


# Mean connectivity as a function of the soft-thresholding power

pdf("figures/Mean_connectivity.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


# automatic network construction and module detection, choosing a soft threshold power of 14
net = blockwiseModules(input,
                       power = 14,
                       networkType = "signed",
                       minModuleSize = 600,
                       numericLabels = FALSE,
                       mergeCutHeight = .25,
                       saveTOMFileBase = "Astrangia_TOM",
                       verbose = 3)

# number of genes associated with each module
table(net$colors)

#black      blue     brown     green      grey      pink       red turquoise    yellow 
# 650      1938      1734       849     23199       608       655      2770      1251 


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

MEList = moduleEigengenes(input, colors = net$unmergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")



# Define numbers of genes and samples
nGenes = ncol(input);
nSamples = nrow(input);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf("figures/Module_trait_relationship.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save.image(file = "Astrangia_WGCNA.RData")


