# Here I attribute GO annotations to specific eigengenes 

library(tidyverse)
library(WGCNA)

# Load workspace from Astranga_WGCNA.R
load("Astrangia_WGCNA.RData")

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(traits$Mapped_Sym);
names(weight) = "Brown Phenotype"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(input, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(input, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# intramodular analysis, pulling out eigengenes of interest

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

pdf("figures/Brown_membershivp_vs_significance.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for brown phenotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#####  GO MWU outputs
modules = sub("MM", "", colnames(geneModuleMembership))
genes = rownames(geneModuleMembership)
inputFiles = c()

for (m in modules){
  moduleGenes = genes[moduleColors==m]
  godf = data.frame(gene=genes,
                    inMod=if_else(genes %in% moduleGenes,
                                  1,
                                  0))
  outname = paste(c('./GO_MWU/', m, '_moduleInput.csv'), collapse='')
  print(paste('module = ', m))
  print(paste('total genes =', length(moduleGenes)))
  write.csv(godf, file=outname, row.names=F, quote=F)
  inputFiles = append(inputFiles, paste(m, '_moduleInput.csv', sep=''))
}

save(inputFiles, file='./GO_MWU/moduleInputFiles.Rdata')
