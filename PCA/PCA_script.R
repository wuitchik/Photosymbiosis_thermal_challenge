### PCAs

# load libraries 
library(tidyverse)
library(DESeq2)
library(vegan)

# load data
load(file = "../DESeq2/Astrangia/sym_dds.RData")

# transform data using vst 
vst = vst(dds.sym, blind = TRUE)

# calculate PCs
pcaData = DESeq2::plotPCA(vst, intgroup = c("Treatment_by_Sym.Status", "Polyp_behaviour", "Treatment", "Sym.Status"), returnData = TRUE)

# extract PC1 and PC2 variance
percentVar = round(100 * attr(pcaData, "percentVar"))

# Make PERMANOVA using adonis

pca <- prcomp(t(assay(vst)), center = TRUE, scale. = FALSE)
host_pca <- adonis2(pca$x ~ 
                     vst$Treatment:vst$Sym.Status +
                     vst$Sym.Status +
                     vst$Treatment,
                   method = 'eu')
host_pca

## Call:
##   adonis(formula = pca$x ~ vst$Treatment * vst$Sym.Status + vst$Sym.Status +      vst$Treatment, method = "eu") 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
## Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## vst$Treatment                 2    100325   50162  5.3256 0.22650  0.001 ***
##   vst$Sym.Status                1     14242   14242  1.5120 0.03215  0.044 *  
##   vst$Treatment:vst$Sym.Status  2     17537    8768  0.9309 0.03959  0.576    
## Residuals                    33    310831    9419         0.70175           
## Total                        38    442934                 1.00000           
## ---
##   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Make the plot

# Make colour pallet
cols = c("Control" = "grey", "Cold" = "#68a2ff", "Heat" = "#ea6227")

# Shape pallet
shapes = c("White" = 1, "Brown" = 19)

# pulling out p-values
pca_pval <- host_pca[["aov.tab"]][["Pr(>F)"]]

treatment_pca_pval <-substitute(italic(P[Treatment])==p, list(p = format(pca_pval[1], digits = 3)))
sym_state_pca_pval <-substitute(italic(P[Phenotype])==p, list(p = format(pca_pval[2], digits = 3)))

pcaData = pcaData %>%
  dplyr::rename(Phenotype = Sym.Status)


ggplot(pcaData, aes(PC1, PC2)) + 
  geom_point(aes(colour = Treatment, shape = Phenotype), size = 3, stroke = 1.25) +
  stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = Treatment)) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  annotate("text", x = 26, y = -43, label = deparse(sym_state_pca_pval), parse = TRUE, size = 3.5) +
  annotate("text", x = 26, y = -40, label = deparse(treatment_pca_pval), parse = TRUE, size = 3.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic()



ggsave("PCA.pdf", plot = last_plot())
ggsave("PCA.png", plot = last_plot())

save.image(file = "PCA.Rdata")
