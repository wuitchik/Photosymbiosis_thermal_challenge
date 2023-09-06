
library(tidyverse)
library(ggpubr)
library(performance)
library(report)

##### Plasticity of PCA space 

# load PCA data
load("../PCA/PCA.Rdata")
source("PCAplast_function.R") # Colleen's handy dandy plasticity function

expDesign = expDesign %>%
  rename(Phenotype = Sym.Status, 
         Phenotype_within_treatment = Treatment_by_Sym.Status)

pcaData = pcaData %>%
  rename(Phenotype_within_treatment = Treatment_by_Sym.Status) %>%
  select(-name)

# Colleen's plot

white_out = 
  PCAplast(pca = pcaData[,1:2],
         data = expDesign, 
         sample_ID = "Sample", 
         num_pca = "all", 
         control_col = "Phenotype_within_treatment",
         control_lvl = "Control_White"
         #group = "Genotype"
         ) %>%
  filter(Phenotype == "White")
  
brown_out = 
  PCAplast(pca = pcaData[,1:2],
           data = expDesign, 
           sample_ID = "Sample", 
           num_pca = "all", 
           control_col = "Phenotype_within_treatment",
           control_lvl = "Control_Brown"
           #group = "Genotype"
  ) %>%
  filter(Phenotype == "Brown")

plast_out = rbind(white_out, brown_out)


# Make colour pallet
cols = c("Control" = "grey", "Cold" = "#68a2ff", "Heat" = "#ea6227")

# Shape pallet
shapes = c("White" = 1, "Brown" = 19)

## Plot the plasticity (PC distances): overlay mean and 1 standard deviation
plast_plot <- ggplot(data = plast_out, aes(x = Phenotype_within_treatment, y = dist)) +
  geom_point(position = position_jitter(width = 0.1), aes(colour = Treatment, shape = Phenotype), size = 3, stroke = 1.25) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = "black") +
  stat_summary(fun = "mean", size = 0.5, colour = "black") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  theme_classic() +
  ylab("Plasticity")

ggsave("Plasticity.pdf", plot = last_plot())
ggsave("Plasticity.png", plot = last_plot())


pca_plot =
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

shapiro.test(plast_out$dist)
bartlett.test(dist ~ Phenotype_within_treatment, data = plast_out)

plasticity.aov = aov(dist ~ Phenotype_within_treatment, data = plast_out)
report(plasticity.aov)
check_model(plasticity.aov)

main_effect = aov(dist ~ Treatment + Phenotype + Phenotype:Treatment + Phenotype_within_treatment, data = plast_out)
report(main_effect)

TukeyHSD(plasticity.aov, conf.level = 0.95)
summary(plasticity.aov)

ggarrange(pca_plot, plast_plot, common.legend = TRUE, legend = "right")

ggsave("Plasticity_PCA.pdf", plot = last_plot())
ggsave("Plasticity_PCA.png", plot = last_plot())


#### Plasticity of Discriminant Function Space

# load PCA data
load("../DAPC/DAPC.Rdata")
source("PCAplast_function.R") # Colleen's handy dandy plasticity function

expDesign = expDesign %>%
  rename(Phenotype = Sym.Status, 
         Phenotype_within_treatment = Treatment_by_Sym.Status)

dapcData = dp$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  rename("PC1" = "LD1",
         "PC2" = "LD2") %>%
  left_join(expDesign)

rownames(dapcData) = dapcData$Sample

# Colleen's plot

white_out = 
  PCAplast(pca = dapcData[,2:3],
           data = expDesign, 
           #sample_ID = "Sample", 
           num_pca = "2", 
           control_col = "Phenotype_within_treatment",
           control_lvl = "Control_White"
           #group = "Genotype"
  ) %>%
  filter(Phenotype == "White")

brown_out = 
  PCAplast(pca = dapcData[,2:3],
           data = expDesign, 
           sample_ID = "Sample", 
           num_pca = "all", 
           control_col = "Phenotype_within_treatment",
           control_lvl = "Control_Brown"
           #group = "Genotype"
  ) %>%
  filter(Phenotype == "Brown")

plast_out = rbind(white_out, brown_out)


# Make colour pallet
cols = c("Control" = "grey", "Cold" = "#68a2ff", "Heat" = "#ea6227")

# Shape pallet
shapes = c("White" = 1, "Brown" = 19)

## Plot the plasticity (PC distances): overlay mean and 1 standard deviation

ggplot(data = dapcData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(x = PC1, y = PC2, colour = Treatment, shape = Phenotype), size = 3, stroke = 1.25) + 
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  theme_classic() 
  


plast_plot <- ggplot(data = plast_out, aes(x = Phenotype_within_treatment, y = dist)) +
  geom_point(position = position_jitter(width = 0.1), aes(colour = Treatment, shape = Phenotype), size = 3, stroke = 1.25) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = "black") +
  stat_summary(fun = "mean", size = 0.5, colour = "black") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  theme_classic() +
  ylab("Plasticity")

ggsave("Plasticity.pdf", plot = last_plot())
ggsave("Plasticity.png", plot = last_plot())


pca_plot =
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

shapiro.test(plast_out$dist)
bartlett.test(dist ~ Phenotype_within_treatment, data = plast_out)

plasticity.aov = aov(dist ~ Phenotype_within_treatment, data = plast_out)
report(plasticity.aov)
check_model(plasticity.aov)

main_effect = aov(dist ~ Treatment + Phenotype + Phenotype:Treatment + Phenotype_within_treatment, data = plast_out)
report(main_effect)

TukeyHSD(plasticity.aov, conf.level = 0.95)
summary(plasticity.aov)

ggarrange(pca_plot, plast_plot, common.legend = TRUE, legend = "right")

ggsave("Plasticity_PCA.pdf", plot = last_plot())
ggsave("Plasticity_PCA.png", plot = last_plot())
