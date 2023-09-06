#Sym state within temp treatment venns 
library(tidyverse)
library(VennDiagram)

#### HOT comparisons ####
heat_brown = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv", row.names = 1)
heat_brown = row.names(heat_brown[heat_brown$padj<0.05 & !is.na(heat_brown$padj),])

heat_white = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv", row.names = 1)
heat_white = row.names(heat_white[heat_white$padj<0.05 & !is.na(heat_white$padj),])

# Make venn diagram and shared DEGs
heat_list = list("Brown" = heat_brown, "White" = heat_white)
heat_sym = venn.diagram(x = heat_list,
                              filename=NULL,
                              col = "transparent",
                              fill = c("orange4", "grey90"),
                              alpha = 0.5,
                              cex = 2.5,
                              fontfamily = "sans",
                              fontface = "bold",
                              cat.default.pos = "text",
                              cat.col = "black",
                              cat.cex = 2.5,
                              cat.fontfamily = "sans",
                              main = "Hot",
                              cat.dist = c(0.08, 0.08),
                              cat.pos = 1);
pdf(file = "heat_sym.pdf",
    height = 4,
    width = 4)
grid.draw(heat_sym)
dev.off()

# Hypergeometric test

a = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv")
h = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv") %>%
  filter(padj < 0.05) 

c = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv") %>%
  filter(padj < 0.05)

overlap = inner_join(h, c, by = "X")

phyper((nrow(overlap)-1), nrow(h), (nrow(a)-nrow(h)), nrow(c), lower.tail = F, log.p = FALSE)

#### COLD comparisons ####
cold_brown = read.csv("../DESeq2/Astrangia/sym_control_vs_cold_results.csv", row.names = 1)
cold_brown = row.names(cold_brown[cold_brown$padj<0.05 & !is.na(cold_brown$padj),])

cold_white = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv", row.names = 1)
cold_white = row.names(cold_white[cold_white$padj<0.05 & !is.na(cold_white$padj),])

# Make venn diagram and shared DEGs
cold_list = list("Brown" = cold_brown, "White" = cold_white)
cold_sym = venn.diagram(x = cold_list,
                        filename=NULL,
                        col = "transparent",
                        fill = c("orange4", "grey90"),
                        alpha = 0.5,
                        cex = 2.5,
                        fontfamily = "sans",
                        fontface = "bold",
                        cat.default.pos = "text",
                        cat.col = "black",
                        cat.cex = 2.5,
                        cat.fontfamily = "sans",
                        main = "Cold",
                        cat.dist = c(0.08, 0.08),
                        cat.pos = 1);
pdf(file = "cold_sym.pdf",
    height = 4,
    width = 4)
grid.draw(cold_sym)
dev.off()

write.csv(cold_list, "cold_venn.csv")


### Pulling out .csv of genes that are shared/different

brown_heat = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv") %>%
  filter(padj < 0.05)
white_heat = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv") %>%
  filter(padj < 0.05)

shared_heat = brown_heat %>%
  semi_join(white_heat, by = "X")

brown_only_heat = brown_heat %>%
  anti_join(white_heat, by = "X")

white_only_heat = white_heat %>%
  anti_join(brown_heat, by = "X")


brown_cold = read.csv("../DESeq2/Astrangia/sym_control_vs_cold_results.csv") %>%
  filter(padj < 0.05)
white_cold = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv") %>%
  filter(padj < 0.05)

shared_cold = brown_cold %>%
  semi_join(white_cold, by = "X")

brown_only_cold = brown_cold %>%
  anti_join(white_cold, by = "X")

white_only_cold = white_cold %>%
  anti_join(brown_cold, by = "X")


### Making .csv with annotations for exploration

annotations = read.delim("astrangia_annotations.txt", sep = "\t")

shared_heat_annotated = shared_heat %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

brown_only_heat_annotated = brown_only_heat %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

white_only_heat_annotated = white_only_heat %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

shared_cold_annotated = shared_cold %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

brown_only_cold_annotated = brown_only_cold %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

white_only_cold_annotated = white_only_cold %>%
  rename(Protein_ID = X) %>%
  left_join(annotations)

write.csv(shared_heat_annotated, "shared_heat.csv", row.names = F)
write.csv(brown_only_heat_annotated, "brown_only_heat.csv", row.names = F)
write.csv(white_only_heat_annotated, "white_only_heat.csv", row.names = F)
write.csv(shared_cold_annotated, "shared_cold.csv", row.names = F)
write.csv(brown_only_cold_annotated, "brown_only_cold.csv", row.names = F)
write.csv(white_only_cold_annotated, "white_only_cold.csv", row.names = F)

### Comparing Apo heat vs control agains apo cold vs control

heat_white = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv", row.names = 1)
heat_white = row.names(heat_white[heat_white$padj<0.05 & !is.na(heat_white$padj),])

cold_white = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv", row.names = 1)
cold_white = row.names(cold_white[cold_white$padj<0.05 & !is.na(cold_white$padj),])


# Make venn diagram and shared DEGs
heat_cold_apo_list = list("Heat" = heat_white, "Cold" = cold_white)
heat_cold_apo = venn.diagram(x = heat_cold_apo_list,
                        filename=NULL,
                        col = "transparent",
                        fill = c("#ea6227", "#68a2ff"),
                        alpha = 0.5,
                        cex = 2.5,
                        fontfamily = "sans",
                        fontface = "bold",
                        cat.default.pos = "text",
                        cat.col = "black",
                        cat.cex = 2.5,
                        cat.fontfamily = "sans",
                        main = "Hot",
                        cat.dist = c(0.08, 0.08),
                        cat.pos = 1);
pdf(file = "heat_cold_apo.pdf",
    height = 4,
    width = 4)
grid.draw(heat_cold_apo)
dev.off()

# Hypergeometric test

a = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv")
h = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv") %>%
  filter(padj < 0.05) 

c = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv") %>%
  filter(padj < 0.05)

overlap = inner_join(h, c, by = "X")

phyper((nrow(overlap)-1), nrow(h), (nrow(a)-nrow(h)), nrow(c), lower.tail = F, log.p = FALSE)
# [1] 1.245824e-40

### Comparing sym heat vs control agains sym cold vs control

heat_brown = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv", row.names = 1)
heat_brown = row.names(heat_brown[heat_brown$padj<0.05 & !is.na(heat_brown$padj),])

cold_brown = read.csv("../DESeq2/Astrangia/sym_control_vs_cold_results.csv", row.names = 1)
cold_brown = row.names(cold_brown[cold_brown$padj<0.05 & !is.na(cold_brown$padj),])


# Make venn diagram and shared DEGs
heat_cold_sym_list = list("Heat" = heat_brown, "Cold" = cold_brown)
heat_cold_sym = venn.diagram(x = heat_cold_sym_list,
                             filename=NULL,
                             col = "transparent",
                             fill = c("#ea6227", "#68a2ff"),
                             alpha = 0.5,
                             cex = 2.5,
                             fontfamily = "sans",
                             fontface = "bold",
                             cat.default.pos = "text",
                             cat.col = "black",
                             cat.cex = 2.5,
                             cat.fontfamily = "sans",
                             main = "Hot",
                             cat.dist = c(0.08, 0.08),
                             cat.pos = 1);
pdf(file = "heat_cold_sym.pdf",
    height = 4,
    width = 4)
grid.draw(heat_cold_sym)
dev.off()

# Hypergeometric test

a = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv")
h = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv") %>%
  filter(padj < 0.05) 

c = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv") %>%
  filter(padj < 0.05)

overlap = inner_join(h, c, by = "X")

phyper((nrow(overlap)-1), nrow(h), (nrow(a)-nrow(h)), nrow(c), lower.tail = F, log.p = FALSE)
# [1] 1.245824e-40
