### Venn Diagrams to compare DEGs between analyses and species

library(tidyverse)
library(VennDiagram)

# load results for temperature comparisons

heat_astrangia = read.csv("../DESeq2/Astrangia/heat_results.csv", row.names = 1)
heat_astrangia = row.names(heat_astrangia[heat_astrangia$padj<0.05 & !is.na(heat_astrangia$padj),])

cold_astrangia = read.csv("../DESeq2/Astrangia/cold_results.csv", row.names = 1) 
cold_astrangia = row.names(cold_astrangia[cold_astrangia$padj<0.05 & !is.na(cold_astrangia$padj),])

heat_oculina = read.csv("../DESeq2/Oculina/heat_results.csv", row.names = 1)
heat_oculina = row.names(heat_oculina[heat_oculina$padj<0.05 & !is.na(heat_oculina$padj),])

cold_oculina = read.csv("../DESeq2/Oculina/cold_results.csv", row.names = 1) 
cold_oculina = row.names(cold_oculina[cold_oculina$padj<0.05 & !is.na(cold_oculina$padj),])

# Astrangia temp comparisons

astrangia_temperature_list = list("Heat" = heat_astrangia, "Cold" = cold_astrangia)
astrangia_temp = venn.diagram(x = astrangia_temperature_list,
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
                              cat.dist = c(0.08, 0.08),
                              cat.pos = 1);
pdf(file = "Astrangia_HeatVsCold.pdf",
    height = 4,
    width = 4)
grid.draw(astrangia_temp)
dev.off()

# Oculina temp comparisons

oculina_temperature_list = list("Heat" = heat_oculina, "Cold" = cold_oculina)
oculina_temp = venn.diagram(x = oculina_temperature_list,
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
                              cat.dist = c(0.08, 0.08),
                              cat.pos = 1);
pdf(file = "Oculina_HeatVsCold.pdf",
    height = 4,
    width = 4)
grid.draw(oculina_temp)
dev.off()

# Heat species comparison

heat_species_comparison = list("Astrangia" = heat_astrangia, "Oculina" = heat_oculina)
heat_species = venn.diagram(x = heat_species_comparison,
                            filename=NULL,
                            col = "transparent",
                            fill = c("goldenrod1", "darkorchid1"),
                            alpha = 0.5,
                            cex = 2.5,
                            fontfamily = "sans",
                            fontface = "bold",
                            cat.default.pos = "text",
                            cat.col = "black",
                            cat.cex = 2.5,
                            cat.fontfamily = "sans",
                            cat.dist = c(0.08, 0.08),
                            cat.pos = 1);
pdf(file = "Heat_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(heat_species)
dev.off()

# Cold species comparison

cold_species_comparison = list("Astrangia" = cold_astrangia, "Oculina" = cold_oculina)
cold_species = venn.diagram(x = cold_species_comparison,
                            filename=NULL,
                            col = "transparent",
                            fill = c("goldenrod1", "darkorchid1"),
                            alpha = 0.5,
                            cex = 2.5,
                            fontfamily = "sans",
                            fontface = "bold",
                            cat.default.pos = "text",
                            cat.col = "black",
                            cat.cex = 2.5,
                            cat.fontfamily = "sans",
                            cat.dist = c(0.08, 0.08),
                            cat.pos = 1);
pdf(file = "Cold_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(cold_species)
dev.off()

# Symbiont state across all treatments (likelihood ratio test)

sym_astrangia = read.csv("../DESeq2/Astrangia/sym_allTreatments_results.csv", row.names = 1)
sym_astrangia = row.names(sym_astrangia[sym_astrangia$padj<0.05 & !is.na(sym_astrangia$padj),])

sym_oculina = read.csv("../DESeq2/Oculina/sym_allTreatments_results.csv", row.names = 1) 
sym_oculina = row.names(sym_oculina[sym_oculina$padj<0.05 & !is.na(sym_oculina$padj),])


sym_species_comparison = list("Astrangia" = sym_astrangia, "Oculina" = sym_oculina)
sym_species = venn.diagram(x = sym_species_comparison,
                            filename=NULL,
                            col = "transparent",
                            fill = c("goldenrod1", "darkorchid1"),
                            alpha = 0.5,
                            cex = 2.5,
                            fontfamily = "sans",
                            fontface = "bold",
                            cat.default.pos = "text",
                            cat.col = "black",
                            cat.cex = 2.5,
                            cat.fontfamily = "sans",
                            cat.dist = c(0.08, 0.08),
                            cat.pos = 1);
pdf(file = "sym_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(sym_species)
dev.off()

# Symbiont state in control

sym_control_astrangia = read.csv("../DESeq2/Astrangia/sym_control_results.csv", row.names = 1)
sym_control_astrangia = row.names(sym_control_astrangia[sym_control_astrangia$padj<0.05 & !is.na(sym_control_astrangia$padj),])

sym_control_oculina = read.csv("../DESeq2/Oculina/sym_control_results.csv", row.names = 1) 
sym_control_oculina = row.names(sym_control_oculina[sym_control_oculina$padj<0.05 & !is.na(sym_control_oculina$padj),])


sym_species_control_comparison = list("Astrangia" = sym_control_astrangia, "Oculina" = sym_control_oculina)
sym_control_species = venn.diagram(x = sym_species_control_comparison,
                           filename=NULL,
                           col = "transparent",
                           fill = c("goldenrod1", "darkorchid1"),
                           alpha = 0.5,
                           cex = 2.5,
                           fontfamily = "sans",
                           fontface = "bold",
                           cat.default.pos = "text",
                           cat.col = "black",
                           cat.cex = 2.5,
                           cat.fontfamily = "sans",
                           cat.dist = c(0.08, 0.08),
                           cat.pos = 1);
pdf(file = "sym_control_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(sym_control_species)
dev.off()

# Symbiont state in heat

sym_heat_astrangia = read.csv("../DESeq2/Astrangia/heat_sym_results.csv", row.names = 1)
sym_heat_astrangia = row.names(sym_heat_astrangia[sym_heat_astrangia$padj<0.05 & !is.na(sym_heat_astrangia$padj),])

sym_heat_oculina = read.csv("../DESeq2/Oculina/heat_sym_results.csv", row.names = 1) 
sym_heat_oculina = row.names(sym_heat_oculina[sym_heat_oculina$padj<0.05 & !is.na(sym_heat_oculina$padj),])


sym_species_heat_comparison = list("Astrangia" = sym_heat_astrangia, "Oculina" = sym_heat_oculina)
sym_heat_species = venn.diagram(x = sym_species_heat_comparison,
                                   filename=NULL,
                                   col = "transparent",
                                   fill = c("goldenrod1", "darkorchid1"),
                                   alpha = 0.5,
                                   cex = 2.5,
                                   fontfamily = "sans",
                                   fontface = "bold",
                                   cat.default.pos = "text",
                                   cat.col = "black",
                                   cat.cex = 2.5,
                                   cat.fontfamily = "sans",
                                   cat.dist = c(0.08, 0.08),
                                   cat.pos = 1);
pdf(file = "sym_heat_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(sym_heat_species)
dev.off()

# Symbiont state in cold

sym_cold_astrangia = read.csv("../DESeq2/Astrangia/cold_sym_results.csv", row.names = 1)
sym_cold_astrangia = row.names(sym_cold_astrangia[sym_cold_astrangia$padj<0.05 & !is.na(sym_cold_astrangia$padj),])

sym_cold_oculina = read.csv("../DESeq2/Oculina/cold_sym_results.csv", row.names = 1) 
sym_cold_oculina = row.names(sym_cold_oculina[sym_cold_oculina$padj<0.05 & !is.na(sym_cold_oculina$padj),])


sym_species_cold_comparison = list("Astrangia" = sym_cold_astrangia, "Oculina" = sym_cold_oculina)
sym_cold_species = venn.diagram(x = sym_species_cold_comparison,
                                filename=NULL,
                                col = "transparent",
                                fill = c("goldenrod1", "darkorchid1"),
                                alpha = 0.5,
                                cex = 2.5,
                                fontfamily = "sans",
                                fontface = "bold",
                                cat.default.pos = "text",
                                cat.col = "black",
                                cat.cex = 2.5,
                                cat.fontfamily = "sans",
                                cat.dist = c(0.08, 0.08),
                                cat.pos = 1);
pdf(file = "sym_cold_AvsO.pdf",
    height = 4,
    width = 4)
grid.draw(sym_cold_species)
dev.off()

