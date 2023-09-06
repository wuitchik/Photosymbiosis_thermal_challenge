# Looking at the oxidative stress response
library(tidyverse)
library(GOfuncR)
library(gplots)
library(report)

# GO:0006979 -- Oxidative stress

children_oxidative = get_child_nodes("GO:0006979") 
children_oxidative_list = children_oxidative$child_go_id

# GO MWU
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_oxidative_list) %>%
  mutate(Treatment = "Heat",
         Phenotype = "Brown")

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Heat",
         Phenotype = "White")

Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "Brown")

Cold_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "White")

# Merged 

all = rbind(Heat_Sym_BP, Heat_apo_BP, Cold_Sym_BP, Cold_apo_BP) %>%
  mutate(Phenotype = as.factor(Phenotype)) %>%
  mutate(Treatment_phenotype = paste(Treatment, Phenotype, sep = " ")) %>%
  mutate(type=ifelse(Phenotype=="Brown","Highlighted","Normal")) 

test = aov(delta.rank ~ Treatment + Phenotype + Treatment:Phenotype, all)
summary(test)
report(test)

Tukey = TukeyHSD(test, conf.level = 0.95)
report(Tukey)

# plot

# Make colour pallet

# Points
shapes = c("White" = 1, "Brown" = 19)

# Colours
fill_cols = c("Cold White" = "white", "Cold Brown" = "#96bcfa","Heat Brown" = "#de835b","Heat White" = "white")
colour_cols = c("Cold White" = "#68a2ff", "Cold Brown" = "#68a2ff","Heat Brown" = "#ea6227","Heat White" = "#ea6227")


ggplot(data = all, aes(Treatment_phenotype, delta.rank, group = term)) +
  geom_point(aes(colour = Treatment_phenotype, shape = Phenotype), size = 3, stroke = 1.25) +
  geom_line() + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = "black") +
  stat_summary(fun = "mean", size = 0.5, colour = "black") +
  #geom_boxplot(aes(fill = Treatment_phenotype, colour = Treatment_phenotype), lwd = 1, alpha =.5) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = colour_cols) +
  scale_shape_manual(values = shapes) +
  ylab("'Oxidative stress' (GO:0006979) \n Delta-Rank") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") 

ggsave("Figure_4.pdf", plot = last_plot())
ggsave("Figure_4.jpg", plot = last_plot())


# Heatmap of oxidative stress

iso2go = read.delim("../GO_DEGs/Astrangia/astrangia_iso2go.tab", sep = "\t") 
list = paste(as.list(children_oxidative_list), collapse = "|")

gene_list = read.csv("../DESeq2/Astrangia/outlier_removed_host_counts.csv") %>%
  rename(Protein_ID = X) %>%
  left_join(iso2go) %>%
  filter(str_detect(GO_Swiss.Prot,list))

load("../DESeq2/Astrangia/sym_dds.RData")

library(DESeq2)

expDesign_grouped = expDesign %>%
  arrange(Treatment, Sym.Status)

order = as.list(expDesign_grouped$Sample)

vst = vst(dds.sym, blind = TRUE) 
vst = as.data.frame(assay(vst)) %>%
  rownames_to_column(var = "Protein_ID") %>%
  filter(Protein_ID %in% gene_list$Protein_ID) %>%
  column_to_rownames(var = "Protein_ID") %>%
  select(paste(order))

vst_means = apply(vst,1,mean)
vst_diff = (vst - vst_means) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pdf("heatmap.pdf")

heatmap.2(as.matrix(vst_diff),
            Rowv = TRUE,
            Colv = FALSE, 
          scale = "row",
          dendrogram = "both",
          trace = "none",
          col = colour,
          labRow = "")

dev.off()


