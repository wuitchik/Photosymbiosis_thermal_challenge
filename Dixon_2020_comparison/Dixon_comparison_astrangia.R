### Comparing the GO Delta Ranks to the stress module by Dixon et al. 2020
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggrepel)
library(GOfuncR)

# set ggrepel to not remove points due to overlap
options(ggrepel.max.overlaps = Inf)

# Read in Dixon's red module
Dixon_BP = read.table("TypeA", sep = " ", header = T)
Dixon_MF = read.table("Dixon_MF.csv", header = T)
Dixon_CC = read.table("Dixon_CC.csv", header = T)

# Read in Astrangia GO, Symbiont comparisons
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)

# Table of stress GO's 

stress = read.csv("stress_list.csv") %>%
  group_by(summary)

children_oxidative = get_child_nodes("GO:0006979") 
children_oxidative_list = children_oxidative$child_go_id



#### Heat Treatment ####


#### Brown Phenotype 
# Biological Processes
heat_sym_goods = intersect(Heat_Sym_BP$term, Dixon_BP$term)
heat_sym = Heat_Sym_BP[Heat_Sym_BP$term %in% heat_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_sym_goods,]

# Combine them
heat_sym_BP_data = merge(heat_sym, Dixon_BP_set, by="term")

#heat_sym_oxidative_data = heat_sym_BP_data[heat_sym_BP_data$term %in% children_oxidative_list,] 

# build hex plot with select Dixon GO's over top

heat_sym_BP_plot = ggplot(heat_sym_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  #geom_point(data = heat_sym_oxidative_data, size = 3) + 
  xlim(-6000,6000) +
  #geom_text_repel(data = heat_sym_oxidative_data, aes(label = name.y)) +
  theme_cowplot()

# Molecular Function
heat_sym_goods_MF = intersect(Heat_Sym_MF$term, Dixon_MF$term)
heat_sym_MF = Heat_Sym_MF[Heat_Sym_MF$term %in% heat_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% heat_sym_goods_MF,]

# Combine them
heat_sym_MF_data = merge(heat_sym_MF, Dixon_MF_set, by="term")

heat_sym_MF_plot = ggplot(heat_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
heat_sym_goods_CC = intersect(Heat_Sym_CC$term, Dixon_CC$term)
heat_sym_CC = Heat_Sym_CC[Heat_Sym_CC$term %in% heat_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% heat_sym_goods_CC,]

# Combine them
heat_sym_CC_data = merge(heat_sym_CC, Dixon_CC_set, by="term")

heat_sym_CC_plot = ggplot(heat_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()


#### White Phenotype 
# Biological Processes
heat_apo_goods = intersect(Heat_apo_BP$term, Dixon_BP$term)
heat_apo = Heat_apo_BP[Heat_apo_BP$term %in% heat_apo_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_apo_goods,]

# Combine them
heat_apo_BP_data = merge(heat_apo, Dixon_BP_set, by="term")

#heat_apo_oxidative_data = heat_apo_BP_data[heat_apo_BP_data$term %in% children_oxidative_list,] 


heat_apo_BP_plot = ggplot(heat_apo_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  xlim(-6000,6000) +
  #geom_point(data = heat_apo_oxidative_data, size = 3) + 
  #geom_text_repel(data = heat_apo_stress_data, aes(label = summary)) +
  theme_cowplot()


# Molecular Function
heat_apo_goods_MF = intersect(Heat_apo_MF$term, Dixon_MF$term)
heat_apo_MF = Heat_apo_MF[Heat_apo_MF$term %in% heat_apo_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% heat_apo_goods_MF,]

# Combine them
heat_apo_MF_data = merge(heat_apo_MF, Dixon_MF_set, by="term")

heat_apo_MF_plot = ggplot(heat_apo_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
heat_apo_goods_CC = intersect(Heat_apo_CC$term, Dixon_CC$term)
heat_apo_CC = Heat_apo_CC[Heat_apo_CC$term %in% heat_apo_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% heat_apo_goods_CC,]

# Combine them
heat_apo_CC_data = merge(heat_apo_CC, Dixon_CC_set, by="term")

heat_apo_CC_plot = ggplot(heat_apo_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

((heat_sym_BP_plot + heat_sym_MF_plot + heat_sym_CC_plot) /
 (heat_apo_BP_plot + heat_apo_MF_plot + heat_apo_CC_plot))

ggsave("Astrangia_heat_within_sym.pdf", 
       lplot(),
       width = 13,
       height = 7,
       units = "in")

#### Cold Treatment ####
# Read in Astrangia GO, Symbiont comparisons
Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_cold_results_modified_pvalues.csv", header = T)

Cold_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_cold_results_modified_pvalues.csv", header = T)

#### Brown Phenotype 

# Biological Processes
Cold_sym_goods = intersect(Cold_Sym_BP$term, Dixon_BP$term)
Cold_sym = Cold_Sym_BP[Cold_Sym_BP$term %in% Cold_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% Cold_sym_goods,]

# Combine them
Cold_sym_BP_data = merge(Cold_sym, Dixon_BP_set, by="term")
#Cold_sym_oxidative_data = Cold_sym_BP_data[Cold_sym_BP_data$term %in% children_oxidative_list,] 


Cold_sym_BP_plot = ggplot(Cold_sym_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  #geom_point(data = Cold_sym_oxidative_data, size = 3) + 
  #geom_text_repel(data = cold_sym_stress_data, aes(label = summary)) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  xlim(-6000,6000) +
  theme_cowplot()

# Molecular Function
Cold_sym_goods_MF = intersect(Cold_Sym_MF$term, Dixon_MF$term)
Cold_sym_MF = Cold_Sym_MF[Cold_Sym_MF$term %in% Cold_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% Cold_sym_goods_MF,]

# Combine them
Cold_sym_MF_data = merge(Cold_sym_MF, Dixon_MF_set, by="term")

Cold_sym_MF_plot = ggplot(Cold_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
Cold_sym_goods_CC = intersect(Cold_Sym_CC$term, Dixon_CC$term)
Cold_sym_CC = Cold_Sym_CC[Cold_Sym_CC$term %in% Cold_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% Cold_sym_goods_CC,]

# Combine them
Cold_sym_CC_data = merge(Cold_sym_CC, Dixon_CC_set, by="term")

Cold_sym_CC_plot = ggplot(Cold_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

#### White phenotype

# Biological Processes
Cold_apo_goods = intersect(Cold_apo_BP$term, Dixon_BP$term)
Cold_apo = Cold_apo_BP[Cold_apo_BP$term %in% Cold_apo_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% Cold_apo_goods,]

# Combine them
Cold_apo_BP_data = merge(Cold_apo, Dixon_BP_set, by="term")
#Cold_apo_oxidative_data = Cold_apo_BP_data[Cold_apo_BP_data$term %in% children_oxidative_list,] 


Cold_apo_BP_plot = ggplot(Cold_apo_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  #geom_point(data = Cold_apo_oxidative_data, size = 3) + 
  #geom_text_repel(data = cold_apo_stress_data, aes(label = summary)) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  xlim(-6000,6000) +
  theme_cowplot()

# Molecular Function
Cold_apo_goods_MF = intersect(Cold_apo_MF$term, Dixon_MF$term)
Cold_apo_MF = Cold_apo_MF[Cold_apo_MF$term %in% Cold_apo_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% Cold_apo_goods_MF,]

# Combine them
Cold_apo_MF_data = merge(Cold_apo_MF, Dixon_MF_set, by="term")

Cold_apo_MF_plot = ggplot(Cold_apo_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
Cold_apo_goods_CC = intersect(Cold_apo_CC$term, Dixon_CC$term)
Cold_apo_CC = Cold_apo_CC[Cold_apo_CC$term %in% Cold_apo_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% Cold_apo_goods_CC,]

# Combine them
Cold_apo_CC_data = merge(Cold_apo_CC, Dixon_CC_set, by="term")

Cold_apo_CC_plot = ggplot(Cold_apo_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon Delta Rank") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  xlim(-6000,6000) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_cowplot()

((Cold_sym_BP_plot + Cold_sym_MF_plot + Cold_sym_CC_plot) /
    (Cold_apo_BP_plot + Cold_apo_MF_plot + Cold_apo_CC_plot))

ggsave("Astrangia_cold_within_sym.pdf", 
       lplot(),
       width = 13,
       height = 7,
       units = "in")


# make biological process plot

Figure_3 = ((heat_sym_BP_plot | Cold_sym_BP_plot) / ( heat_apo_BP_plot | Cold_apo_BP_plot)) 

Figure_3 + plot_annotation(tag_levels = 'A')

ggsave("Figure_3.pdf", plot = last_plot(), width = 8, units = "in")
ggsave("Figure_3.jpg", plot = last_plot(), width = 8, units = "in")

