## Delta ranks between cold and hot analyses

library(tidyverse)
library(plotly)
library(viridis)
library(scales)

# Load in both apo and sym comparisons
Hot_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Hot_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Hot_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Hot_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Hot_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Hot_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)

Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_cold_results_modified_pvalues.csv", header = T)

Cold_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Cold_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_cold_results_modified_pvalues.csv", header = T)


### Biological Processes
# merge into one data set

Hot_BP = merge(Hot_Sym_BP, Hot_apo_BP, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_MF = merge(Hot_Sym_MF, Hot_apo_MF, by = "term")%>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_CC = merge(Hot_Sym_CC, Hot_apo_CC, by = "term")%>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))


## experimental plot

scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))


Hot_BP_plot = 
  ggplot(data = Hot_BP, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "Biological Process") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
scale_color_viridis(option = "G", direction = -1) +
theme_cowplot()

ggplotly(Hot_BP_plot)

Hot_MF_plot = 
  ggplot(data = Hot_MF, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "Molecular Functions") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_cowplot()

ggplotly(Hot_MF_plot)

Hot_CC_plot = 
  ggplot(data = Hot_CC, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "Cellular Components") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_cowplot()

ggplotly(Hot_CC_plot)



  

