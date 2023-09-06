### Log fold change analysis 

# libraries
library(tidyverse)
library(patchwork)
library(ggpubr) 

# Load data
hot_sym = read.csv("../DESeq2/Astrangia/sym_control_vs_hot_results.csv") 
hot_apo = read.csv("../DESeq2/Astrangia/apo_control_vs_hot_results.csv")

cold_sym = read.csv("../DESeq2/Astrangia/sym_control_vs_cold_results.csv")
cold_apo = read.csv("../DESeq2/Astrangia/apo_control_vs_cold_results.csv")


# merge data sets
hot = hot_sym %>%
  merge(hot_apo, by="X") %>%
  rename(GeneID = X,
                Brown = "log2FoldChange.x",
                White = "log2FoldChange.y") %>%
  mutate_if(is.numeric , replace_na, replace = 1) %>%
  mutate(Colours = case_when(padj.x < 0.1 & padj.y < 0.1 ~ "both",
                             padj.x < 0.1 & padj.y > 0.1  ~ "brown",
                             padj.x > 0.1 & padj.y < 0.1 ~ "white",
                             padj.x > 0.1 & padj.y > 0.1 ~ "none"))

non_sig = hot %>%
  filter(Colours == "none")

brown_only = hot %>%
  filter(Colours == "brown")

white_only = hot %>%
  filter (Colours == "white")

both = hot %>%
  filter (Colours == "both")



ggplot(data = non_sig, aes(x = Brown, y = White)) +
  geom_point(alpha = 0.1) +
  geom_point(data = brown_only, colour = "brown") +
  geom_point(data = white_only, colour = "green") +
  geom_point(data = both, colour = "red") +
  labs(x = "Brown LFC",
       y = "White LFC") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_cowplot()


+
  #scale_fill_continuous() +
  labs(x = "Brown LFC",
       y = "White LFC") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
 # scale_color_viridis(option = "G", direction = -1) +
  theme_cowplot()

