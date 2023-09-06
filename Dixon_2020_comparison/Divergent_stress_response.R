## How does the stress response change between apo and sym colonies

# First we define the stress response as that found by the red module in Dixon 2020
Dixon_BP = read.table("Dixon_BP.csv", header = T)
Dixon_MF = read.table("Dixon_MF.csv", header = T)
Dixon_CC = read.table("Dixon_CC.csv", header = T)

# Load in both apo and sym comparisons
Ast_Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Ast_Heat_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Ast_Heat_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Ast_Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Ast_Heat_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Ast_Heat_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)

### Biological Processes
# Filter those comparisons to only contain red module
Heat_sym_BP_stress = Ast_Heat_Sym_BP %>%
  filter(term %in% Dixon_BP$term)

Heat_apo_BP_stress = Ast_Heat_apo_BP %>%
  filter(term %in% Dixon_BP$term)

# Organize by most divergent response

BP_divergent = merge(Heat_sym_BP_stress, Heat_apo_BP_stress, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))


## experimental plot

BP =
  ggplot(data = BP_divergent, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "BP Dixon Stress Response") +
  scale_colour_continuous(type = "viridis") +
  theme_cowplot()

ggplotly(BP)

### Molecular functions

# Filter those comparisons to only contain red module
Heat_sym_MF_stress = Ast_Heat_Sym_MF %>%
  filter(term %in% Dixon_MF$term)

Heat_apo_MF_stress = Ast_Heat_apo_MF %>%
  filter(term %in% Dixon_MF$term)

# Organize by most divergent response

MF_divergent = merge(Heat_sym_MF_stress, Heat_apo_MF_stress, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))


## experimental plot

MF =
  ggplot(data = MF_divergent, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "MF Dixon Stress Response") +
  scale_colour_continuous(type = "viridis") +
  theme_cowplot()

ggplotly(MF)

### Cellular components

# Filter those comparisons to only contain red module
Heat_sym_CC_stress = Ast_Heat_Sym_CC %>%
  filter(term %in% Dixon_CC$term)

Heat_apo_CC_stress = Ast_Heat_apo_CC %>%
  filter(term %in% Dixon_CC$term)

# Organize by most divergent response

CC_divergent = merge(Heat_sym_CC_stress, Heat_apo_CC_stress, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))


## experimental plot

CC =
  ggplot(data = CC_divergent, aes(x = delta.rank.x, y = delta.rank.y, color = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Brown Delta-Rank",
       y = "White Delta-Rank",
       title = "CC Dixon Stress Response") +
  scale_colour_continuous(type = "viridis") +
  theme_cowplot()

ggplotly(CC)

