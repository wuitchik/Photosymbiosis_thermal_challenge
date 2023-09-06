#### Feeding behaviour analyses and figures ####

# libraries 
library(reshape2)
library(cowplot)
library(ggpubr)
library(tidyverse)


# data loading and tidying 
as_feeding = read.csv("Astrangia_files/Feeding_data.csv") %>%
  select(-Tank) %>%
  melt() %>%
  mutate(Day = as.integer(str_replace(variable, "Day.", ""))) %>%
  select(-X, -nubbin, -variable) %>%
  rename(system = treat) %>%
  rename(polyp_behaviour = value) %>%
  mutate(Sym.Status = ifelse(Sym.Status == "B", "Brown", "White")) %>%
  group_by(polyp_behaviour, Day, Treatment, Sym.Status, Genotype, system) %>%
  mutate(Genotype = substr(Genotype, 2,2)) %>%
  summarise(Frequency = n()) 

oc_feeding = read.csv("Oculina_files/Oculina Feeding Behavior - Sheet1.csv") %>%
  mutate(Genotype = substr(Coral.ID, 1,1)) %>%
  mutate(Sym.Status = ifelse(Symbiotic.status == "Sym", "Brown", "White")) %>%
  mutate(Treatment = case_when( 
    substr(Tank.ID, 1, 4) == "heat" ~ "Heat",
    substr(Tank.ID, 1, 4) == "cold" ~ "Cold",
    substr(Tank.ID, 1, 7) == "control" ~ "Control")) %>%
  select(-Symbiotic.status, -Coral.ID) %>%
  dplyr::rename(polyp_behaviour = Behavior) %>%
  dplyr::rename(system = Tank.ID) %>%
  group_by(polyp_behaviour, Day, Treatment, Sym.Status, Genotype, system) %>%
  dplyr::summarise(Frequency = n()) 

#### Plot it ####

#add colour pallet 
#extracting hex codes from different pallets
brewer.pal(10, "Blues")
brewer.pal(10, "Greys")
brewer.pal(5, "Reds")

as_color = as_feeding %>%
  mutate(colors = case_when(
    Treatment == "Heat" & polyp_behaviour == "1" ~ "#FEE5D9",
    Treatment == "Heat" & polyp_behaviour == "2" ~ "#FCAE91",
    Treatment == "Heat" & polyp_behaviour == "3" ~ "#FB6A4A",
    Treatment == "Heat" & polyp_behaviour == "4" ~ "#DE2D26",
    Treatment == "Heat" & polyp_behaviour == "5" ~ "#A50F15",
    Treatment == "Control" & polyp_behaviour == "5" ~ "#525252",
    Treatment == "Control" & polyp_behaviour == "4" ~ "#969696",
    Treatment == "Control" & polyp_behaviour == "3" ~ "#BDBDBD",
    Treatment == "Control" & polyp_behaviour == "2" ~ "#D9D9D9",
    Treatment == "Control" & polyp_behaviour == "1" ~ "#F0F0F0",
    Treatment == "Cold" & polyp_behaviour == "5" ~ "#08519C",
    Treatment == "Cold" & polyp_behaviour == "4" ~ "#3182BD",
    Treatment == "Cold" & polyp_behaviour == "3" ~ "#6BAED6",
    Treatment == "Cold" & polyp_behaviour == "2" ~ "#BDD7E7",
    Treatment == "Cold" & polyp_behaviour == "1" ~ "#C6DBEF",
  ))

as_color_codes <- as_color$colors
names(as_color_codes) <- as_color$colors

oc_color = oc_feeding %>%
  mutate(colors = case_when(
    Treatment == "Heat" & polyp_behaviour == "5" ~ "#A50F15",
    Treatment == "Heat" & polyp_behaviour == "4" ~ "#DE2D26",
    Treatment == "Heat" & polyp_behaviour == "3" ~ "#FB6A4A",
    Treatment == "Heat" & polyp_behaviour == "2" ~ "#FCAE91",
    Treatment == "Heat" & polyp_behaviour == "1" ~ "#FEE5D9",
    Treatment == "Control" & polyp_behaviour == "5" ~ "#525252",
    Treatment == "Control" & polyp_behaviour == "4" ~ "#969696",
    Treatment == "Control" & polyp_behaviour == "3" ~ "#BDBDBD",
    Treatment == "Control" & polyp_behaviour == "2" ~ "#D9D9D9",
    Treatment == "Control" & polyp_behaviour == "1" ~ "#F0F0F0",
    Treatment == "Cold" & polyp_behaviour == "5" ~ "#08519C",
    Treatment == "Cold" & polyp_behaviour == "4" ~ "#3182BD",
    Treatment == "Cold" & polyp_behaviour == "3" ~ "#6BAED6",
    Treatment == "Cold" & polyp_behaviour == "2" ~ "#BDD7E7",
    Treatment == "Cold" & polyp_behaviour == "1" ~ "#C6DBEF",
  ))
oc_color_codes <- oc_color$colors
names(oc_color_codes) <- oc_color$colors

as_behaviour_plot = 
  ggplot(as_color, aes(y = Frequency, x = Day, fill = as.factor(colors))) + 
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  facet_grid(Treatment ~ Sym.Status) + 
  labs(fill = "Behavioural Score") +
  scale_fill_manual(values = as_color_codes) +
  theme_cowplot()

oc_behaviour_plot = 
  ggplot(oc_color, aes(y = Frequency, x = Day, fill = as.factor(colors))) + 
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  facet_grid(Treatment ~ Sym.Status) + 
  labs(fill = "Behavioural Score") +
  scale_fill_manual(values = oc_color_codes) +
  theme_cowplot()

behaviour = ggarrange(as_behaviour_plot,
                      oc_behaviour_plot,
                      common.legend = TRUE,
                      legend = "none",
                      labels = c("A. Astrangia poculata", "B. Oculina arbuscula"),
                      vjust = 0.5,
                      hjust = -.25)

ggsave("Astrangia_Behaviour.png", plot = as_behaviour_plot, width = 8, height = 5, units = "in")
ggsave("Astrangia_Behaviour.pdf", plot = as_behaviour_plot, width = 8, height = 5, units = "in")

# note, legend is added in illustrator as it is a pain to make in R
