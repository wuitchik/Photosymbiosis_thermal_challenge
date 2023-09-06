## Comparing GO Delta Ranks

# load libraries
library(tidyverse)
library(patchwork)
library(ggpubr)

# CHOOSE SET OF INPUT FILES TO LOOP RUN -----------------------------------
astrangia_data_path = "Astrangia/MWU_files/"   # path to the astrangia data
oculina_data_path = "Oculina/MWU_files/" # Oculina path

astrangia_results_files <- dir(astrangia_data_path, pattern = ".csv*") # get file names
oculina_results_files <- dir(oculina_data_path, pattern = ".csv*") # get file names
astrangia_names <- sub("\\.csv.*", "", astrangia_results_files)
oculina_names <- sub("\\.csv.*", "", oculina_results_files)

# Load all result files from folder
for(i in astrangia_names){
  filepath = file.path(astrangia_data_path,paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.table(filepath, header = T))
}

for(i in oculina_names){
  filepath = file.path(oculina_data_path,paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.table(filepath, header = T))
}

#### Biological Processes and sym state

MF_heat_sym_goods=intersect(astrangia_MF_heat_sym$term,oculina_MF_heat_sym$term)
astrangia = astrangia_MF_heat_sym[astrangia_MF_heat_sym$term %in% MF_heat_sym_goods,]
oculina = oculina_MF_heat_sym[oculina_MF_heat_sym$term %in% MF_heat_sym_goods,]

# Combine them
merged_data = merge(astrangia,oculina,by="term")


# This is to manually look for interesting go terms, and you can play with it in excel
MF_intersect = merged_data %>%
  filter(p.adj.x <0.1) %>%
  filter(p.adj.y < 0.1)

write.csv(MF_intersect, "MF_heat_sym_intersect.csv")

# Here is the actual plot,

MF_Heat_Sym =
  ggplot(merged_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  geom_point(alpha = 0.75) + 
  labs( x = "Astrangia Delta Rank",
        y = "Oculina Delta Rank") +
  labs(title = "MF Heat Sym State") +
  stat_regline_equation(label.y = 5.5, aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_smooth(method = "lm", color = "red") +
  theme_classic()

library(ggpubr)

(BP_Heat_Sym + MF_Heat_Sym + CC_Heat_Sym + plot_annotation(tag_levels = "A"))

