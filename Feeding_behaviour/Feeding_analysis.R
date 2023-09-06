#### Feeding behaviour analyses using an ordinal logistic regression####

# libraries 
library(reshape2)
library(cowplot)
library(ordinal)
library(report)
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

#### Ordinal logistic regression treating genotype and system as random effects. ####

# Astrangia Full Model 
as_feeding_data = as_feeding %>%
  mutate(polyp_behaviour = factor(polyp_behaviour, levels = 1:5, ordered = T),
         Treatment = as.factor(Treatment),
         Treatment = relevel(Treatment, ref = "Control"),
         Sym.Status = as.factor(Sym.Status),
         Sym.Status = relevel(Sym.Status, ref = "White"),
         Genotype = as.factor(Genotype), 
         system = as.factor(system))

as_model = clmm(polyp_behaviour ~ Treatment * Sym.Status * Day + (1 | system) + (1 | Genotype) , data = as_feeding_data)
as_summary = summary(as_model) # note, object is saved for inspection at end of script

# post hoc test
library(emmeans)

post_hoc = emmeans(as_model, ~ Day+Sym.Status+Treatment)

pairs(post_hoc)

# Astrangia reduced model, no interaction
as_reduced_model = clmm(polyp_behaviour ~ Treatment + Sym.Status + (1 | system) + (1 | Genotype), data = as_feeding_data)
summary(as_reduced_model)



# Astrangia interaction by likelihood ratio test
as_comparison = anova(as_model, as_reduced_model)
summary(as_comparison)


## Oculina full Model 
oc_feeding_data = oc_feeding %>%
  mutate(polyp_behaviour = factor(polyp_behaviour, levels = 1:5, ordered = T),
         Treatment = as.factor(Treatment),
         Treatment = relevel(Treatment, ref = "Control"),
         Sym.Status = as.factor(Sym.Status),
         Sym.Status = relevel(Sym.Status, ref = "White"),
         Genotype = as.factor(Genotype), 
         system = as.factor(system))

oc_model = clmm(polyp_behaviour ~ Treatment * Sym.Status + (1 | system) + (1 | Genotype), data = oc_feeding_data)
oc_summary = summary(oc_model)

# Oculina reduced model
oc_reduced_model = clmm(polyp_behaviour ~ Treatment + Sym.Status + (1 | system) + (1 | Genotype), data = oc_feeding_data)
summary(oc_reduced_model)

# Oculina likelihood ratio test, comparing models
oc_comparison = anova(oc_model, oc_reduced_model)


# Write csv's of model results
write.csv(as_summary$coefficients, "Astrangia_files/astrangia_coefficients.csv")
write.csv(as_summary$info, "Astrangia_files/astrangia_info.csv")
write.csv(as_comparison, "Astrangia_files/astrangia_likelihood_ratio.csv")

write.csv(oc_summary$coefficients, "Oculina_files/oculina_coefficients.csv")
write.csv(oc_summary$info, "Oculina_files/oculina_info.csv")
write.csv(oc_comparison, "Oculina_files/oculina_likelihood_ratio.csv")






