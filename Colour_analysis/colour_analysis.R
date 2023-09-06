# Colour analysis
library(tidyverse)
library(performance)
library(cowplot)
library(ggpubr)
library(report)
library(ggpubr)


# load data
colour_data = read.csv("astrangia_relative_counts.csv")
load("../DESeq2/Astrangia/sym_dds.RData")

# add meta data to colour spreadsheet
data = colour_data %>%
  right_join(expDesign) %>%
  rename(Phenotype = Sym.Status) %>%
  select(-System, -Polyp_behaviour, -Day, -Genotype, -Treatment_by_Sym.Status)

# make sym count to host count ratio
data = data %>%
  mutate(Sym.Percent = ((Sym / Holobiont)*100))

# model selection
gaussian = glm(Sym.Percent ~ Treatment + Phenotype + Treatment:Phenotype, family = gaussian, data = data)
check_normality(gaussian) # Warning: Non-normality of residuals detected (p = 0.037).
check_heteroscedasticity(gaussian) # OK: Error variance appears to be homoscedastic (p = 0.152)

gaussian_sqrt_transform = glm(sqrt(Sym.Percent) ~ Treatment + Phenotype + Treatment:Phenotype, data = data)
check_normality(gaussian_sqrt_transform) # OK: residuals appear as normally distributed (p = 0.197).
check_heteroscedasticity(gaussian_sqrt_transform) # OK: Error variance appears to be homoscedastic (p = 0.934).
check_model(gaussian_sqrt_transform)
summary(gaussian_sqrt_transform)

# all:
#   glm(formula = sqrt(Sym.Percent) ~ Treatment + Phenotype + Treatment:Phenotype, 
#       data = data)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.0204  -0.2849  -0.1197   0.4092   1.0876  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      1.60767    0.17050   9.429  6.9e-11 ***
#   TreatmentControl                -0.22681    0.27492  -0.825 0.415283    
# TreatmentHeat                    0.02881    0.24958   0.115 0.908805    
# PhenotypeWhite                  -0.96638    0.24958  -3.872 0.000483 ***
#   TreatmentControl:PhenotypeWhite  0.45915    0.37131   1.237 0.224976    
# TreatmentHeat:PhenotypeWhite     0.37697    0.39198   0.962 0.343196    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.2325512)
# 
# Null deviance: 13.6315  on 38  degrees of freedom
# Residual deviance:  7.6742  on 33  degrees of freedom
# AIC: 61.275
# 
# Number of Fisher Scoring iterations: 2

report(gaussian_sqrt_transform)

# plot it
cols = c("White" = "grey", "Brown" = "orange4")

phenotype_plot = ggplot(data = data, aes(Phenotype, sqrt(Sym.Percent), fill = Phenotype)) +
  geom_jitter(aes(color = Phenotype), size = 2.5) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Phenotype",
       y = expression(sqrt(paste("% Symbiont Reads"))))+
  annotate("text", x = 1, y = 0.25, label = expression(P[paste("phenotype")] < 0.001)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_cowplot()

# link with carlos's colour analysis
carlos_data = read.csv("Astrangia_color.csv") %>%
  select(-Phenotype) %>%
  right_join(data, by = c("Sample")) 

# build linear model 
color_lm = lm(sqrt(Sym.Percent) ~ Random, data = carlos_data)
check_normality(color_lm) #OK: residuals appear as normally distributed (p = 0.546).
check_heteroscedasticity(color_lm) # OK: Error variance appears to be homoscedastic (p = 0.464).

check_model(color_lm)
summary(color_lm)

report(color_lm)

# plot it
correlation_plot = 
  ggplot(data = carlos_data, aes(Random, sqrt(Sym.Percent))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(aes(color = Phenotype), size = 2.5) +
  stat_cor(label.y = 2.3, label.x = 100,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 2.19, label.x = 100) +
  labs(x = "Red Channel Intensity (au)",
       y = expression(sqrt(paste("% Symbiont Reads"))))+
  scale_color_manual(values = cols) +
  theme_cowplot()

# combine plots

combined = phenotype_plot + rremove("legend") + correlation_plot + rremove("legend") 
ggsave("color_fig.pdf", combined, width = 8, height = 4, units = "in")
ggsave("color_fig.jpg", combined, width = 8, height = 4, units = "in")


