## Experimental design figure

# libraries
library(tidyverse)
library(cowplot)
library(patchwork)

# Data organization, Astrangia poculata = as, Oculina arbuscula = oc

as_temp = read.csv("Astrangia_Temperature.csv", header = T) %>%
  mutate(Treatment = case_when(
    treatment == "hot" ~ "Heat",
    treatment == "cold" ~ "Cold", 
    treatment == "control" ~ "Control")) %>%
  group_by(Treatment, day) %>%
dplyr::summarize(avg=mean(temperature))

oc_temp = read.csv("Oculina_temperature.csv", header = T) %>%
  mutate(Treatment = case_when(
    treatment == "hot" ~ "Heat",
    treatment == "cold" ~ "Cold", 
    treatment == "control" ~ "Control")) %>%
  group_by(Treatment, day) %>%
  select(-salinity, -X) %>%
  rename(temperature = temp_C) %>%
  summarize(avg=mean(temperature))

# Plot it

cols = c("Control" = "black",
         "Cold" = "dodgerblue1",
         "Heat" = "orangered1")

as_temp_plot = ggplot(as_temp, aes(x=day, y=avg, group=Treatment, colour = Treatment)) + 
  geom_line(size = 1.5) +
  scale_color_manual(values = cols) +
  ylim(-5,35) +
  scale_x_continuous(name = "Day", breaks = seq(0,16,2)) +
  labs(y="Temperature (°C)", x = "Day") +
  ggtitle("A) Astrangia poculata") +
  theme_cowplot()

oc_temp_plot = ggplot(oc_temp, aes(x=day, y=avg, group=Treatment, colour = Treatment)) + 
  geom_line(size = 1.5) +
  scale_color_manual(values = cols) +
  ylim(-5,35) +
  scale_x_continuous(name = "Day", breaks = seq(0,16,2)) +
  labs(y="Temperature (°C)", x = "Day") +
  ggtitle("B) Oculina arbuscula") +
  theme_cowplot()

ramp = as_temp_plot + oc_temp_plot + plot_layout(guides = "collect")
ggsave("Temperature_ramp.pdf", plot = ramp, width = 8, height = 4, units = "in")
            
            
            
            
                                                                     