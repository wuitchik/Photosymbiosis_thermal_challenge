### Fv/Fm, for both Astrangia and Oculina Fv/Fm

library(tidyverse)
library(readxl)
library(patchwork)

cols = c("Brown" = "orange4", "White" = "grey90")

astrangia = read_excel("Astrangia_PAM.xlsx") %>%
  filter(day == 14) %>%
  mutate(FvFm = as.numeric(avgfvfm),
         sym_state = ifelse(sym_state == "B", "Brown", "White"))

a = ggplot(data = astrangia, aes(treatment, FvFm, fill = sym_state)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("A) Astrangia poculata") +
  ylim(0, 0.75) +
  theme_classic()

oculina = read_excel("Oculina_PAM.xlsx") %>%
  filter(day == 14) %>%
  mutate(FvFm = as.numeric(avgfvfm),
         sym_state = ifelse(sym_state == "Sym", "Brown", "White"))

o = ggplot(data = oculina, aes(treatment, FvFm, fill = sym_state)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("B) Oculina arbuscula") +
  ylim(0, 0.75) +
  theme_classic()

a + o + plot_layout(guides = "collect")

ggsave(plot = last_plot(), "FvFm.pdf")

# final - initial

astrangia2 = read_excel("Astrangia_PAM.xlsx") %>%
  filter(day == 14 | day == 5) %>%
  mutate(FvFm = as.numeric(avgfvfm),
         sym_state = ifelse(sym_state == "B", "Brown", "White")) %>%
  select(day, coral_id, sym_state, FvFm, treatment) %>%
  group_by(coral_id, sym_state, treatment) %>%
  summarize(delta_FvFm = last(FvFm) - first(FvFm))

a2 = ggplot(data = astrangia2, aes(treatment, delta_FvFm, fill = sym_state)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("A) Astrangia poculata") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylim(-0.4, 0.4) +
  theme_classic()

oculina2 = read_excel("Oculina_PAM.xlsx") %>%
  filter(day == 14 | day == 5) %>%
  mutate(FvFm = as.numeric(avgfvfm),
         sym_state = ifelse(sym_state == "Sym", "Brown", "White")) %>%
  select(day, coral_id, sym_state, FvFm, treatment) %>%
  group_by(coral_id, sym_state, treatment) %>%
  summarize(delta_FvFm = last(FvFm) - first(FvFm))

o2 = ggplot(data = oculina2, aes(treatment, delta_FvFm, fill = sym_state)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("B) Oculina arbuscula") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylim(-0.4, 0.4) +
  theme_classic()

a2 + o2 + plot_layout(guides = "collect")
ggsave(plot = last_plot(), "deltaFvFm.pdf")



