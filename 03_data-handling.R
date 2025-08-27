################################################################################
#' Data Handling and Plotting for Simulation Model
#'
#' Author: Sam Appleyard-Sanderson
#'
#' Date: 2025-07-08
#' 
###############################################################################
# Referencing 
knitr::write_bib(c(.packages(), "bookdown"), "packages.bib")
packageVersion("ggh4x")
citation(package = "stringr", lib.loc = NULL, auto = NULL)
toBibtex(citation("stringr"))

# --- Load packages ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(viridis)
library(readr)
library(stringr)
library(patchwork)
library(ggh4x)
library(ragg)
# --- Theme for all graphs in the thesis ---------------------------------------
# Publication-ready theme for A4 single plots
theme_pub_a4 <- theme_minimal(base_size = 16) +
  theme(
    # Titles and labels
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 10, color = "black"), # tick labels
    
    # Gridlines
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # Background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Legend
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = NA),
    
    # Margins (good for A4 printing)
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

theme_set(theme_pub_a4)

guild_colors <- c("1" = "#440154FF", "2" = "#21908CFF", "3" = "#FDE725FF")

# -- Labels --------------------------------------------------------------------
life_stage_labels <- c(
  "Pat" = "Adult Pat",
  "Pupating" = "Pupating",
  "Field" = "Adult Field",
  "Feeding" = "Larval Feeding"
)

guild_labels <- c(
  "1" = "Fly",
  "2" = "Endocoprid beetle",
  "3" = "Paracoprid beetle"
)

stocking_density_labels <- c(
  "1" = "5",
  "2" = "10",
  "3" = "15"
)

# --- Figure 1 (Average community size) -------------------------------------
## Setup
plot_data_1 <- treatment_results_44_adults_die %>%
  filter(t_sim > 500, tst_level == 1, pats_per_day == 3) %>%
  select(t_sim, tst_level, total_pat_number, pats_per_day, starts_with("total_adult_pat_abun_")) %>%
  mutate(
    repeat_id = ((row_number() - 1) %/% 4500) + 1,
    across(starts_with("total_adult_pat_abun_"), ~ .x / total_pat_number, .names = "average_{.col}")
  ) %>%
  pivot_longer(
    cols = starts_with("average_total_adult_pat_abun_"),
    names_to = "Guild",
    values_to = "Average_abun"
  ) %>%
  mutate(Guild = str_remove(Guild, "average_total_adult_pat_abun_")) %>%
  group_by(t_sim, pats_per_day, Guild) %>%
  summarise(mean_average_abun = mean(Average_abun), .groups = "drop") %>%
  group_by(t_sim, pats_per_day) %>%
  mutate(total_mean_abun = sum(mean_average_abun)) 

## Plot
figure_1 <- ggplot(plot_data_1, 
       aes(x = t_sim,
           y = total_mean_abun)) +
  geom_smooth(method = "loess", span = 0.05, color = "black", size = 1.2) +
  labs(
    x = "Days",
    y = "Mean Adult Abundance Per Pat"
  ) +
  theme_pub_a4

agg_png("figure_1.png", width = 6.27, height = 9.69, units = "in", res = 300)
print(figure_1)
dev.off()

# --- Figure 2 (Guild log abundance over TST) ----------------------------------
## Setup
plot_data_2 <- treatment_results_44_adults_die %>%
  filter(t_sim > 500) %>% # only last 500 days
  select(-c(total_pat_number, average_pat_biomass, iteration_id)) %>%
  select(-starts_with("average_adult_pat_abun_")) %>%
  rename_with(~ paste0(
    case_when(
      str_detect(., "total_adult_pat_abun") ~ "Pat",
      str_detect(., "total_adult_field_abun") ~ "Field",
      str_detect(., "total_larvae_feeding_abun") ~ "Feeding",
      str_detect(., "total_larvae_pupating_abun") ~ "Pupating",
      TRUE ~ .
    ),
    "_",
    str_extract(., "\\d+$")
  ), starts_with("total_")) %>%
  pivot_longer(
    cols = starts_with(c("Pat_", "Field_", "Feeding_", "Pupating_")),
    names_to = "pop_name",
    values_to = "abun"
  ) %>%
  mutate(log_mean_abundance = log10(abun + 1)) %>%
  select(-abun) %>%
  mutate(repeat_id = ((row_number() - 1) %/% 198000) + 1) %>%
  relocate(repeat_id, tst_level, pats_per_day, .before = t_sim) %>%
  separate(pop_name, into = c("Life_stage", "Guild"), sep = "_") %>%
  group_by(tst_level, pats_per_day, Guild, Life_stage) %>%
  summarise(mean_log_abun = mean(log_mean_abundance), .groups = "drop")

## Plot
figure_2_guild_1 <- plot_data_2 %>%
    filter(Guild == 1) %>%
    ggplot(aes(x = tst_level,
               y = mean_log_abun,
               color = as.factor(pats_per_day),
               group = pats_per_day)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~ Life_stage, labeller = labeller(Life_stage = life_stage_labels)) +
    scale_color_viridis_d(name = "Stocking density (AU/ha)",
                          labels = stocking_density_labels) +
    labs(
      x = "TST Level",
      y = "Log10(x + 1) mean abundance",
      title = "A. Fly"
    ) +
    scale_y_continuous(limits = c(0, 6.5)) +
    theme_pub_a4 +
    theme(
      plot.title = element_text(face = "plain", size = 14, hjust = 0),
      axis.title = element_text(face = "bold", size = 12)
    )

figure_2_guild_1

figure_2_guild_2 <- plot_data_2 %>%
  filter(Guild == 2) %>%
  ggplot(aes(x = tst_level,
             y = mean_log_abun,
             color = as.factor(pats_per_day),
             group = pats_per_day)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Life_stage, labeller = labeller(Life_stage = life_stage_labels)) +
  scale_color_viridis_d(name = "Stocking density (AU/ha)",
                        labels = stocking_density_labels) +
  labs(
    x = "TST Level",
    y = "Log10(x + 1) mean abundance",
    title = "B. Endocoprid Beetle"
  ) +
  scale_y_continuous(limits = c(0, 6.5)) +
  theme_pub_a4 +
  theme(
    plot.title = element_text(face = "plain", size = 14, hjust = 0),
    axis.title = element_text(face = "bold", size = 12)
  )

figure_2_guild_2

figure_2_guild_3 <- plot_data_2 %>%
  filter(Guild == 3) %>%
  ggplot(aes(x = tst_level,
             y = mean_log_abun,
             color = as.factor(pats_per_day),
             group = pats_per_day)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Life_stage, labeller = labeller(Life_stage = life_stage_labels)) +
  scale_color_viridis_d(name = "Stocking density (AU/ha)",
                        labels = stocking_density_labels) +
  labs(
    x = "TST Level",
    y = "Log10(x + 1) mean abundance",
    title = "C. Paracoprid beetle"
  ) +
  scale_y_continuous(limits = c(0, 6.5)) +
  theme_pub_a4 +
  theme(
    plot.title = element_text(face = "plain", size = 14, hjust = 0),
    axis.title = element_text(face = "bold", size = 12)
  )

figure_2_guild_3

agg_png("figure_2.png", 
        width = 2480, height = 3508, res = 300)

(figure_2_guild_1 / figure_2_guild_2 / figure_2_guild_3) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

dev.off()

# --- Figure 3 (Faceted abundance plot) ----------------------------------------
## Setup
plot_data_3 <- treatment_results_44_adults_die %>%
  filter(t_sim > 500) %>% # only want last 500 days of data
  select(-c(total_pat_number, average_pat_biomass, iteration_id)) %>% # remove these columns
  select(-starts_with("average_adult_pat_abun_")) %>% # remove these columns
  rename_with(~ paste0( # rename these columns
    case_when(
      str_detect(., "total_adult_pat_abun") ~ "Pat",
      str_detect(., "total_adult_field_abun") ~ "Field",
      str_detect(., "total_larvae_feeding_abun") ~ "Feeding",
      str_detect(., "total_larvae_pupating_abun") ~ "Pupating",
      TRUE ~ .
    ),
    "_",
    str_extract(., "\\d+$")
  ), starts_with("total_")) %>%
  pivot_longer(
    cols = starts_with(c("Pat_", "Field_", "Feeding_", "Pupating_")), # take columns with these names
    names_to = "pop_name", # make a new column with the population names
    values_to = "abun" # name a new column with the abundance values
  ) %>%
  mutate(log_mean_abundance = log10(abun +1)) %>% # new column of logged abundances
  select(-abun) %>% # delete the abun column
  mutate(repeat_id = ((row_number() - 1) %/% 198000) + 1) %>% # reset repeat id to fix mistake
  relocate(repeat_id, tst_level, pats_per_day, .before = t_sim) %>% # organised columns
  separate(pop_name, into = c("Life_stage", "Guild"), sep = "_") %>% # split population names by guild
  filter(Life_stage == "Pat") %>%
  mutate(Guild = factor(Guild, levels = c("1", "2", "3"), labels = c("Fly", "Endocoprid Beetle", "Paracopriod Beetle")))

## Plot
figure_3 <- ggplot(plot_data_3, aes(x = t_sim, y = log_mean_abundance, color = Guild)) +
  geom_line(linewidth = 1) +  # slightly thicker for clarity in print
  ggh4x::facet_grid2(
    tst_level ~ pats_per_day,
    scales = "fixed",
    independent = "none",
    axes = "margins"
  ) +
  labs(
    x = "Time (days)",
    y = "Log10(x + 1) mean abundance",
    color = "Guild"
  ) +
  scale_color_viridis_d(name = "Stocking density (AU/ha)",
                        labels = stocking_density_labels) +
  coord_fixed(ratio = 35) +  # keep consistent scaling
  theme_pub_a4 +
  theme(
    strip.background = element_blank(),      # remove grey facet background
    strip.text = element_text(face = "bold", size = 14), # control facet labels
    legend.position = "bottom",              # move legend if you want under all plots
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12)
  )

agg_png("figure_3.png", 
        width = 2480, height = 3508, res = 300)
print(figure_3)
dev.off()

# Adding top and right axis labels, and removing most of the axis ticks in ppt. 

# --- Figure 4 (Faceted pat biomass plot) --------------------------------------
## Setup

## Plot
figure_4 <- ggplot(plot_data_biomass, aes(x = t_sim, y = mean_pat_biomass)) +
  geom_line(linewidth = 1, color = "black") +
  ggh4x::facet_grid2(
    rows = vars(tst_level),
    cols = vars(pats_per_day_label),
    scales = "fixed",
    independent = "none",
    axes = "margins"
  ) +
  labs(
    x = "Time (days)",
    y = "Mean pat biomass (g)"
  ) +
 theme_pub_a4 +
 theme(
    strip.background = element_blank(),      # remove grey facet background
    strip.text = element_text(face = "bold", size = 14), # control facet labels
    legend.position = "bottom",              # move legend if you want under all plots
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12)
  )

agg_png("figure_4.png", 
        width = 2480, height = 3508, res = 300)
print(figure_4)
dev.off()
