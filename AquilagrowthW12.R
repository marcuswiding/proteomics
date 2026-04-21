library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(zoo)
# ---------------------------------
# Preamble for correct formatting
# ---------------------------------
strain_colors <- c(
  "Wild Type"         = "#3C5488", # NPG Dark Blue (Often used for controls)
  "Evolved Wild Type" = "#8491B4", # NPG Muted Blue/Grey
  "ΔGFD2"             = "#E64B35", # NPG Red
  "ΔYAP1"             = "#4DBBD5", # NPG Light Blue
  "ΔSSN2"             = "#00A087", # NPG Teal/Green
  "ΔSSN8"             = "#F39B7F", # NPG Peach/Orange
  "ΔSSN8/SSN2"        = "#7E6148", # NPG Brown
  "ΔINM1"             = "#DC0000"  # NPG Dark Red
)
#Nature-style ggplot theme
theme_marcus<- function() {
  theme_classic(base_size = 12, base_family = "serif") %+replace%
    theme(
      # Solid black L-shaped axes and ticks
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Black text for axes
      axis.text = element_text(color = "black", size = 20),
      axis.title = element_text(color = "black", size = 20, face = "bold"),
      
      # Clean legend with no boxes
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 20, color = "black"),
      legend.title = element_text(size = 20, face = "bold", color = "black"),
      
      # Remove all grid lines and background colors
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Center the plot title
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      
      # Small margin around the plot
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}
# ----------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-Chalmers/Masters thesis/Results/SSH cultivation")
folder_w12 <- 'SSHcultivationW12'

files_w12 <- list.files(folder_w12, pattern = "\\.xlsx$", full.names = TRUE)

data_list_w12 <- lapply(files_w12, read_excel,
                          sheet = "Biomass-Data",
                          skip = 4)
names(data_list_w12) <- c('A', 'AA', 'AAA', 'B', 'BB', 'BBB', 'C', 'CC', 'CCC', 'D', 'DD', 'DDD')

alldata_w12 <- bind_rows(lapply(names(data_list_w12), function(n) {
  df <- data_list_w12[[n]]
  tibble(
    time = df[[2]], 
    biomass = df[[3]],
    strain = n
  )
}))

alldata_w12$time <- alldata_w12$time / 60 / 60 #Converts time to hours

alldata_w12 <- alldata_w12 %>%
  filter(!(strain %in% c('B', 'BB', 'BBB') & time < 25))

alldata_w12 <- alldata_w12 %>%
  group_by(strain) %>%
  mutate(
    time = ifelse(strain %in% c('B', 'BB', 'BBB'), 
                  time - min(time, na.rm = TRUE),  # shift E.WT to start at 0
                  time)
  ) %>%
  ungroup()

replicate_groups <- c(
  'A'   = 'Wild Type', 'AA'  = 'Wild Type', 'AAA' = 'Wild Type',
  'B'   = 'Evolved Wild Type',     'BB'  = 'Evolved Wild Type',     'BBB' = 'Evolved Wild Type',
  'C'   = 'ΔGFD2',     'CC'  = 'ΔGFD2',     'CCC' = 'ΔGFD2',
  'D'   = 'ΔYAP1',     'DD'  = 'ΔYAP1',     'DDD' = 'ΔYAP1'
)

# -----------------------------------------------------------------------------
# Data Cleaning
# -----------------------------------------------------------------------------

# The data is cleaned using a rolling median

#This code visualises what outliers are removed using the rolling median
alldata_w12 %>%
  group_by(strain) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    group       = replicate_groups[strain],          # map to group name
    rolling_med = rollmedian(biomass, k = 31, fill = NA, align = "center"),
    deviation   = abs(biomass - rolling_med) / rolling_med,
    is_outlier  = !is.na(deviation) & deviation > 0.1
  ) %>%
  mutate(group = factor(group, levels = c(
    "Wild Type",
    "Evolved Wild Type",
    "ΔGFD2",
    "ΔYAP1"
  ))) %>% #Forces order 
  ggplot(aes(x = time, y = biomass)) +
  geom_point(aes(color = group), size = 1) +                          # color by group
  geom_point(data = ~ filter(.x, is_outlier), color = "red", size = 1) + #outliers get red color
  scale_color_manual(values = strain_colors) +
  labs(x = "Time (hours)",
       y = "Biomass (a.u)",
       color = "Strain") +
  facet_wrap(~strain, scales = "free_y", ncol = 3) +
  theme_marcus()
ggsave('w12outliers.png', width = 10, height = 6, dpi = 300)


alldata_w12_clean <- alldata_w12 %>%
  group_by(strain) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    rolling_med = rollmedian(biomass, k = 31, fill = NA, align = "center"),
    deviation   = abs(biomass - rolling_med) / rolling_med
  ) %>%
  filter(is.na(deviation) | deviation <= 0.1) %>%  # how much the data can deviate
  select(time, biomass, strain) %>%
  ungroup()

summary_w12 <- alldata_w12_clean %>%
  mutate(group = replicate_groups[strain]) %>%
  group_by(strain) %>%
  mutate(timepoint = row_number()) %>%   # align by measurement order, not exact time
  group_by(group, timepoint) %>%
  summarise(
    time         = mean(time, na.rm = TRUE),
    mean_biomass = mean(biomass, na.rm = TRUE),
    sd_biomass   = sd(biomass,   na.rm = TRUE),
    n            = n(),
    .groups = "drop"
  ) %>%
  filter(n == 3) %>%   # remove any timepoint where a replicate was cleaned out
  
  mutate(group = factor(group, levels = c(
    "Wild Type",
    "Evolved Wild Type",
    "ΔGFD2",
    "ΔYAP1"
 ))) #Forces order

#Sampling lines dataframe
vlines_all <- data.frame(
  xintercept = c(28.01, 28.01, 46.14-25, 46.14, 50.38, 50.38, 69.84-25, 69.84),
  label      = c("Diauxic shift", "", "", "Diauxic shift", 'Post-diauxic growth', 'Post-diauxic growth', '', 'Post-diauxic growth'),
  hjust      = -0.1,
  vjust      = 1.5,
  lcolor      = c(
    strain_colors["Wild Type"],
    strain_colors["ΔGFD2"],
    strain_colors["Evolved Wild Type"],
    strain_colors["ΔYAP1"],
    strain_colors["Wild Type"],
    strain_colors["ΔGFD2"],
    strain_colors["Evolved Wild Type"],
    strain_colors["ΔYAP1"]
  ),
  color     = 'black',
  group      = c("Wild Type", "ΔGFD2", "Evolved Wild Type", "ΔYAP1", "Wild Type", "ΔGFD2", "Evolved Wild Type", "ΔYAP1")  # set directly
)

# Calculate y_val by filtering group first, then finding closest time
vlines_all <- vlines_all %>%
  rowwise() %>%
  mutate(
    y_val = {
      sub <- summary_w12[summary_w12$group == group, ]        # filter to correct strain
      sub$mean_biomass[which.min(abs(sub$time - xintercept))] # find closest timepoint
    }
  ) %>%
  ungroup()

ggplot(summary_w12, aes(x = time, y = mean_biomass, color = group, fill = group)) +
  
  geom_ribbon(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), 
              alpha = 0.2, color = NA) +
  
  geom_line(linewidth = 1) +
  
#  geom_vline(data = vlines_all, aes(xintercept = xintercept),
#             color = vlines_all$lcolor,
#             linewidth = 0.3, linetype = "dashed") +
  
  geom_point(data = vlines_all, aes(x = xintercept, y = y_val),
             shape = 25, size = 3,
             color = vlines_all$lcolor,
             fill  = vlines_all$lcolor,
             inherit.aes = FALSE) +
  
  geom_text(data = vlines_all,
            aes(x = xintercept, y = y_val, label = label),
            hjust = 0.5,    # centered above triangle
            vjust = -1,   # just above triangle, increase to move higher
            color = vlines_all$color,
            size = 5, lineheight = 0.9,
            inherit.aes = FALSE) +
  
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  labs(x = "Time (hours)",
       y = "Biomass Mean ± SD (a.u)",
       color = "Strain",
       fill = "Strain") +
  theme_marcus() +
  theme(plot.margin = margin(t = 40, r = 5, b = 5, l = 5)) +  # increased top margin
  coord_cartesian(clip = "off")   # prevents text being clipped at plot edge

ggsave('AquilaW12.png', width = 10, height = 6, dpi = 300)