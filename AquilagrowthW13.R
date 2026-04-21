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
  "ΔSSN2/SSN8"        = "#7E6148", # NPG Brown
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
folder_w13 <- 'SSHcultivationW13'

files_w13 <- list.files(folder_w13, pattern = "\\.xlsx$", full.names = TRUE)

data_list_w13 <- lapply(files_w13, read_excel,
                          sheet = "Biomass-Data",
                          skip = 4)
names(data_list_w13) <- c('E', 'EE', 'EEE', 'F', 'FF', 'FFF', 'G', 'GG', 'GGG', 'H', 'HH', 'HHH')

alldata_w13 <- bind_rows(lapply(names(data_list_w13), function(n) {
  df <- data_list_w13[[n]]
  tibble(
    time = df[[2]], 
    biomass = df[[3]],
    strain = n
  )
}))

alldata_w13$time <- alldata_w13$time / 60 / 60 #Converts time to hours

replicate_groups <- c(
  'E'   = 'ΔSSN2', 'EE'  = 'ΔSSN2', 'EEE' = 'ΔSSN2',
  'F'   = 'ΔSSN8',     'FF'  = 'ΔSSN8',     'FFF' = 'ΔSSN8',
  'G'   = 'ΔSSN2/SSN8',     'GG'  = 'ΔSSN2/SSN8',     'GGG' = 'ΔSSN2/SSN8',
  'H'   = 'ΔINM1',     'HH'  = 'ΔINM1',     'HHH' = 'ΔINM1'
)

# -----------------------------------------------------------------------------
# Data Cleaning
# -----------------------------------------------------------------------------

# The data is cleaned using a rolling median

#This code visualises what outliers are removed using the rolling median
alldata_w13 %>%
  group_by(strain) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    group       = replicate_groups[strain],          # map to group name
    rolling_med = rollmedian(biomass, k = 31, fill = NA, align = "center"),
    deviation   = abs(biomass - rolling_med) / rolling_med,
    is_outlier  = !is.na(deviation) & deviation > 0.1
  ) %>%
  mutate(group = factor(group, levels = c(
    "ΔSSN2",
    "ΔSSN8",
    "ΔSSN2/SSN8",
    "ΔINM1"
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
ggsave('w13outliers.png', width = 10, height = 6, dpi = 300)


alldata_w13_clean <- alldata_w13 %>%
  group_by(strain) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    rolling_med = rollmedian(biomass, k = 31, fill = NA, align = "center"),
    deviation   = abs(biomass - rolling_med) / rolling_med
  ) %>%
  filter(is.na(deviation) | deviation <= 0.1) %>%  # how much the data can deviate
  select(time, biomass, strain) %>%
  ungroup()

summary_w13 <- alldata_w13_clean %>%
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
    "ΔSSN2",
    "ΔSSN8",
    "ΔSSN2/SSN8",
    "ΔINM1"
 ))) #Forces order

#Sampling lines dataframe
vlines_all <- data.frame(
  xintercept = c(31, 31, 31, 31, 50.25, 50.25, 50.25, 50.25),
  label      = c("Diauxic shift", '', '', '', 'Post-diauxic growth', '', '', ''),
  hjust      = -0.1,
  vjust      = 1.5,
  lcolor      = c(
    strain_colors["ΔSSN2"],
    strain_colors["ΔSSN2/SSN8"],
    strain_colors["ΔSSN8"],
    strain_colors["ΔINM1"],
    strain_colors["ΔSSN2"],
    strain_colors["ΔSSN2/SSN8"],
    strain_colors["ΔSSN8"],
    strain_colors["ΔINM1"]
  ),
  color     = 'black',
  group      = c("ΔSSN2", "ΔSSN2/SSN8", "ΔSSN8", "ΔINM1", "ΔSSN2", "ΔSSN2/SSN8", "ΔSSN8", "ΔINM1")  # set directly
)

# Calculate y_val by filtering group first, then finding closest time
vlines_all <- vlines_all %>%
  rowwise() %>%
  mutate(
    y_val = {
      sub <- summary_w13[summary_w13$group == group, ]        # filter to correct strain
      sub$mean_biomass[which.min(abs(sub$time - xintercept))] # find closest timepoint
    }
  ) %>%
  ungroup()

ggplot(summary_w13, aes(x = time, y = mean_biomass, color = group, fill = group)) +
  
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

ggsave('AquilaW13.png', width = 10, height = 6, dpi = 300)