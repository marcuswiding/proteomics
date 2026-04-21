library(readxl)
library(ggplot2)
library(tidyverse)
library(patchwork)
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
  theme_classic(base_size = 12, base_family = "Times New Roman") %+replace%
    theme(
      # Solid black L-shaped axes and ticks
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Black text for axes
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 16, face = "bold"),
      
      # Clean legend with no boxes
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, face = "bold", color = "black"),
      
      # Remove all grid lines and background colors
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Center the plot title
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      
      # Native tag handling
      plot.tag = element_text(size = 18, face = "bold", color = "black", margin = margin(r = 15, b = 15)),
      plot.tag.position = "topleft",
      
      # Increased Top and Left margins so the tags have room to exist without being cut off
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}
# ------------------------------------------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-Chalmers/Masters thesis/Results/96 well plate cultivation")
processresults <- function(dataframe, columns){
  dataframe %>%
    pivot_longer(cols = all_of(columns), 
                 names_to = "Well", 
                 values_to = "Value") %>%
    
    # str_extract pulls just the letters (e.g., "A" from "A1", "A2", "A3")
    mutate(Condition = str_extract(Well, "[A-Za-z]+")) %>%
    
    group_by(time, Condition) %>%
    mutate(Condition = case_when(
      Condition == "A" ~ "Wild Type",
      Condition == "B" ~ "Evolved Wild Type",
      Condition == "C" ~ "ΔGFD2",
      Condition == "D" ~ "ΔYAP1",
      Condition == "E" ~ "ΔSSN2",
      Condition == "F" ~ "ΔSSN8",
      Condition == "G" ~ "ΔSSN2/SSN8",
      Condition == "H" ~ "ΔINM1"))  %>% 
    
    # 4. Calculate Mean and SD for each concentration at each time point
    summarise(
      Mean_Value = mean(Value, na.rm = TRUE),
      SD_Value = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    mutate(Condition = factor(Condition, levels = c(
      "Wild Type",
      "Evolved Wild Type",
      "ΔGFD2",
      "ΔYAP1",
      "ΔSSN2",
      "ΔSSN8",
      "ΔSSN2/SSN8",
      "ΔINM1"))) #Forces order
}

wells <- paste0(rep(LETTERS[1:8], each = 12), 1:12)

colnames <- c("time", wells)

Results <- read_excel('/Users/marcuswiding/Library/CloudStorage/OneDrive-Chalmers/Masters thesis/Results/96 well plate cultivation/Cellusecmarch2026.xlsx', 
                      sheet = 2, col_names = colnames)
Results <- Results %>%
          mutate(time = time/60)

SSH20 <- c(paste0(rep(LETTERS[1:8], each = 3), 1:3))
SSH40 <- c(paste0(rep(LETTERS[1:8], each = 3), 4:6))
SSH60 <- c(paste0(rep(LETTERS[1:8], each = 3), 7:9))
SSH80 <- c(paste0(rep(LETTERS[1:8], each = 3), 10:12))

longresults_SSH20 <- processresults(Results, SSH20)
longresults_SSH40 <- processresults(Results, SSH40)
longresults_SSH60 <- processresults(Results, SSH60)
longresults_SSH80 <- processresults(Results, SSH80)


#Plots created using ggplot
SSH20plot <- ggplot(longresults_SSH20, aes(x = time, y = Mean_Value, color = Condition, fill = Condition)) +
    
    # Add the standard deviation as a transparent shaded ribbon (alpha = 0.2)
    # color = NA removes the harsh outline from the ribbon itself
    geom_ribbon(aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value), 
                alpha = 0.2, color = NA) +
    
    # Add the solid lines for the averages
    geom_line(linewidth = 1) +
  
    scale_color_manual(values = strain_colors) +
    scale_fill_manual(values = strain_colors) +
    
    # Clean up the labels and theme
    labs(title = "Synthetic Spruce Hydrolysate 20%",
         tag = "A",
      x = "Time (hours)",
         y = "Green value (Mean ± SD)",
         color = "Strain",
         fill = "Strain") +
    theme_marcus() 

SSH40plot <- ggplot(longresults_SSH40, aes(x = time, y = Mean_Value, color = Condition, fill = Condition)) +
  
  # Add the standard deviation as a transparent shaded ribbon (alpha = 0.2)
  # color = NA removes the harsh outline from the ribbon itself
  geom_ribbon(aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value), 
              alpha = 0.2, color = NA) +
  
  # Add the solid lines for the averages
  geom_line(linewidth = 1) +
  
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  # Clean up the labels and theme
  labs(title = "Synthetic Spruce Hydrolysate 40%",
       tag = "B",
        x = "Time (hours)",
       y = "Green value (Mean ± SD)",
       color = "Strain",
       fill = "Strain") +
  theme_marcus() 

SSH60plot <- ggplot(longresults_SSH60, aes(x = time, y = Mean_Value, color = Condition, fill = Condition)) +
  
  # Add the standard deviation as a transparent shaded ribbon (alpha = 0.2)
  # color = NA removes the harsh outline from the ribbon itself
  geom_ribbon(aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value), 
              alpha = 0.2, color = NA) +
  
  # Add the solid lines for the averages
  geom_line(linewidth = 1) +
  
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  # Clean up the labels and theme
  labs(title = "Synthetic Spruce Hydrolysate 60%",
       tag = "C",
    x = "Time (hours)",
       y = "Green value (Mean ± SD)",
       color = "Strain",
       fill = "Strain") +
  theme_marcus() 


SSH80plot <- ggplot(longresults_SSH80, aes(x = time, y = Mean_Value, color = Condition, fill = Condition)) +
  
  # Add the standard deviation as a transparent shaded ribbon (alpha = 0.2)
  # color = NA removes the harsh outline from the ribbon itself
  geom_ribbon(aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value), 
              alpha = 0.2, color = NA) +
  
  # Add the solid lines for the averages
  geom_line(linewidth = 1) +
  
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  # Clean up the labels and theme
  labs(title = "Synthetic Spruce Hydrolysate 80%",
       tag = "D",
        x = "Time (hours)",
       y = "Green value (Mean ± SD)",
       color = "Strain",
       fill = "Strain") +
  theme_marcus() 

(SSH20plot + SSH40plot + SSH60plot + SSH80plot) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.position = "right")
  

# Save the plot
ggsave("SSHmultiwell.png", width = 10, height = 6, dpi = 300)
