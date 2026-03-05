library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
#This script automatically takes Aquila data, extracts a.u values and plots graphs with all samples

#Choose folder with excel files. If working directory is same as source file, it is enough to specify folder names (entire path not needed)
folder_delft <- 'Delft cultivation'
folder_SSH <- 'SSH'

files_delft <- list.files(folder_delft, pattern = "\\.xlsx$", full.names = TRUE)
files_SSH <- list.files(folder_SSH, pattern = "\\.xlsx$", full.names = TRUE)

#combines all data by applying read_excel and specifying correct sheet, skips first 4 rows of sheet
data_list_delft <- lapply(files_delft, read_excel,
                    sheet = "Biomass-Data",
                    skip = 4)
data_list_SSH <- lapply(files_SSH, read_excel,
                          sheet = "Biomass-Data",
                          skip = 4)

#Sets names later used in legend to names of excel files
names(data_list_delft) <- tools::file_path_sans_ext(basename(files_delft))
names(data_list_SSH) <- tools::file_path_sans_ext(basename(files_SSH))

#binds rows and defines time and biomass columns
alldata_delft <- bind_rows(lapply(names(data_list_delft), function(n) {
  df <- data_list_delft[[n]]
  tibble(
    time = df[[1]], #column 1 for absolute time, column 2 for relative time
    biomass = df[[3]],
    strain = n
  )
}))
alldata_SSH <- bind_rows(lapply(names(data_list_SSH), function(n) {
  df <- data_list_SSH[[n]]
  tibble(
    time = df[[1]], #column 1 for absolute time, column 2 for relative time
    biomass = df[[3]],
    strain = n
  )
}))

# Define start and end times for the plots
start_time <- as.POSIXct("2026-01-25 16:21:20") 
end_time   <- as.POSIXct("2026-02-28 00:14:58")

#The data is filtered based on times above
alldata_delft_filtered <- alldata_delft %>% 
  filter(time >= start_time & time <= end_time)
alldata_SSH_filtered <- alldata_SSH %>% 
  filter(time >= start_time & time <= end_time)

#Plots created using ggplot
delftplot <- ggplot(alldata_delft_filtered, aes(x = time, y = biomass, color = strain)) +
  geom_line(linewidth = 1) +
  labs(x = "Time", y = "Biomass (a.u.)", color = "Strain") +
  theme_minimal() + theme(legend.position = "bottom") 
SSHplot <- ggplot(alldata_SSH_filtered, aes(x = time, y = biomass, color = strain)) +
  geom_line(linewidth = 1) +
  labs(x = "Time", y = "Biomass (a.u.)", color = "Strain") +
  theme_minimal() + theme(legend.position = "bottom") 


# Placed side by side using package patchwork
delftplot + SSHplot & plot_annotation(tag_levels = 'A') & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

