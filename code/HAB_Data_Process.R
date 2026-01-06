### require ####
rm(list = ls())
library(tidyverse)
library(dplyr) 
library(sf)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(readxl)
library(terra)
library(lubridate)
library(tidyr)
library(haven)
library(geosphere)
##############################################################################
### Import data, modify the path according to the location of your dataset ###
##############################################################################
###       Part 1.Harmful algae cell concentration data processing          ###
##############################################################################
### Step 1
algae_2014_after <- st_read("/data/HAB/Recent_Harmful_Algal_Bloom_(HAB)_Events/Recent_Harmful_Algal_Bloom_(HAB)_Events.shp")
algae_2014_after <- filter(algae_2014_after, COUNT_ != 0)
algae_2014_after <- st_transform(algae_2014_after, crs = 4326)
algae_2014_after$SAMPLE_DAT <- as.Date(algae_2014_after$SAMPLE_DAT, format = "%Y-%m-%d")
algae_2019_2023 <- algae_2014_after[
  algae_2014_after$SAMPLE_DAT >= as.Date("2019-01-01") &
    algae_2014_after$SAMPLE_DAT <= as.Date("2023-12-31"),]

### Step 2
#county_coast
county_coast <- st_read("/data/HAB/county_coast/county_coast.shp") 
county_coast1 <- st_transform(county_coast, 4326)
#tract
tract_fl <- st_read("/data/HAB/tl_2019_12_tract/tl_2019_12_tract.shp") 
tract_fl1 <- st_transform(tract_fl, 4326)
tract_coast <- tract_fl1[tract_fl1$COUNTYFP %in% county_coast1$COUNTYFP, ]
tract_coast1 <- tract_coast[tract_coast$ALAND > 0, ]

#create buffer
crs_fl_albers <- 3086
tract_coast_proj <- st_transform(tract_coast1, crs = crs_fl_albers)
algae_proj <- st_transform(algae_2019_2023, crs = crs_fl_albers)

buffer_dist_2m <- 2 * 1609.344  # 1 mile = 1609.344 米
buffer_dist_1m <- 1 * 1609.344  # 1 mile = 1609.344 米
buffer_dist_3m <- 3 * 1609.344  # 1 mile = 1609.344 米
tract_buffers_2 <- st_buffer(tract_coast_proj, dist = buffer_dist_2m)
tract_buffers_1 <- st_buffer(tract_coast_proj, dist = buffer_dist_1m)
tract_buffers_3 <- st_buffer(tract_coast_proj, dist = buffer_dist_3m)

tract_buffers_2 <- tract_buffers_2[, c("GEOID")]  
tract_buffers_1 <- tract_buffers_1[, c("GEOID")]  
tract_buffers_3 <- tract_buffers_3[, c("GEOID")]  

algae_in_buffer_2 <- st_join(
  x = algae_proj,
  y = tract_buffers_2,
  join = st_within,
  left = FALSE  
)

algae_in_buffer_1 <- st_join(
  x = algae_proj,
  y = tract_buffers_1,
  join = st_within,
  left = FALSE 
)

algae_in_buffer_3 <- st_join(
  x = algae_proj,
  y = tract_buffers_3,
  join = st_within,
  left = FALSE  
)

### Step 3

algae_in_buffer_2 <- algae_in_buffer_2 %>%
  mutate(
    week = floor_date(SAMPLE_DAT, "week", week_start = 1)
  )

#删除2018-12-31
algae_in_buffer_2 <- algae_in_buffer_2 %>%
  filter(week != as.Date("2018-12-31"))

algae_in_buffer_1 <- algae_in_buffer_1 %>%
  mutate(
    week = floor_date(SAMPLE_DAT, "week", week_start = 1)
  )

#删除2018-12-31
algae_in_buffer_1 <- algae_in_buffer_1 %>%
  filter(week != as.Date("2018-12-31"))

algae_in_buffer_3 <- algae_in_buffer_3 %>%
  mutate(
    week = floor_date(SAMPLE_DAT, "week", week_start = 1)
  )

#删除2018-12-31
algae_in_buffer_3 <- algae_in_buffer_3 %>%
  filter(week != as.Date("2018-12-31"))


algae_in_buffer_week2 <- algae_in_buffer_2 %>%
  group_by(GEOID, week) %>%
  summarise(
    count_event = n(),
    average_cell = mean(COUNT_, na.rm = TRUE),
    total_cell = sum(COUNT_, na.rm = TRUE),
    .groups = "drop"
  ) %>% st_drop_geometry()

# 2mile
write_dta(algae_in_buffer_week2, "/data/HAB/algae_tract_2mile.dta")

algae_in_buffer_week1 <- algae_in_buffer_1 %>%
  group_by(GEOID, week) %>%
  summarise(
    count_event = n(),
    average_cell = mean(COUNT_, na.rm = TRUE),
    total_cell = sum(COUNT_, na.rm = TRUE),
    .groups = "drop"
  ) %>% st_drop_geometry()

# 1mile
write_dta(algae_in_buffer_week1, "/data/HAB/algae_tract_1mile.dta")

algae_in_buffer_week3 <- algae_in_buffer_3 %>%
  group_by(GEOID, week) %>%
  summarise(
    count_event = n(),
    average_cell = mean(COUNT_, na.rm = TRUE),
    total_cell = sum(COUNT_, na.rm = TRUE),
    .groups = "drop"
  ) %>% st_drop_geometry()

# 3mile
write_dta(algae_in_buffer_week3, "/data/HAB/algae_tract_3mile.dta")

##############################################################################
###                        Part 2.Merge the STU data                       ###
##############################################################################

options(scipen = 999)

df_2019 <- read_csv("/data/STU/TU_tract_Florida/TU_tract_Florida_2019.csv") 
df_2019 <- df_2019 %>%
  mutate(GEOID = as.character(GEOID)) %>%  
  mutate(GEOID = as.numeric(GEOID))        

df_2020 <- read_csv("/data/STU/TU_tract_Florida/TU_tract_Florida_2020.csv") 
df_2020 <- df_2020 %>%
  mutate(GEOID = as.character(GEOID)) %>%  
  mutate(GEOID = as.numeric(GEOID))       

df_2021 <- read_csv("/data/STU/TU_tract_Florida/TU_tract_Florida_2021.csv") 
df_2021 <- df_2021 %>%
  mutate(GEOID = as.character(GEOID)) %>%  
  mutate(GEOID = as.numeric(GEOID))        

df_2022 <- read_csv("/data/STU/TU_tract_Florida/TU_tract_Florida_2022.csv") 
df_2022 <- df_2022 %>%
  mutate(GEOID = as.character(GEOID)) %>%  
  mutate(GEOID = as.numeric(GEOID))        

df_2023 <- read_csv("/data/STU/TU_tract_Florida/TU_tract_Florida_2023.csv") 
df_2023 <- df_2023 %>%
  mutate(GEOID = as.character(GEOID)) %>%  
  mutate(GEOID = as.numeric(GEOID))        

combined_data <- bind_rows(
  df_2019,
  df_2020,
  df_2021,
  df_2022,
  df_2023
)
combined_data_coast <- combined_data[combined_data$GEOID %in% tract_coast$GEOID, ]
write_dta(combined_data_coast, "/data/TU_tract_Florida_2019_2023_coast.dta")
















