# Harmful Algae Blooms Decrease Time Use at Social Infrastructure Places

This repository contains replication code and data availability for “Harmful Algae Blooms Decrease Time Use at Social Infrastructure Places”.
Please direct questions to Haoluan Wang at haoluan.wang@miami.edu, who may pass them to the appropriate coauthor(s).

## Repository organization

-	data: contains raw data and processed data used for models.
-	data/HAB: contains raw HAB data and processed 1, 2, and 3-mile-week-tract data.
-	data/STU: contains raw STU data and processed week-tract data.
-	data/Weather: contains raw weather data and processed week-tract level data for temperature, precipitation, humidity, and wind speed.
-	code: contains scripts for analysis, figures, and tables
-	figures: contains figures and tables

## Data

**Harmful algae cell concentration**

We obtained geo-referenced HAB event data for 2019–2023 from the Florida Fish and Wildlife Conservation Commission. Each event records the sampling date, geographic coordinates (latitude and longitude), and the concentration of Karenia brevis (the toxic algal species responsible for red tides) reported in cells per liter. The raw data are available at
https://geodata.myfwc.com/maps/eecd72261da2412192123f5a96c4150c/about.

Data processing involved several steps. First, we removed observations with zero cell counts. Second, because most sampling sites were located in nearshore waters, we restricted our analysis to 35 coastal Florida counties. For each census tract in these counties, we used the centroid to create a 2-mile buffer and then retained all HAB events within that buffer. As robustness checks, we varied the buffer radius to 1 mile and 3 miles. Finally, we aggregated all HAB events to the tract-by-week level and computed both the average and total cell concentration (million per liter). Data are processed using ./code/XXX.R, and this script produces ./data/HAB/XXX_2mile.csv, XXX_1mile.csv, and XXX_3mile.csv.

**Social-infrastructure time use (STU)**

We obtained the STU data from Wang and Guo (2025). To measure the duration and intensity of time spent at social infrastructure places, the authors derived weekly STU patterns from Advan Weekly Patterns data. This dataset, based on 32 million devices across the U.S. (approximately 15% of U.S. phone users), captures weekly visit information for each point of interest. The aggregated data, ranging from census tracts to metropolitan areas, include both the total number of visits and the total visit time (minutes). The raw data are available at https://doi.org/10.1038/s41597-025-05504-9. We processed these data using ./code/XXX.XX, which outputs ./data/STU/XXX.XX.

**Weather variables**

The gridMET dataset from the Climatology Lab provides daily, high-resolution weather data (1/24° arc-degree, approximately 4 km) from 1979 to the present (Abatzoglou, 2013). For this study, we selected key weather variables, including minimum air temperature (°C), maximum air temperature (°C), precipitation (mm), wind speed (m/s), and relative humidity (%). Daily observations were aggregated to the weekly level based on the STU dataset’s week definitions, and mean values for each variable were computed for every census tract included in the analysis. The raw data are available at https://www.climatologylab.org/gridmet.html. We processed these data using ./code/XXX.XX, which outputs ./data/Weather/XXX.XX. 

Finally, we used weekly STU data for 2019–2023 at the census tract level, aligned to the 2019 census tract boundaries, and matched them to the harmful algae cell concentration data and weather data, using ./code/XXX.XX, which outputs ./data/algae_reg_2mile.dta, algae_reg_1mile.dta and algae_reg_3mile.dta.

**Census tract variables**
To further investigate the socio-demographic characteristics associated with the census tracts in this study, we selected several variables from the American Community Survey (ACS) 5-Year-Estimates 2015-2019 (aligned to the 2019 census tract boundaries), including total population, age groups (0-18, 18-34, 35-64, and over 65 years old), race and ethnicity (White, Black, Asian, and Hispanic populations), and median household income. The raw data are available at https://www.socialexplorer.com/.

## Figures

Figure 1 XXXXXX

Figure 2: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure 3: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure 4: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure 5: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure 6: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Table S1: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure S1: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure S2: Created from algae_reg_1mile.dta, using ./code/Stu_Algae.R

Figure S3: Created from algae_reg_3mile.dta, using ./code/Stu_Algae.R

Figure S4: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure S5: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R

Figure S6: Created from algae_reg_2mile.dta, using ./code/Stu_Algae.R
