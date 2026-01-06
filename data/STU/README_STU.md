## COLUMNS
- `'GEOID'`: The 11-digit FIPS code for Census Tracts
- `'week'`: The Monday's date for each week
- `'Total_time'`, `'mean_time_all'`, `'visit_od_all'`: Total social-infrastructure visit time (in minutes), and the number of visits, and the average visiting time. These 3 columns do not have missing values.
- `'mean_time_Art'`, `'visit_od_Art'`: The average visiting time (in minutes) for Art activities and the number of visits. These columns may have missing values.
- `'mean_time_Event'`, `'visit_od_Event'`: Same as above
- `'mean_time_Civic'`, `'visit_od_Civic'`
- `'mean_time_Eating'`, `'visit_od_Eating'`
- `'mean_time_Sports'`, `'visit_od_Sports'`
- `'mean_time_Consume'`, `'visit_od_Consume'`
- `'mean_time_Grocery'`, `'visit_od_Grocery'`
- `'mean_time_Religious'`, `'visit_od_Religious'`
- `'weekly_smoke_days'`: The number of smoke days in each week, ranging from 0-7
- `'weekly_mean_smoke'`: The weekly average PM2.5 level

## Note:
- There are actually visits in `visit_od_all` that are not categorized as one of the 8 activities, so when you add up `visit_od_XXX`, the total number does not match `visit_od_all`. (Pending confirmation from the paper.)
- Some tracts for some weeks are excluded in the original dataset, due to data privacy (e.g., there might be too few visits in that tract during that week).
