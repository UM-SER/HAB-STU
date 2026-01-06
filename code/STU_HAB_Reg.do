clear

//STU
use "/data/TU_tract_Florida_2019_2023_coast.dta", clear
format GEOID %15.0f
duplicates list GEOID week
save "/data/TU_tract_Florida_2019_2023_coast.dta" ,replace

//weather variables
import delimited "/data/Weather/all_variables.csv",clear
rename geoid GEOID
format GEOID %15.0f
gen week_date = date(week, "YMD")
format week_date %td
drop week
rename week_date week
save all_variables.dta,replace

//STU_Algae_reg_data
foreach dist in 1mile 2mile 3mile {

    use "/data/HAB/algae_tract_`dist'.dta", clear

    destring GEOID,replace
	format GEOID %15.0f
	duplicates list GEOID week
	merge 1:1 GEOID week using TU_tract_Florida_2019_2023_coast.dta //merge STU 
	drop if _merge==1
	drop _merge
	
	replace count_event=0 if count_event==.
	replace average_cell=0 if average_cell==.
	replace total_cell=0 if total_cell==.
	merge 1:1 GEOID week using all_variables.dta //merge control variable
	drop if _merge==2
	
	gen year = year(week)
	gen month = month(week)
	gen week1=week(week)
	gen GEOID_str = string(GEOID, "%011.0f")
	gen county = substr(GEOID_str, 3, 3)
	destring county,replace
	drop GEOID_str
	gen ln_visit_od_all=ln(visit_od_all)
	gen ln_Total_time=ln(Total_time)
	gen ln_mean_time_all=ln(mean_time_all)
	gen average_cell_m=average_cell/1000000
	gen total_cell_m=total_cell/1000000
	egen county_year=group(county year)
	egen county_month=group(county month)
	egen month_year=group(month year)
	egen tract_week=group(GEOID week1)
	egen county_week=group(county week1)

    save "/data/algae_reg_`dist'.dta", replace
}


