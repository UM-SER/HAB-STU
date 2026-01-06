### require ####
library(haven)
library(dplyr)
library(fixest)
library(ggplot2)
library(patchwork)
library(ggtext)
library(cowplot)
library(grid)
library(gtable)
library(broom)

##############################################################################
### Import data, modify the path according to the location of your dataset ###
##############################################################################
df <- read_dta("/data/algae_reg_2mile.dta") 
df$GEOID <- format(df$GEOID, scientific = FALSE)

socio <- read.csv("/data/ACS2019_FL_tract.csv")

df_1mile <- read_dta("/data/algae_reg_1mile.dta") 
df_1mile$GEOID <- format(df_1mile$GEOID, scientific = FALSE)

df_3mile <- read_dta("/data/algae_reg_3mile.dta") 
df_3mile$GEOID <- format(df_3mile$GEOID, scientific = FALSE)

##############################################################################
###                   Figure 2 Effects of HAB on time use                  ###
##############################################################################

##############################################################################
###                               Figure 2a                                ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df, 
                        cluster = ~GEOID)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df, 
                        cluster = ~GEOID)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

time_ave_cell5 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df, 
                        cluster = ~GEOID)
summary(time_ave_cell5)

time_ave_cell6 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df, 
                        cluster = ~tract_week)
summary(time_ave_cell6)

time_ave_cell7 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df, 
                        cluster = ~GEOID)
summary(time_ave_cell7)

time_ave_cell8 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df, 
                        cluster = ~tract_week)
summary(time_ave_cell8)

models <- list(
  time_ave_cell1,
  time_ave_cell2,
  time_ave_cell3,
  time_ave_cell4,
  time_ave_cell5,
  time_ave_cell6,
  time_ave_cell7,
  time_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "a"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)

figure2a <- wrap_elements(g2)  
figure2a

##############################################################################
###                               Figure 2b                                ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df, 
                          cluster = ~GEOID)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df, 
                          cluster = ~GEOID)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df, 
                          cluster = ~tract_week)
summary(time_total_cell4)

time_total_cell5 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df, 
                          cluster = ~GEOID)
summary(time_total_cell5)

time_total_cell6 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df, 
                          cluster = ~tract_week)
summary(time_total_cell6)

time_total_cell7 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df, 
                          cluster = ~GEOID)
summary(time_total_cell7)

time_total_cell8 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df, 
                          cluster = ~tract_week)
summary(time_total_cell8)

models <- list(
  time_total_cell1,
  time_total_cell2,
  time_total_cell3,
  time_total_cell4,
  time_total_cell5,
  time_total_cell6,
  time_total_cell7,
  time_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15, 
      "p < 0.05"  = 17,
      "p > 0.05"  = 18 
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),  
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw) 

figure2b <- wrap_elements(g2)   
figure2b

##############################################################################
###                               Figure 2c                                ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df, 
                         cluster = ~GEOID)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df, 
                         cluster = ~GEOID)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

visit_ave_cell5 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df, 
                         cluster = ~GEOID)
summary(visit_ave_cell5)

visit_ave_cell6 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df, 
                         cluster = ~tract_week)
summary(visit_ave_cell6)

visit_ave_cell7 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df, 
                         cluster = ~GEOID)
summary(visit_ave_cell7)

visit_ave_cell8 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df, 
                         cluster = ~tract_week)
summary(visit_ave_cell8)

models <- list(
  visit_ave_cell1,
  visit_ave_cell2,
  visit_ave_cell3,
  visit_ave_cell4,
  visit_ave_cell5,
  visit_ave_cell6,
  visit_ave_cell7,
  visit_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015, 0),  
    breaks = c(0, -0.005, -0.010, -0.015),
    labels = function(x) ifelse(
      x == 0,
      "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  )+
  labs(
    x = NULL,
    y = NULL,
    tag = "c" 
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figure2c <- wrap_elements(g2) 
figure2c

##############################################################################
###                               Figure 2d                                ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df, 
                           cluster = ~GEOID)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df, 
                           cluster = ~GEOID)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

visit_total_cell5 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df, 
                           cluster = ~GEOID)
summary(visit_total_cell5)

visit_total_cell6 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df, 
                           cluster = ~tract_week)
summary(visit_total_cell6)

visit_total_cell7 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df, 
                           cluster = ~GEOID)
summary(visit_total_cell7)

visit_total_cell8 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df, 
                           cluster = ~tract_week)
summary(visit_total_cell8)

models <- list(
  visit_total_cell1,
  visit_total_cell2,
  visit_total_cell3,
  visit_total_cell4,
  visit_total_cell5,
  visit_total_cell6,
  visit_total_cell7,
  visit_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) + 
  labs(
    x = NULL,
    y = NULL,
    tag = "d"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figure2d <- wrap_elements(g2)    
figure2d

figure2_withoutlegend <- (figure2a | figure2b) /
  (figure2c | figure2d) +
  plot_layout(
    widths  = c(1, 1),   
    heights = c(1, 1)   
  )

figure2_withoutlegend

df_shapes <- data.frame(
  x = c(3, 3.6, 4.2, 4.75), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.3, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE 
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figure2 <- figure2_withoutlegend / p_legend +
  plot_layout(heights = c(6,6, 0.5))   

figure2

##############################################################################
###             Figure 3 Heterogeneous effects of HAB by season            ###
##############################################################################
df <- df %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "Spring",  
    month %in% 6:8 ~ "Summer",
    month %in% 9:11 ~ "Fall",  
    month %in% c(12, 1, 2) ~ "Winter",  
    TRUE ~ NA_character_ 
  ))

df_spring <- df %>% 
  filter(season == "Spring")
df_summer <- df %>% 
  filter(season == "Summer")
df_fall <- df %>% 
  filter(season == "Fall")
df_winter <- df %>% 
  filter(season == "Winter")

##############################################################################
###                               Figure 3a                                ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_spring, 
                        cluster = ~tract_week)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_summer, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_fall, 
                        cluster = ~tract_week)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_winter, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

models_season <- list(
  Spring = time_ave_cell1,
  Summer = time_ave_cell2,
  Fall   = time_ave_cell3,
  Winter = time_ave_cell4
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    season   = name,
    estimate = ct["average_cell_m", "Estimate"],
    se       = ct["average_cell_m", "Std. Error"],
    p        = ct["average_cell_m", "Pr(>|t|)"]
  )
}

coef_df_season <- bind_rows(
  lapply(names(models_season), \(nm) extract_coef(models_season[[nm]], nm))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    ),
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
  ) %>%
  arrange(season)

figure3a <- ggplot(coef_df_season, aes(x = season, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = NULL,
    y = "ln(Total visit time)", 
    title = "Average cell concentration",
    tag = "a"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure3a

##############################################################################
###                               Figure 3b                                ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_spring, 
                          cluster = ~tract_week)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_summer, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_fall, 
                          cluster = ~tract_week)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_winter, 
                          cluster = ~tract_week)
summary(time_total_cell4)

models_season <- list(
  Spring = time_total_cell1,
  Summer = time_total_cell2,
  Fall   = time_total_cell3,
  Winter = time_total_cell4
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    season   = name,
    estimate = ct["total_cell_m", "Estimate"],
    se       = ct["total_cell_m", "Std. Error"],
    p        = ct["total_cell_m", "Pr(>|t|)"]
  )
}

coef_df_season <- bind_rows(
  lapply(names(models_season), \(nm) extract_coef(models_season[[nm]], nm))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    ),
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
  ) %>%
  arrange(season)

figure3b <- ggplot(coef_df_season, aes(x = season, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = NULL,
    y = "ln(Total visit time)",       
    title = "Total cell concentration",
    tag = "b"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure3b

##############################################################################
###                               Figure 3c                                ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_spring, 
                         cluster = ~tract_week)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_summer, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_fall, 
                         cluster = ~tract_week)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_winter, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

models_season <- list(
  Spring = visit_ave_cell1,
  Summer = visit_ave_cell2,
  Fall   = visit_ave_cell3,
  Winter = visit_ave_cell4
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    season   = name,
    estimate = ct["average_cell_m", "Estimate"],
    se       = ct["average_cell_m", "Std. Error"],
    p        = ct["average_cell_m", "Pr(>|t|)"]
  )
}

coef_df_season <- bind_rows(
  lapply(names(models_season), \(nm) extract_coef(models_season[[nm]], nm))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    ),
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
  ) %>%
  arrange(season)

figure3c <- ggplot(coef_df_season, aes(x = season, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = NULL,
    y = "ln(Number of visits)", 
    title = "Average cell concentration",
    tag = "c"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure3c

##############################################################################
###                               Figure 3d                                ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_spring, 
                           cluster = ~tract_week)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_summer, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_fall, 
                           cluster = ~tract_week)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_winter, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

models_season <- list(
  Spring = visit_total_cell1,
  Summer = visit_total_cell2,
  Fall   = visit_total_cell3,
  Winter = visit_total_cell4
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    season   = name,
    estimate = ct["total_cell_m", "Estimate"],
    se       = ct["total_cell_m", "Std. Error"],
    p        = ct["total_cell_m", "Pr(>|t|)"]
  )
}

coef_df_season <- bind_rows(
  lapply(names(models_season), \(nm) extract_coef(models_season[[nm]], nm))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    ),
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
  ) %>%
  arrange(season)

figure3d <- ggplot(coef_df_season, aes(x = season, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.010,0.005),
    breaks = c(-0.010,-0.005,0,0.005),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = NULL,
    y = "ln(Number of visits)", 
    title = "Total cell concentration",
    tag = "d"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure3d

figure3_1 <- figure3a + plot_spacer() + figure3b +
  plot_layout(widths = c(1, 0.03, 1)) 

figure3_2 <- figure3c + plot_spacer() + figure3d +
  plot_layout(widths = c(1, 0.03, 1)) 

figure3_withoutlegend <- figure3_1/figure3_2
figure3_withoutlegend

df_shapes <- data.frame(
  x = c(3.3, 4.0),
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.4, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE  
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figure3 <- figure3_withoutlegend / p_legend +
  plot_layout(heights = c(5,5, 0.5))  

figure3

##############################################################################
###      Figure 4 Heterogeneous effects of HAB by concentration level      ###
##############################################################################
df <- df %>%
  mutate(
    background1 = ifelse(average_cell > 0 & average_cell <= 1000, 1, 0),
    very_low1 = ifelse(average_cell > 1000 & average_cell <= 10000, 1, 0),
    low1 = ifelse(average_cell > 10000 & average_cell <= 100000, 1, 0),
    medium1 = ifelse(average_cell > 100000 & average_cell <= 1000000, 1, 0),
    high1 = ifelse(average_cell > 1000000, 1, 0)
  )

time_ave_cell <- feols(ln_Total_time ~ background1 + very_low1 + low1 + medium1 + high1 + air_temperature_avg + precipitation + humidity + wind_speed
                       |  GEOID + week1 + month_year, 
                       data = df, 
                       cluster = ~tract_week)
summary(time_ave_cell)

visit_ave_cell <- feols(ln_visit_od_all ~ background1 + very_low1 + low1 + medium1 + high1 + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df, 
                        cluster = ~tract_week)
summary(visit_ave_cell)

df <- df %>%
  mutate(
    background2 = ifelse(total_cell > 0 & total_cell <= 1000, 1, 0),
    very_low2 = ifelse(total_cell > 1000 & total_cell <= 10000, 1, 0),
    low2 = ifelse(total_cell > 10000 & total_cell <= 100000, 1, 0),
    medium2 = ifelse(total_cell > 100000 & total_cell <= 1000000, 1, 0),
    high2 = ifelse(total_cell > 1000000, 1, 0)
  )

time_total_cell <- feols(ln_Total_time ~ background2 + very_low2 + low2 + medium2 + high2 + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df, 
                         cluster = ~tract_week)
summary(time_total_cell)

visit_total_cell <- feols(ln_visit_od_all ~ background2 + very_low2 + low2 + medium2 + high2 + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df, 
                          cluster = ~tract_week)
summary(visit_total_cell)

extract_coef <- function(model, var_name) {
  ct <- summary(model)$coeftable
  data.frame(
    var_name = var_name,
    estimate = ct[var_name, "Estimate"],
    se = ct[var_name, "Std. Error"],
    p = ct[var_name, "Pr(>|t|)"]
  )
}

coef_time_ave_cell <- bind_rows(
  extract_coef(time_ave_cell, "background1"),
  extract_coef(time_ave_cell, "very_low1"),
  extract_coef(time_ave_cell, "low1"),
  extract_coef(time_ave_cell, "medium1"),
  extract_coef(time_ave_cell, "high1")
)

coef_visit_ave_cell <- bind_rows(
  extract_coef(visit_ave_cell, "background1"),
  extract_coef(visit_ave_cell, "very_low1"),
  extract_coef(visit_ave_cell, "low1"),
  extract_coef(visit_ave_cell, "medium1"),
  extract_coef(visit_ave_cell, "high1")
)

coef_time_total_cell <- bind_rows(
  extract_coef(time_total_cell, "background2"),
  extract_coef(time_total_cell, "very_low2"),
  extract_coef(time_total_cell, "low2"),
  extract_coef(time_total_cell, "medium2"),
  extract_coef(time_total_cell, "high2")
)

coef_visit_total_cell <- bind_rows(
  extract_coef(visit_total_cell, "background2"),
  extract_coef(visit_total_cell, "very_low2"),
  extract_coef(visit_total_cell, "low2"),
  extract_coef(visit_total_cell, "medium2"),
  extract_coef(visit_total_cell, "high2")
)

calculate_confidence_intervals <- function(coef_df) {
  coef_df %>%
    mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      sig_cat = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
        right = TRUE,
        include.lowest = TRUE
      )
    )
}

coef_time_ave_cell <- calculate_confidence_intervals(coef_time_ave_cell)
coef_visit_ave_cell <- calculate_confidence_intervals(coef_visit_ave_cell)
coef_time_total_cell <- calculate_confidence_intervals(coef_time_total_cell)
coef_visit_total_cell <- calculate_confidence_intervals(coef_visit_total_cell)

coef_time_ave_cell$var_name <- factor(coef_time_ave_cell$var_name, 
                                      levels = c("background1","very_low1", "low1", "medium1", "high1"),
                                      labels = c("Background", "Very Low", "Low", "Medium", "High"))

figure4a <- ggplot(coef_time_ave_cell, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(
    x == 0, "0",
    scales::label_number(accuracy = 0.01)(x))) +
  labs(
    x = "Cell concentration group",
    y = "ln(Total visit time)",
    title = "Average cell concentration",
    tag = "a"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure4a

coef_time_total_cell$var_name <- factor(coef_time_total_cell$var_name, 
                                        levels = c("background2","very_low2", "low2", "medium2", "high2"),
                                        labels = c("Background","Very Low", "Low", "Medium", "High"))

figure4b <- ggplot(coef_time_total_cell, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.1,0.1),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x))) +
  labs(
    x = "Cell concentration group",
    y = "ln(Total visit time)",
    title = "Total cell concentration",
    tag = "b"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),  
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure4b

coef_visit_ave_cell$var_name <- factor(coef_visit_ave_cell$var_name, 
                                       levels = c("background1", "very_low1", "low1", "medium1", "high1"),
                                       labels = c("Background","Very Low", "Low", "Medium", "High"))

figure4c <- ggplot(coef_visit_ave_cell, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.10,0.1),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x))) +
  labs(
    x = "Cell concentration group",
    y = "ln(Number of visits)",
    title = "Average cell concentration",
    tag = "c"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure4c

coef_visit_total_cell$var_name <- factor(coef_visit_total_cell$var_name, 
                                         levels = c("background2","very_low2", "low2", "medium2", "high2"),
                                         labels = c("Background","Very Low", "Low", "Medium", "High"))

figure4d <- ggplot(coef_visit_total_cell, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01" = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.1,0.1),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x))) +
  labs(
    x = "Cell concentration group",
    y = "ln(Number of visits)",
    title = "Total cell concentration",
    tag = "d"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),  
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figure4d

figure4_1 <- figure4a + plot_spacer() + figure4b +
  plot_layout(widths = c(1, 0.03, 1))   

figure4_2 <- figure4c + plot_spacer() + figure4d +
  plot_layout(widths = c(1, 0.03, 1)) 

figure4_withoutlegend <- figure4_1/figure4_2
figure4_withoutlegend

df_shapes <- data.frame(
  x = c(2.8, 3.5, 4.1),   
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.1, xend = 1.4, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.5, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figure4 <- figure4_withoutlegend / p_legend +
  plot_layout(heights = c(5, 5, 0.5))   

figure4

##############################################################################
###         Figure 5 Heterogeneous effects of HAB by activity type         ###
##############################################################################
df <- df %>%
  mutate(ln_visit_od_Art = log(visit_od_Art))
df <- df %>%
  mutate(ln_visit_od_Eating = log(visit_od_Eating))
df <- df %>%
  mutate(ln_visit_od_Sports = log(visit_od_Sports))
df <- df %>%
  mutate(ln_visit_od_Consume = log(visit_od_Consume))
df <- df %>%
  mutate(ln_visit_od_Grocery = log(visit_od_Grocery))
df <- df %>%
  mutate(ln_visit_od_Religious = log(visit_od_Religious))

art_ave <- feols(ln_visit_od_Art ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                 |  GEOID + week1 + month_year, 
                 data = df, 
                 cluster = ~tract_week)
summary(art_ave)

sport_ave <- feols(ln_visit_od_Sports ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                   |  GEOID + week1 + month_year, 
                   data = df, 
                   cluster = ~tract_week)
summary(sport_ave)

eating_ave <- feols(ln_visit_od_Eating ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                    |  GEOID + week1 + month_year, 
                    data = df, 
                    cluster = ~tract_week)
summary(eating_ave)

consume_ave <- feols(ln_visit_od_Consume ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                     |  GEOID + week1 + month_year, 
                     data = df, 
                     cluster = ~tract_week)
summary(consume_ave)

grocery_ave <- feols(ln_visit_od_Grocery ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                     |  GEOID + week1 + month_year, 
                     data = df, 
                     cluster = ~tract_week)
summary(grocery_ave)

religious_ave <- feols(ln_visit_od_Religious ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                       |  GEOID + week1 + month_year, 
                       data = df, 
                       cluster = ~tract_week)
summary(religious_ave)

models <- list(
  Art       = art_ave,
  Sports    = sport_ave,
  Dining    = eating_ave,
  Consume   = consume_ave,
  Grocery   = grocery_ave,
  Religious = religious_ave
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    type     = name,
    estimate = ct["average_cell_m", "Estimate"],
    se       = ct["average_cell_m", "Std. Error"],
    p        = ct["average_cell_m", "Pr(>|t|)"]
  )
}

coef_df <- bind_rows(
  lapply(names(models), \(nm) extract_coef(models[[nm]], nm))
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%  
  mutate(type = factor(type, levels = type))   

figure5a <-ggplot(coef_df, aes(x = type, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.02, 0.01),               
    breaks = c(-0.02, -0.01, 0, 0.01),      
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = NULL,
    y = "ln(Number of visits)",
    title = "Average cell concentration",
    tag = "a"  
  ) +
  theme(    
    panel.background = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
    plot.title = element_text(hjust = 0.5, size=12),  
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.01, 0.98)  
  )
figure5a


art_total <- feols(ln_visit_od_Art ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                   |  GEOID + week1 + month_year, 
                   data = df, 
                   cluster = ~tract_week)
summary(art_total)

sport_total <- feols(ln_visit_od_Sports ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                     |  GEOID + week1 + month_year, 
                     data = df, 
                     cluster = ~tract_week)
summary(sport_total)

eating_total <- feols(ln_visit_od_Eating ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                      |  GEOID + week1 + month_year, 
                      data = df, 
                      cluster = ~tract_week)
summary(eating_total)

consume_total <- feols(ln_visit_od_Consume ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                       |  GEOID + week1 + month_year, 
                       data = df, 
                       cluster = ~tract_week)
summary(consume_total)

grocery_total <- feols(ln_visit_od_Grocery ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                       |  GEOID + week1 + month_year, 
                       data = df, 
                       cluster = ~tract_week)
summary(grocery_total)

religious_total <- feols(ln_visit_od_Religious ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df, 
                         cluster = ~tract_week)
summary(religious_total)

models <- list(
  Art       = art_total,
  Sports   = sport_total,
  Dining    = eating_total,
  Consume   = consume_total,
  Grocery   = grocery_total,
  Religious = religious_total
)

extract_coef <- function(model, name) {
  ct <- summary(model)$coeftable
  data.frame(
    type     = name,
    estimate = ct["total_cell_m", "Estimate"],
    se       = ct["total_cell_m", "Std. Error"],
    p        = ct["total_cell_m", "Pr(>|t|)"]
  )
}

coef_df <- bind_rows(
  lapply(names(models), \(nm) extract_coef(models[[nm]], nm))
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%  
  mutate(type = factor(type, levels = type))  

figure5b <- ggplot(coef_df, aes(x = type, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) + 
  labs(
    x = NULL,
    y = "ln(Number of visits)",
    title = "Total cell concentration",
    tag = "b"  
  ) +
  theme(    
    panel.background = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
    plot.title = element_text(hjust = 0.5, size=12),  
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.01, 0.98)  
  )
figure5b

figure5_withoutlegend <- figure5a + plot_spacer() + figure5b +
  plot_layout(widths = c(1, 0.03, 1))   
figure5_withoutlegend

df_shapes <- data.frame(
  x = c(2.8, 3.5, 4.1, 4.7),   
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1, xend = 1.3, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.4, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE  
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1),  
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figure5 <- figure5_withoutlegend / p_legend +
  plot_layout(heights = c(5, 0.5))  

figure5

##############################################################################
###    Figure 6 Heterogeneous effects of HAB by socio-demographic group    ###
##############################################################################
df_socio <- df %>%
  left_join(
    socio %>%
      mutate(GEOID = as.character(GEOID)) %>% 
      select(GEOID, Total_population, Age_0.18, Age_65., White, Black, Median_household_income),
    by = "GEOID"
  )

df_missing <- df %>%
  anti_join(socio %>% mutate(GEOID = as.character(GEOID)), by = "GEOID")
nrow(df_missing)

##############################################################################
###                            Figure 6a & 6e                              ###
##############################################################################
population_ave1 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population > median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_ave1)

population_ave2 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population < median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_ave2)

compare_coef_two_feols <- function(mod1, mod2, var, cluster_var) {
  b1 <- coef(mod1)[var]
  b2 <- coef(mod2)[var]
  v1 <- vcov(mod1)[var, var]
  v2 <- vcov(mod2)[var, var]
  diff_beta <- b1 - b2
  se_diff <- sqrt(v1 + v2) 
  df <- min(mod1$nobs, mod2$nobs) - max(length(coef(mod1)), length(coef(mod2)))
  t_stat <- diff_beta / se_diff
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  data.frame(
    var = var,
    beta1 = b1,
    beta2 = b2,
    diff = diff_beta,
    se_diff = se_diff,
    t_stat = t_stat,
    p_value = p_value,
    df = df
  )
}

res <- compare_coef_two_feols(population_ave1, population_ave2, "average_cell_m")

t1 <- tidy(population_ave1, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(population_ave2, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                            
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

diff_res <- compare_coef_two_feols(population_ave1, population_ave2, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2  

figure6a <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), 
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), 
             size = 4, color = "#C84747") +
  scale_shape_manual(values = c(
    "p < 0.001" = 16,  
    "p < 0.01" = 15, 
    "p < 0.05" = 17,  
    "p > 0.05" = 18 
  )) +  
  scale_y_continuous(
    limits = c(-0.022,0.0065),
    breaks = c(-0.02, -0.01, 0),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.01)(x))
  ) +
  labs(
    x = "Total population",
    y = "ln(Total visit time)",
    shape = NULL,      
    title = NULL,
    tag = "a"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.02, 0.98)  
  ) +
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.2), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max * 1.3, 
           label = label_text, size = 3.5) 

figure6a

population_ave3 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population > median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_ave3)

population_ave4 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population < median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_ave4)

t1 <- tidy(population_ave3, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(population_ave4, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                            
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(population_ave3, population_ave4, "average_cell_m")

diff_res <- compare_coef_two_feols(population_ave3, population_ave4, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2  

figure6e <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16, 
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.02,0.0065),
    breaks = c(-0.02, -0.01, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Total population",
    y = "ln(Number of visits)",
    shape = NULL,     
    title = NULL,
    tag = "e"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max * 1.3, 
           label = label_text, size = 3.5) 

figure6e

##############################################################################
###                            Figure 6b & 6f                              ###
##############################################################################
df_socio <- df_socio %>%
  mutate(senior = Age_65./Total_population)

senior_ave1 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior > median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_ave1)

senior_ave2 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior < median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_ave2)

t1 <- tidy(senior_ave1, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(senior_ave2, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                        
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(senior_ave1, senior_ave2, "average_cell_m")

diff_res <- compare_coef_two_feols(senior_ave1, senior_ave2, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2  

figure6b <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.058,0.0098),
    breaks = c(-0.04, -0.02, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Elderly share",
    y = "ln(Total visit time)",
    shape = NULL,      
    title = NULL,
    tag = "b" 
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.55), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.55), linewidth = 0.5) + 
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max * 1.7, 
           label = label_text, size = 3.5) 

figure6b

senior_ave3 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior > median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_ave3)

senior_ave4 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior < median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_ave4)

t1 <- tidy(senior_ave3, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(senior_ave4, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,              
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(senior_ave3, senior_ave4, "average_cell_m")

diff_res <- compare_coef_two_feols(senior_ave3, senior_ave4, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2 

figure6f <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.04,0.008),
    breaks = c(-0.04, -0.02, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Elderly share",
    y = "ln(Number of visits)",
    shape = NULL,     
    title = NULL,
    tag = "f" 
  ) +
  theme(
    legend.position = "none", 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.45), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.45), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.6, 
           label = label_text, size = 3.5) 

figure6f

##############################################################################
###                            Figure 6c & 6g                              ###
##############################################################################
df_socio <- df_socio %>%
  mutate(black1 = Black/Total_population)

summary(df_socio$black1)

black_ave1 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 > median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_ave1)

black_ave2 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 < median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_ave2)

t1 <- tidy(black_ave1, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(black_ave2, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                    
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(black_ave1, black_ave2, "average_cell_m")

diff_res <- compare_coef_two_feols(black_ave1, black_ave2, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2  

figure6c <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15, 
      "p < 0.05"  = 17, 
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.066,0.0098),
    breaks = c(-0.06, -0.03, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Black share",
    y = "ln(Total visit time)",
    shape = NULL,     
    title = NULL,
    tag = "c"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.6), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.6), linewidth = 0.5) +   
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.8, 
           label = label_text, size = 3.5) 

figure6c



black_ave3 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 > median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_ave3)

black_ave4 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 < median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_ave4)


t1 <- tidy(black_ave3, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(black_ave4, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                      
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(black_ave3, black_ave4, "average_cell_m")

diff_res <- compare_coef_two_feols(black_ave3, black_ave4, "average_cell_m")

diff_beta <- diff_res$diff
p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2 


figure6g <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15, 
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.05,0.0085),
    breaks = c(-0.04, -0.02, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Black share",
    y = "ln(Number of visits)",
    shape = NULL,    
    title = NULL,
    tag = "g"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.6), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.6), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.6, 
           label = label_text, size = 3.5) 

figure6g

##############################################################################
###                            Figure 6d & 6h                              ###
##############################################################################
summary(df_socio$Median_household_income)

income_ave1 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income > median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_ave1)

income_ave2 <- feols(
  ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income < median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_ave2)

t1 <- tidy(income_ave1, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(income_ave2, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,              
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(income_ave1, income_ave2, "average_cell_m")

diff_res <- compare_coef_two_feols(income_ave1, income_ave2, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2  

figure6d <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.022,0.0065),
    breaks = c(-0.02, -0.01, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Median household income",
    y = "ln(Total visit time)",
    shape = NULL,      
    title = NULL,
    tag = "d"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +   
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.3, 
           label = label_text, size = 3.5) 

figure6d

income_ave3 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income > median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_ave3)

income_ave4 <- feols(
  ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income < median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_ave4)

t1 <- tidy(income_ave3, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(income_ave4, conf.int = TRUE) %>%
  filter(term == "average_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,             
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(income_ave3, income_ave4, "average_cell_m")

diff_res <- compare_coef_two_feols(income_ave3, income_ave4, "average_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

label_text <- paste0(format(round(diff_beta, 3), big.mark = ","), " (", p_label, ")")

y_max <- 0.004

x_pos <- 1  
x_neg <- 2

figure6h <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.02,0.0065),
    breaks = c(-0.02, -0.01, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.01)(x)
    )
  ) +
  labs(
    x = "Median household income",
    y = "ln(Numner of visits)",
    shape = NULL,    
    title = NULL,
    tag = "h" 
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.25), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.3, 
           label = label_text, size = 3.5) 

figure6h

figure6_withoutlegend <- (figure6a | figure6b | figure6c | figure6d) /
  (figure6e | figure6f | figure6g | figure6h)
figure6_withoutlegend


df_shapes <- data.frame(
  x = c(2.8, 3.45, 4.05, 4.65), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1, xend = 1.3, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.4, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE    
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1),   
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figure6 <- figure6_withoutlegend / p_legend +
  plot_layout(heights = c(5, 5, 0.5))   

figure6

##############################################################################
###                     Table S1 Summary statistics                        ###
##############################################################################
summary_visit <- df %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "Spring",   
    month %in% 6:8 ~ "Summer",   
    month %in% 9:11 ~ "Fall",    
    month %in% c(12, 1, 2) ~ "Winter",
    TRUE ~ NA_character_  
  )) %>%
  group_by(year, season) %>%
  summarise(mean_y = round(mean(visit_od_all, na.rm = TRUE), 0)) %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))) %>%
  arrange(year, season) 

print(summary_visit)

summary_time <- df %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "Spring",   
    month %in% 6:8 ~ "Summer",   
    month %in% 9:11 ~ "Fall",    
    month %in% c(12, 1, 2) ~ "Winter",
    TRUE ~ NA_character_  
  )) %>%
  group_by(year, season) %>%
  summarise(mean_y = round(mean(Total_time, na.rm = TRUE), 0)) %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))) %>%
  arrange(year, season) 

print(summary_time)

summary_averagecell <- df %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "Spring",  
    month %in% 6:8 ~ "Summer",   
    month %in% 9:11 ~ "Fall",   
    month %in% c(12, 1, 2) ~ "Winter",  
    TRUE ~ NA_character_  
  )) %>%
  group_by(year, season) %>%
  summarise(mean_y = round(mean(average_cell_m, na.rm = TRUE), 3)) %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))) %>%
  arrange(year, season) 
print(summary_averagecell)

summary_totalcell <- df %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "Spring",  
    month %in% 6:8 ~ "Summer",   
    month %in% 9:11 ~ "Fall",   
    month %in% c(12, 1, 2) ~ "Winter",  
    TRUE ~ NA_character_  
  )) %>%
  group_by(year, season) %>%
  summarise(mean_y = round(mean(total_cell_m, na.rm = TRUE), 3)) %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))) %>%
  arrange(year, season) 
print(summary_totalcell)

##############################################################################
###          Figure S1 Effects of HAB using the restricted sample          ###
##############################################################################
GEOID_to_remove <- df %>%
  group_by(GEOID) %>%
  summarise(all_zero = all(count_event == 0)) %>%
  filter(all_zero == TRUE) %>%
  pull(GEOID)

df_restricted <- df %>%
  filter(!GEOID %in% GEOID_to_remove)

num_unique_geoid <- df_restricted %>%
  summarise(unique_geoid_count = n_distinct(GEOID))

num_unique_geoid

##############################################################################
###                               Figure S1a                               ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_restricted, 
                        cluster = ~GEOID)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_restricted, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_restricted, 
                        cluster = ~GEOID)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_restricted, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

time_ave_cell5 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_restricted, 
                        cluster = ~GEOID)
summary(time_ave_cell5)

time_ave_cell6 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_restricted, 
                        cluster = ~tract_week)
summary(time_ave_cell6)

time_ave_cell7 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_restricted, 
                        cluster = ~GEOID)
summary(time_ave_cell7)

time_ave_cell8 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_restricted, 
                        cluster = ~tract_week)
summary(time_ave_cell8)

models <- list(
  time_ave_cell1,
  time_ave_cell2,
  time_ave_cell3,
  time_ave_cell4,
  time_ave_cell5,
  time_ave_cell6,
  time_ave_cell7,
  time_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "a"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)

figureS1a <- wrap_elements(g2)  
figureS1a

##############################################################################
###                               Figure S1b                               ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_restricted, 
                          cluster = ~GEOID)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_restricted, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_restricted, 
                          cluster = ~GEOID)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_restricted, 
                          cluster = ~tract_week)
summary(time_total_cell4)

time_total_cell5 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_restricted, 
                          cluster = ~GEOID)
summary(time_total_cell5)

time_total_cell6 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_restricted, 
                          cluster = ~tract_week)
summary(time_total_cell6)

time_total_cell7 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_restricted, 
                          cluster = ~GEOID)
summary(time_total_cell7)

time_total_cell8 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_restricted, 
                          cluster = ~tract_week)
summary(time_total_cell8)

models <- list(
  time_total_cell1,
  time_total_cell2,
  time_total_cell3,
  time_total_cell4,
  time_total_cell5,
  time_total_cell6,
  time_total_cell7,
  time_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15, 
      "p < 0.05"  = 17,
      "p > 0.05"  = 18 
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),  
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw) 

figureS1b <- wrap_elements(g2)   
figureS1b

##############################################################################
###                               Figure S1c                               ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_restricted, 
                         cluster = ~GEOID)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_restricted, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_restricted, 
                         cluster = ~GEOID)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_restricted, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

visit_ave_cell5 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_restricted, 
                         cluster = ~GEOID)
summary(visit_ave_cell5)

visit_ave_cell6 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_restricted, 
                         cluster = ~tract_week)
summary(visit_ave_cell6)

visit_ave_cell7 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_restricted, 
                         cluster = ~GEOID)
summary(visit_ave_cell7)

visit_ave_cell8 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_restricted, 
                         cluster = ~tract_week)
summary(visit_ave_cell8)

models <- list(
  visit_ave_cell1,
  visit_ave_cell2,
  visit_ave_cell3,
  visit_ave_cell4,
  visit_ave_cell5,
  visit_ave_cell6,
  visit_ave_cell7,
  visit_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015, 0),
    breaks = c(0, -0.005, -0.010),
    labels = function(x) ifelse(
      x == 0,
      "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  )+
  labs(
    x = NULL,
    y = NULL,
    tag = "c" 
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS1c <- wrap_elements(g2) 
figureS1c

##############################################################################
###                               Figure S1d                               ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_restricted, 
                           cluster = ~GEOID)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_restricted, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_restricted, 
                           cluster = ~GEOID)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_restricted, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

visit_total_cell5 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_restricted, 
                           cluster = ~GEOID)
summary(visit_total_cell5)

visit_total_cell6 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_restricted, 
                           cluster = ~tract_week)
summary(visit_total_cell6)

visit_total_cell7 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_restricted, 
                           cluster = ~GEOID)
summary(visit_total_cell7)

visit_total_cell8 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_restricted, 
                           cluster = ~tract_week)
summary(visit_total_cell8)

models <- list(
  visit_total_cell1,
  visit_total_cell2,
  visit_total_cell3,
  visit_total_cell4,
  visit_total_cell5,
  visit_total_cell6,
  visit_total_cell7,
  visit_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) + 
  labs(
    x = NULL,
    y = NULL,
    tag = "d"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS1d <- wrap_elements(g2)    
figureS1d

figureS1_withoutlegend <- (figureS1a | figureS1b) /
  (figureS1c | figureS1d) +
  plot_layout(
    widths  = c(1, 1),   
    heights = c(1, 1)   
  )

figureS1_withoutlegend

df_shapes <- data.frame(
  x = c(3, 3.6, 4.2, 4.75), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.3, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE 
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS1 <- figureS1_withoutlegend / p_legend +
  plot_layout(heights = c(6,6, 0.5))   

figureS1


##############################################################################
###             Figure S2 Effects of HAB using a 1-mile buffer             ###
##############################################################################

##############################################################################
###                               Figure S2a                               ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_1mile, 
                        cluster = ~GEOID)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_1mile, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_1mile, 
                        cluster = ~GEOID)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_1mile, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

time_ave_cell5 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_1mile, 
                        cluster = ~GEOID)
summary(time_ave_cell5)

time_ave_cell6 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_1mile, 
                        cluster = ~tract_week)
summary(time_ave_cell6)

time_ave_cell7 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_1mile, 
                        cluster = ~GEOID)
summary(time_ave_cell7)

time_ave_cell8 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_1mile, 
                        cluster = ~tract_week)
summary(time_ave_cell8)

models <- list(
  time_ave_cell1,
  time_ave_cell2,
  time_ave_cell3,
  time_ave_cell4,
  time_ave_cell5,
  time_ave_cell6,
  time_ave_cell7,
  time_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "a"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)

figureS2a <- wrap_elements(g2)  
figureS2a

##############################################################################
###                               Figure S2b                               ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_1mile, 
                          cluster = ~GEOID)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_1mile, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_1mile, 
                          cluster = ~GEOID)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_1mile, 
                          cluster = ~tract_week)
summary(time_total_cell4)

time_total_cell5 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_1mile, 
                          cluster = ~GEOID)
summary(time_total_cell5)

time_total_cell6 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_1mile, 
                          cluster = ~tract_week)
summary(time_total_cell6)

time_total_cell7 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_1mile, 
                          cluster = ~GEOID)
summary(time_total_cell7)

time_total_cell8 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_1mile, 
                          cluster = ~tract_week)
summary(time_total_cell8)

models <- list(
  time_total_cell1,
  time_total_cell2,
  time_total_cell3,
  time_total_cell4,
  time_total_cell5,
  time_total_cell6,
  time_total_cell7,
  time_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15, 
      "p < 0.05"  = 17,
      "p > 0.05"  = 18 
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0023,0.001),
    breaks = c(-0.002, -0.001, 0, 0.001),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),  
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw) 

figureS2b <- wrap_elements(g2)   
figureS2b

##############################################################################
###                               Figure S2c                               ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_1mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_1mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_1mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_1mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

visit_ave_cell5 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_1mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell5)

visit_ave_cell6 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_1mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell6)

visit_ave_cell7 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_1mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell7)

visit_ave_cell8 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_1mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell8)

models <- list(
  visit_ave_cell1,
  visit_ave_cell2,
  visit_ave_cell3,
  visit_ave_cell4,
  visit_ave_cell5,
  visit_ave_cell6,
  visit_ave_cell7,
  visit_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015, 0),  
    breaks = c(0, -0.005, -0.010, -0.015),
    labels = function(x) ifelse(
      x == 0,
      "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  )+
  labs(
    x = NULL,
    y = NULL,
    tag = "c" 
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS2c <- wrap_elements(g2) 
figureS2c

##############################################################################
###                               Figure S2d                               ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_1mile, 
                           cluster = ~GEOID)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_1mile, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_1mile, 
                           cluster = ~GEOID)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_1mile, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

visit_total_cell5 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_1mile, 
                           cluster = ~GEOID)
summary(visit_total_cell5)

visit_total_cell6 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_1mile, 
                           cluster = ~tract_week)
summary(visit_total_cell6)

visit_total_cell7 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_1mile, 
                           cluster = ~GEOID)
summary(visit_total_cell7)

visit_total_cell8 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_1mile, 
                           cluster = ~tract_week)
summary(visit_total_cell8)

models <- list(
  visit_total_cell1,
  visit_total_cell2,
  visit_total_cell3,
  visit_total_cell4,
  visit_total_cell5,
  visit_total_cell6,
  visit_total_cell7,
  visit_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0023, 0.001),
    breaks = c(-0.002, -0.001, 0, 0.001),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "d"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS2d <- wrap_elements(g2)    
figureS2d

figureS2_withoutlegend <- (figureS2a | figureS2b) /
  (figureS2c | figureS2d) +
  plot_layout(
    widths  = c(1, 1),   
    heights = c(1, 1)   
  )

figureS2_withoutlegend

df_shapes <- data.frame(
  x = c(3, 3.6, 4.2, 4.75), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.3, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE 
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS2 <- figureS2_withoutlegend / p_legend +
  plot_layout(heights = c(6,6, 0.5))   

figureS2

##############################################################################
###             Figure S3 Effects of HAB using a 3-mile buffer             ###
##############################################################################

##############################################################################
###                               Figure S3a                               ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_3mile, 
                        cluster = ~GEOID)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_3mile, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_3mile, 
                        cluster = ~GEOID)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_3mile, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

time_ave_cell5 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_3mile, 
                        cluster = ~GEOID)
summary(time_ave_cell5)

time_ave_cell6 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_3mile, 
                        cluster = ~tract_week)
summary(time_ave_cell6)

time_ave_cell7 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_3mile, 
                        cluster = ~GEOID)
summary(time_ave_cell7)

time_ave_cell8 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_3mile, 
                        cluster = ~tract_week)
summary(time_ave_cell8)

models <- list(
  time_ave_cell1,
  time_ave_cell2,
  time_ave_cell3,
  time_ave_cell4,
  time_ave_cell5,
  time_ave_cell6,
  time_ave_cell7,
  time_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "a"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)

figureS3a <- wrap_elements(g2)  
figureS3a

##############################################################################
###                               Figure S3b                               ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_3mile, 
                          cluster = ~GEOID)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_3mile, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_3mile, 
                          cluster = ~GEOID)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_3mile, 
                          cluster = ~tract_week)
summary(time_total_cell4)

time_total_cell5 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_3mile, 
                          cluster = ~GEOID)
summary(time_total_cell5)

time_total_cell6 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_3mile, 
                          cluster = ~tract_week)
summary(time_total_cell6)

time_total_cell7 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_3mile, 
                          cluster = ~GEOID)
summary(time_total_cell7)

time_total_cell8 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_3mile, 
                          cluster = ~tract_week)
summary(time_total_cell8)

models <- list(
  time_total_cell1,
  time_total_cell2,
  time_total_cell3,
  time_total_cell4,
  time_total_cell5,
  time_total_cell6,
  time_total_cell7,
  time_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15, 
      "p < 0.05"  = 17,
      "p > 0.05"  = 18 
    )
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),  
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw) 

figureS3b <- wrap_elements(g2)   
figureS3b

##############################################################################
###                               Figure S3c                               ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_3mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_3mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_3mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_3mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

visit_ave_cell5 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_3mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell5)

visit_ave_cell6 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_3mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell6)

visit_ave_cell7 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_3mile, 
                         cluster = ~GEOID)
summary(visit_ave_cell7)

visit_ave_cell8 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_3mile, 
                         cluster = ~tract_week)
summary(visit_ave_cell8)

models <- list(
  visit_ave_cell1,
  visit_ave_cell2,
  visit_ave_cell3,
  visit_ave_cell4,
  visit_ave_cell5,
  visit_ave_cell6,
  visit_ave_cell7,
  visit_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015, 0),
    breaks = c(0, -0.005, -0.010, -0.015),
    labels = function(x) ifelse(
      x == 0,
      "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  )+
  labs(
    x = NULL,
    y = NULL,
    tag = "c" 
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS3c <- wrap_elements(g2) 
figureS3c

##############################################################################
###                               Figure S3d                               ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_3mile, 
                           cluster = ~GEOID)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_3mile, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_3mile, 
                           cluster = ~GEOID)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_3mile, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

visit_total_cell5 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_3mile, 
                           cluster = ~GEOID)
summary(visit_total_cell5)

visit_total_cell6 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_3mile, 
                           cluster = ~tract_week)
summary(visit_total_cell6)

visit_total_cell7 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_3mile, 
                           cluster = ~GEOID)
summary(visit_total_cell7)

visit_total_cell8 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_3mile, 
                           cluster = ~tract_week)
summary(visit_total_cell8)

models <- list(
  visit_total_cell1,
  visit_total_cell2,
  visit_total_cell3,
  visit_total_cell4,
  visit_total_cell5,
  visit_total_cell6,
  visit_total_cell7,
  visit_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    breaks = c(-0.0009, -0.0006, -0.0003, 0),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.0001)(x))) + 
  labs(
    x = NULL,
    y = NULL,
    tag = "d"
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS3d <- wrap_elements(g2)    
figureS3d

figureS3_withoutlegend <- (figureS3a | figureS3b) /
  (figureS3c | figureS3d) +
  plot_layout(
    widths  = c(1, 1),   
    heights = c(1, 1)   
  )

figureS3_withoutlegend

df_shapes <- data.frame(
  x = c(3, 3.6, 4.2, 4.75), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.3, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE 
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS3 <- figureS3_withoutlegend / p_legend +
  plot_layout(heights = c(6,6, 0.5))   

figureS3

##############################################################################
###    Figure S4 Effects of HAB using the sample excluding the year 2020   ###
##############################################################################
df_e2020 <- df %>%
  filter(year != 2020)

##############################################################################
###                               Figure S4a                               ###
##############################################################################
time_ave_cell1 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_e2020, 
                        cluster = ~GEOID)
summary(time_ave_cell1)

time_ave_cell2 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + year, 
                        data = df_e2020, 
                        cluster = ~tract_week)
summary(time_ave_cell2)

time_ave_cell3 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_e2020, 
                        cluster = ~GEOID)
summary(time_ave_cell3)

time_ave_cell4 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month_year, 
                        data = df_e2020, 
                        cluster = ~tract_week)
summary(time_ave_cell4)

time_ave_cell5 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_e2020, 
                        cluster = ~GEOID)
summary(time_ave_cell5)

time_ave_cell6 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + month + county_year, 
                        data = df_e2020, 
                        cluster = ~tract_week)
summary(time_ave_cell6)

time_ave_cell7 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_e2020, 
                        cluster = ~GEOID)
summary(time_ave_cell7)

time_ave_cell8 <- feols(ln_Total_time ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                        |  GEOID + week1 + year + county_month, 
                        data = df_e2020, 
                        cluster = ~tract_week)
summary(time_ave_cell8)

models <- list(
  time_ave_cell1,
  time_ave_cell2,
  time_ave_cell3,
  time_ave_cell4,
  time_ave_cell5,
  time_ave_cell6,
  time_ave_cell7,
  time_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +  
  labs(
    x = NULL,
    y = NULL,
    tag = "a"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)

figureS4a <- wrap_elements(g2)  
figureS4a

##############################################################################
###                               Figure S4b                               ###
##############################################################################
time_total_cell1 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_e2020, 
                          cluster = ~GEOID)
summary(time_total_cell1)

time_total_cell2 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + year, 
                          data = df_e2020, 
                          cluster = ~tract_week)
summary(time_total_cell2)

time_total_cell3 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_e2020, 
                          cluster = ~GEOID)
summary(time_total_cell3)

time_total_cell4 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month_year, 
                          data = df_e2020, 
                          cluster = ~tract_week)
summary(time_total_cell4)

time_total_cell5 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_e2020, 
                          cluster = ~GEOID)
summary(time_total_cell5)

time_total_cell6 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + month + county_year, 
                          data = df_e2020, 
                          cluster = ~tract_week)
summary(time_total_cell6)

time_total_cell7 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_e2020, 
                          cluster = ~GEOID)
summary(time_total_cell7)

time_total_cell8 <- feols(ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                          |  GEOID + week1 + year + county_month, 
                          data = df_e2020, 
                          cluster = ~tract_week)
summary(time_total_cell8)

models <- list(
  time_total_cell1,
  time_total_cell2,
  time_total_cell3,
  time_total_cell4,
  time_total_cell5,
  time_total_cell6,
  time_total_cell7,
  time_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15, 
      "p < 0.05"  = 17,
      "p > 0.05"  = 18 
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0023,0.001),
    breaks = c(-0.002, -0.001, 0, 0.001),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),  
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Total visit time)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw) 

figureS4b <- wrap_elements(g2)   
figureS4b

##############################################################################
###                               Figure S4c                               ###
##############################################################################
visit_ave_cell1 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_e2020, 
                         cluster = ~GEOID)
summary(visit_ave_cell1)

visit_ave_cell2 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + year, 
                         data = df_e2020, 
                         cluster = ~tract_week)
summary(visit_ave_cell2)

visit_ave_cell3 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_e2020, 
                         cluster = ~GEOID)
summary(visit_ave_cell3)

visit_ave_cell4 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month_year, 
                         data = df_e2020, 
                         cluster = ~tract_week)
summary(visit_ave_cell4)

visit_ave_cell5 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_e2020, 
                         cluster = ~GEOID)
summary(visit_ave_cell5)

visit_ave_cell6 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + month + county_year, 
                         data = df_e2020, 
                         cluster = ~tract_week)
summary(visit_ave_cell6)

visit_ave_cell7 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_e2020, 
                         cluster = ~GEOID)
summary(visit_ave_cell7)

visit_ave_cell8 <- feols(ln_visit_od_all ~ average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                         |  GEOID + week1 + year + county_month, 
                         data = df_e2020, 
                         cluster = ~tract_week)
summary(visit_ave_cell8)

models <- list(
  visit_ave_cell1,
  visit_ave_cell2,
  visit_ave_cell3,
  visit_ave_cell4,
  visit_ave_cell5,
  visit_ave_cell6,
  visit_ave_cell7,
  visit_ave_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "average_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "average_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015, 0),  
    breaks = c(0, -0.005, -0.010, -0.015),
    labels = function(x) ifelse(
      x == 0,
      "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  )+
  labs(
    x = NULL,
    y = NULL,
    tag = "c" 
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"), 
    plot.tag.position = c(0.3, 0.98) 
  ) 
p_top

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "Average cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS4c <- wrap_elements(g2) 
figureS4c

##############################################################################
###                               Figure S4d                               ###
##############################################################################
visit_total_cell1 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_e2020, 
                           cluster = ~GEOID)
summary(visit_total_cell1)

visit_total_cell2 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + year, 
                           data = df_e2020, 
                           cluster = ~tract_week)
summary(visit_total_cell2)

visit_total_cell3 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_e2020, 
                           cluster = ~GEOID)
summary(visit_total_cell3)

visit_total_cell4 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month_year, 
                           data = df_e2020, 
                           cluster = ~tract_week)
summary(visit_total_cell4)

visit_total_cell5 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_e2020, 
                           cluster = ~GEOID)
summary(visit_total_cell5)

visit_total_cell6 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + month + county_year, 
                           data = df_e2020, 
                           cluster = ~tract_week)
summary(visit_total_cell6)

visit_total_cell7 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_e2020, 
                           cluster = ~GEOID)
summary(visit_total_cell7)

visit_total_cell8 <- feols(ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                           |  GEOID + week1 + year + county_month, 
                           data = df_e2020, 
                           cluster = ~tract_week)
summary(visit_total_cell8)

models <- list(
  visit_total_cell1,
  visit_total_cell2,
  visit_total_cell3,
  visit_total_cell4,
  visit_total_cell5,
  visit_total_cell6,
  visit_total_cell7,
  visit_total_cell8
)
names(models) <- paste0("Model ", 1:8)

get_coef_info <- function(m, var = "total_cell_m") {
  s  <- summary(m)
  ct <- s$coeftable
  est <- ct[var, "Estimate"]
  se  <- ct[var, "Std. Error"]
  p   <- ct[var, "Pr(>|t|)"]
  data.frame(estimate = est, se = se, p = p)
}

coef_df <- bind_rows(
  lapply(seq_along(models), function(i) {
    tmp <- get_coef_info(models[[i]], "total_cell_m")
    tmp$model     <- i
    tmp$model_lab <- names(models)[i]
    tmp
  })
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig_cat = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      right = TRUE,
      include.lowest = TRUE
    )
  )

model_order <- coef_df %>%
  arrange(desc(estimate)) %>%
  pull(model)

coef_df <- coef_df %>%
  arrange(desc(estimate)) %>%
  mutate(
    model_f = factor(model, levels = model_order)
  )

p_top <- ggplot(coef_df, aes(x = model_f, y = estimate)) +
  geom_linerange( aes( x = model_f, ymin = lower, ymax = upper ), linewidth = 0.6, color = "#C84747" ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = sig_cat),
    size = 3,
    color = "#C84747"
  ) +
  scale_shape_manual(
    name   = "Point estimate",
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01"  = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0023, 0.001),
    breaks = c(-0.002, -0.001, 0, 0.001),
    labels = function(x) ifelse(x == 0, "0", scales::label_number(accuracy = 0.001)(x))) +
  labs(
    x = NULL,
    y = NULL,
    tag = "d"  
  ) +
  theme(
    panel.background = element_blank(),   
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),   
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.3, 0.98)  
  ) 

feature_levels <- c(
  "dep_title", "dep_value",
  "indep_title", "indep_value",
  "control_title", "control_value",
  "fe_title", "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_title", "cl_tract", "cl_tractweek"
)

feature_labels <- c(
  dep_title   = "**Dependent variable**",
  dep_value   = "ln(Number of visits)",
  indep_title = "**Independent variable**",
  indep_value = "<span style='color:white'>aso</span>Total cell concentration (million per liter)",
  #indep_value = "asoTotal cell concentration (million per liter)",
  control_title = "**Control variables**",
  control_value = "Weather controls",
  fe_title    = "**Fixed effects**",
  fe_type1    = "Tract, Week, Month, Year",
  fe_type2    = "Tract, Week, Month × Year",
  fe_type3    = "Tract, Week, Month, County × Year",
  fe_type4    = "Tract, Week, Year, County × Month",
  cl_title    = "**Standard error clustering**",
  cl_tract    = "Tract",
  cl_tractweek = "Tract × Week"
)

common_ids_with_points <- c("dep_value", "indep_value", "control_value")
common_df <- expand.grid(
  feature_id = common_ids_with_points,
  model     = model_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fe_ids_list <- list(
  "fe_type1",  # model 1
  "fe_type1",  # model 2
  "fe_type2",  # model 3
  "fe_type2",  # model 4
  "fe_type3",  # model 5
  "fe_type3",  # model 6
  "fe_type4",  # model 7
  "fe_type4"   # model 8
)

fe_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = fe_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

cluster_ids_list <- list(
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek",
  "cl_tract", "cl_tractweek"
)

cluster_df <- bind_rows(
  lapply(1:8, function(m) {
    data.frame(
      feature_id = cluster_ids_list[[m]],
      model      = m,
      stringsAsFactors = FALSE
    )
  })
)

header_ids <- c("dep_title", "indep_title", "control_title",
                "fe_title", "cl_title")

header_df <- data.frame(
  feature_id = header_ids,
  model      = NA_integer_,
  stringsAsFactors = FALSE
)

spec_df <- bind_rows(common_df, fe_df, cluster_df, header_df) %>%
  mutate(
    model_f    = factor(model, levels = model_order),
    feature_id = factor(feature_id, levels = rev(feature_levels))
  )

value_ids_with_grid <- c(
  "dep_value",
  "indep_value",
  "control_value",
  "fe_type1", "fe_type2", "fe_type3", "fe_type4",
  "cl_tract", "cl_tractweek"
)

grid_df <- data.frame(
  y = as.numeric(factor(value_ids_with_grid,
                        levels = rev(feature_levels)))
)

p_bottom <- ggplot() +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    color = "grey95"
  ) +
  geom_point(
    data = subset(spec_df, !is.na(model_f)),
    aes(x = model_f, y = feature_id),
    size = 2, color = "#4596CD"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = feature_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_blank(),
    plot.margin        = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text.y        = element_markdown(lineheight = 1.5)
  )

add_super_panel_border <- function(pw, col = "black", lwd = 2) {
  g <- patchwork:::patchworkGrob(pw)
  panel_ids <- grep("panel", g$layout$name)
  
  top_row    <- min(g$layout$t[panel_ids])
  bottom_row <- max(g$layout$b[panel_ids])
  left_col   <- min(g$layout$l[panel_ids])
  right_col  <- max(g$layout$r[panel_ids])
  
  gtable::gtable_add_grob(
    g,
    grobs = grid::rectGrob(
      gp = grid::gpar(col = col, fill = NA, lwd = lwd)
    ),
    t = top_row, b = bottom_row,
    l = left_col, r = right_col,
    z = Inf, name = "super_panel_border"
  )
}

p_top2 <- p_top + theme(panel.border = element_blank(),
                        plot.margin  = margin(5, 5, 0, 5))
p_bottom2 <- p_bottom + theme(panel.border = element_blank(),
                              plot.margin  = margin(0, 5, 5, 5))

pw <- p_top2 / p_bottom2 +
  plot_layout(heights = c(1.5, 1.8)) &
  theme(plot.margin = margin(-1, 2, 2, 2))

g2 <- add_super_panel_border(pw)  

figureS4d <- wrap_elements(g2)    
figureS4d

figureS4_withoutlegend <- (figureS4a | figureS4b) /
  (figureS4c | figureS4d) +
  plot_layout(
    widths  = c(1, 1),   
    heights = c(1, 1)   
  )

figureS4_withoutlegend

df_shapes <- data.frame(
  x = c(3, 3.6, 4.2, 4.75), 
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1.4, xend = 1.7, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.8, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2.3, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE 
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1), 
    ylim = c(0.8, 1.1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS4 <- figureS4_withoutlegend / p_legend +
  plot_layout(heights = c(6,6, 0.5))   

figureS4

##############################################################################
###              Figure S5 Dynamic effects of HAB on time use              ###
##############################################################################
df <- df %>%
  arrange(GEOID, year, month, week1) %>%  
  group_by(GEOID) %>%  
  mutate(
    lag1_average_cell_m = lag(average_cell_m, 1),  
    lag2_average_cell_m = lag(average_cell_m, 2),  
    lag3_average_cell_m = lag(average_cell_m, 3)  
  ) %>%
  ungroup() 

df <- df %>%
  arrange(GEOID, year, month, week1) %>%  
  group_by(GEOID) %>%  
  mutate(
    lead1_average_cell_m = lead(average_cell_m, 1),  
    lead2_average_cell_m = lead(average_cell_m, 2),  
    lead3_average_cell_m = lead(average_cell_m, 3)  
  ) %>%
  ungroup()  

time_ave_cell_lag_lead <- feols(ln_Total_time ~ lag1_average_cell_m + lag2_average_cell_m + lag3_average_cell_m  +
                                  lead1_average_cell_m + lead2_average_cell_m + lead3_average_cell_m  +
                                  average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                                |  GEOID + week1 + month_year, 
                                data = df, 
                                cluster = ~tract_week)

summary(time_ave_cell_lag_lead)

visit_ave_cell_lag_lead <- feols(ln_visit_od_all ~ lag1_average_cell_m + lag2_average_cell_m + lag3_average_cell_m  +
                                   lead1_average_cell_m + lead2_average_cell_m + lead3_average_cell_m  +
                                   average_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
                                 |  GEOID + week1 + month_year, 
                                 data = df, 
                                 cluster = ~tract_week)

summary(visit_ave_cell_lag_lead)

extract_coef <- function(model, var_name) {
  ct <- summary(model)$coeftable
  data.frame(
    var_name = var_name,
    estimate = ct[var_name, "Estimate"],
    se = ct[var_name, "Std. Error"],
    p = ct[var_name, "Pr(>|t|)"]
  )
}

coef_time_lag_lead <- bind_rows(
  extract_coef(time_ave_cell_lag_lead, "lag3_average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "lag2_average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "lag1_average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "lead1_average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "lead2_average_cell_m"),
  extract_coef(time_ave_cell_lag_lead, "lead3_average_cell_m")
)

coef_visit_lag_lead <- bind_rows(
  extract_coef(visit_ave_cell_lag_lead, "lag3_average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "lag2_average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "lag1_average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "lead1_average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "lead2_average_cell_m"),
  extract_coef(visit_ave_cell_lag_lead, "lead3_average_cell_m")
)


calculate_confidence_intervals <- function(coef_df) {
  coef_df %>%
    mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      sig_cat = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
        right = TRUE,
        include.lowest = TRUE
      )
    )
}

coef_time_lag_lead  <- calculate_confidence_intervals(coef_time_lag_lead)
coef_visit_lag_lead <- calculate_confidence_intervals(coef_visit_lag_lead)

coef_time_lag_lead$var_name <- factor(
  coef_time_lag_lead$var_name,
  levels = c(
    "lag3_average_cell_m",
    "lag2_average_cell_m",
    "lag1_average_cell_m",
    "average_cell_m",
    "lead1_average_cell_m",
    "lead2_average_cell_m",
    "lead3_average_cell_m"
  ),
  labels = c("-3", "-2", "-1", "0", "1", "2", "3")
)

coef_visit_lag_lead$var_name <- coef_time_lag_lead$var_name

figureS5a <- ggplot(coef_time_lag_lead, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01"  = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x))) +
  labs(
    x = "Week",
    y = "ln(Total visit time)",
    title = "Average cell concentration",
    tag = "a"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figureS5a

figureS5c <- ggplot(coef_visit_lag_lead, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01"  = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.015,0.003),
    breaks = c(-0.015, -0.010, -0.005, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x))) +
  labs(
    x = "Week",
    y = "ln(Number of visits)",
    title = "Average cell concentration",
    tag = "c"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figureS5c

df <- df %>%
  arrange(GEOID, year, month, week1) %>%   
  group_by(GEOID) %>%
  mutate(
    lag1_total_cell_m = lag(total_cell_m, 1),
    lag2_total_cell_m = lag(total_cell_m, 2),
    lag3_total_cell_m = lag(total_cell_m, 3)
  ) %>%
  ungroup()

df <- df %>%
  arrange(GEOID, year, month, week1) %>%
  group_by(GEOID) %>%
  mutate(
    lead1_total_cell_m = lead(total_cell_m, 1),
    lead2_total_cell_m = lead(total_cell_m, 2),
    lead3_total_cell_m = lead(total_cell_m, 3)
  ) %>%
  ungroup()

time_total_cell_lag_lead <- feols(
  ln_Total_time ~ lag1_total_cell_m + lag2_total_cell_m + lag3_total_cell_m +
    lead1_total_cell_m + lead2_total_cell_m + lead3_total_cell_m +
    total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df,
  cluster = ~tract_week
)

summary(time_total_cell_lag_lead)

visit_total_cell_lag_lead <- feols(
  ln_visit_od_all ~ lag1_total_cell_m + lag2_total_cell_m + lag3_total_cell_m +
    lead1_total_cell_m + lead2_total_cell_m + lead3_total_cell_m +
    total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df,
  cluster = ~tract_week
)

summary(visit_total_cell_lag_lead)

extract_coef_safe <- function(model, var_name) {
  ct <- summary(model)$coeftable
  if (!(var_name %in% rownames(ct))) {
    return(data.frame(
      var_name = var_name,
      estimate = NA_real_,
      se = NA_real_,
      p = NA_real_
    ))
  }
  data.frame(
    var_name = var_name,
    estimate = ct[var_name, "Estimate"],
    se = ct[var_name, "Std. Error"],
    p = ct[var_name, "Pr(>|t|)"]
  )
}

coef_time_total_lag_lead <- bind_rows(
  extract_coef_safe(time_total_cell_lag_lead, "lag3_total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "lag2_total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "lag1_total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "lead1_total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "lead2_total_cell_m"),
  extract_coef_safe(time_total_cell_lag_lead, "lead3_total_cell_m")
)

coef_visit_total_lag_lead <- bind_rows(
  extract_coef_safe(visit_total_cell_lag_lead, "lag3_total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "lag2_total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "lag1_total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "lead1_total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "lead2_total_cell_m"),
  extract_coef_safe(visit_total_cell_lag_lead, "lead3_total_cell_m")
)

calculate_confidence_intervals <- function(coef_df) {
  coef_df %>%
    mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      sig_cat = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
        right = TRUE,
        include.lowest = TRUE
      )
    )
}

coef_time_total_lag_lead  <- calculate_confidence_intervals(coef_time_total_lag_lead)
coef_visit_total_lag_lead <- calculate_confidence_intervals(coef_visit_total_lag_lead)

coef_time_total_lag_lead$var_name <- factor(
  coef_time_total_lag_lead$var_name,
  levels = c(
    "lag3_total_cell_m",
    "lag2_total_cell_m",
    "lag1_total_cell_m",
    "total_cell_m",
    "lead1_total_cell_m",
    "lead2_total_cell_m",
    "lead3_total_cell_m"
  ),
  labels = c("-3", "-2", "-1", "0", "1", "2", "3")
)

coef_visit_total_lag_lead$var_name <- coef_time_total_lag_lead$var_name

figureS5b <- ggplot(coef_time_total_lag_lead, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01"  = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0026, 0.0002),
    breaks = c(-0.0024, -0.0016, -0.0008, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.0001)(x)
    )
  ) +
  labs(
    x = "Week",
    y = "ln(Total visit time)",
    title = "Total cell concentration",
    tag = "b"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figureS5b

figureS5d <- ggplot(coef_visit_total_lag_lead, aes(x = var_name, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat), size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,
      "p < 0.01"  = 15,
      "p < 0.05"  = 17,
      "p > 0.05"  = 18
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0024, 0.0002),
    breaks = c(-0.0024, -0.0016, -0.0008, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.0001)(x)
    )
  ) +
  labs(
    x = "Week",
    y = "ln(Number of visits)",
    title = "Total cell concentration",
    tag = "d"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

figureS5d

figureS5_1 <- figureS5a + plot_spacer() + figureS5b +
  plot_layout(widths = c(1, 0.03, 1))   

figureS5_2 <- figureS5c + plot_spacer() + figureS5d +
  plot_layout(widths = c(1, 0.03, 1))   

figureS5_withoutlegend <- figureS5_1/figureS5_2

df_shapes <- data.frame(
  x = c(2.6, 3.35, 4, 4.66),  
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 0.8, xend = 1.1, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.2, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 1.7, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE   
  ) +
  
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  
  coord_cartesian(
    xlim = c(0, 6.1),   
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS5 <- figureS5_withoutlegend / p_legend +
  plot_layout(heights = c(5, 5, 0.5))   

figureS5

##############################################################################
###    Figure S6 Heterogeneous effects of HAB by socio-demographic group   ###
##############################################################################

##############################################################################
###                           Figure S6a & S6e                             ###
##############################################################################
population_total1 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population > median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_total1)

population_total2 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population < median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_total2)

compare_coef_two_feols <- function(mod1, mod2, var, cluster_var) {
  b1 <- coef(mod1)[var]
  b2 <- coef(mod2)[var]
  v1 <- vcov(mod1)[var, var]
  v2 <- vcov(mod2)[var, var]
  diff_beta <- b1 - b2
  se_diff <- sqrt(v1 + v2)  
  df <- min(mod1$nobs, mod2$nobs) - max(length(coef(mod1)), length(coef(mod2)))
  t_stat <- diff_beta / se_diff
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  data.frame(
    var = var,
    beta1 = b1,
    beta2 = b2,
    diff = diff_beta,
    se_diff = se_diff,
    t_stat = t_stat,
    p_value = p_value,
    df = df
  )
}

t1 <- tidy(population_total1, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(population_total2, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                         
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(population_total1, population_total2, "total_cell_m")

diff_res <- compare_coef_two_feols(population_total1, population_total2, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.0004

x_pos <- 1  
x_neg <- 2

figureS6a <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18    
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0022,0.00065),
    breaks = c(-0.002, -0.001, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Total population",
    y = "ln(Total visit time)",
    shape = NULL,     
    title = NULL,
    tag = "a" 
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.2), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.2), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.3, 
           label = label_text, size = 3.5) 

figureS6a

population_total3 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population > median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_total3)

population_total4 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Total_population < median(Total_population, na.rm = TRUE)),
  cluster = ~tract_week)
summary(population_total4)

t1 <- tidy(population_total3, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(population_total4, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                           
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(population_total3, population_total4,  "total_cell_m")

diff_res <- compare_coef_two_feols(population_total3, population_total4, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.0006

x_pos <- 1  
x_neg <- 2

figureS6e <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16, 
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.002,0.00085),
    breaks = c(-0.002,-0.001,0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Total population",
    y = "ln(Number of visits)",
    shape = NULL,      
    title = NULL,
    tag = "e"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.15), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.15), linewidth = 0.5) +   
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.2, 
           label = label_text, size = 3.5) 

figureS6e

##############################################################################
###                           Figure S6b & S6f                             ###
##############################################################################
senior_total1 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior > median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_total1)

senior_total2 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior < median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_total2)

t1 <- tidy(senior_total1, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(senior_total2, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                           
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(senior_total1, senior_total2, "total_cell_m")

diff_res <- compare_coef_two_feols(senior_total1, senior_total2, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.003

x_pos <- 1  
x_neg <- 2

figureS6b <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16, 
      "p < 0.01" = 15,  
      "p < 0.05"  = 17, 
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.016,0.0048),
    breaks = c(-0.012, -0.006, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Elderly share",
    y = "ln(Total visit time)",
    shape = NULL,     
    title = NULL,
    tag = "b"  
  ) +
  theme(
    legend.position = "none", 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98) 
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.2), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.2), linewidth = 0.5) +   
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.3, 
           label = label_text, size = 3.5) 

figureS6b

senior_total3 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior > median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_total3)

senior_total4 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(senior < median(senior, na.rm = TRUE)),
  cluster = ~tract_week)
summary(senior_total4)

t1 <- tidy(senior_total3, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(senior_total4, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                            
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(senior_total3, senior_total4, "total_cell_m")

diff_res <- compare_coef_two_feols(senior_total3, senior_total4, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.003

x_pos <- 1  
x_neg <- 2 

figureS6f <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.012,0.0046),
    breaks = c(-0.01, -0.005, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Elderly share",
    y = "ln(Number of visits)",
    shape = NULL,      
    title = NULL,
    tag = "f"  
  ) +
  theme(
    legend.position = "none", 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) + 
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.15), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.15), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.25, 
           label = label_text, size = 3.5) 

figureS6f

##############################################################################
###                           Figure S6c & S6g                             ###
##############################################################################
black_total1 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 > median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_total1)

black_total2 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 < median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_total2)

t1 <- tidy(black_total1, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(black_total2, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                             
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(black_total1, black_total2, "total_cell_m")

diff_res <- compare_coef_two_feols(black_total1, black_total2, "total_cell_m")

diff_beta <- diff_res$diff
p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.003

x_pos <- 1  
x_neg <- 2  

figureS6c <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18    
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0062,0.0039),
    breaks = c(-0.006, -0.003, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Black share",
    y = "ln(Total visit time)",
    shape = NULL,      
    title = NULL,
    tag = "c"  
  ) +
  theme(
    legend.position = "none", 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) + 
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.15, 
           label = label_text, size = 3.5) 

figureS6c

black_total3 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 > median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_total3)

black_total4 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(black1 < median(black1, na.rm = TRUE)),
  cluster = ~tract_week)
summary(black_total4)

t1 <- tidy(black_total3, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(black_total4, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                 
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(black_total3, black_total4, "total_cell_m")

diff_res <- compare_coef_two_feols(black_total3, black_total4, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.003

x_pos <- 1  
x_neg <- 2  


figureS6g <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17, 
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.0055,0.0038),
    breaks = c(-0.004, -0.002, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Black share",
    y = "ln(Number of visits)",
    shape = NULL,     
    title = NULL,
    tag = "g"  
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.08), linewidth = 0.5) + 
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.08), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.12, 
           label = label_text, size = 3.5) 

figureS6g

##############################################################################
###                           Figure S6d & S6h                             ###
##############################################################################
income_total1 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income > median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_total1)

income_total2 <- feols(
  ln_Total_time ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income < median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_total2)

t1 <- tidy(income_total1, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(income_total2, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                             
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(income_total1, income_total2, "total_cell_m")

diff_res <- compare_coef_two_feols(income_total1, income_total2, "total_cell_m")

diff_beta <- diff_res$diff

p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.0015

x_pos <- 1 
x_neg <- 2  

figureS6d <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16,  
      "p < 0.01" = 15,  
      "p < 0.05"  = 17, 
      "p > 0.05"  = 18  
    )
  ) +
  scale_y_continuous(
    limits = c(-0.005, 0.0021),
    breaks = c(-0.004, -0.002, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Median household income",
    y = "ln(Total visit time)",
    shape = NULL,     
    title = NULL,
    tag = "d" 
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) +  
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) +  
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.2, 
           label = label_text, size = 3.5) 


figureS6d

income_total3 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income > median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_total3)

income_total4 <- feols(
  ln_visit_od_all ~ total_cell_m + air_temperature_avg + precipitation + humidity + wind_speed
  | GEOID + week1 + month_year,
  data = df_socio %>% 
    filter(Median_household_income < median(Median_household_income, na.rm = TRUE)),
  cluster = ~tract_week)
summary(income_total4)


t1 <- tidy(income_total3, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "High")

t2 <- tidy(income_total4, conf.int = TRUE) %>%
  filter(term == "total_cell_m") %>%
  mutate(group = "Low")

coef_df <- bind_rows(t1, t2)
coef_df

coef_df <- coef_df %>%
  mutate(
    lower = conf.low,
    upper = conf.high,
    sig_cat = cut(
      p.value,                         
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

res <- compare_coef_two_feols(income_total3, income_total4, "total_cell_m")

diff_res <- compare_coef_two_feols(income_total3, income_total4, "total_cell_m")

diff_beta <- diff_res$diff
p_value <- diff_res$p_value

p_label <- ifelse(p_value < 0.001, "P < 0.001", paste0("P=", sprintf("%.3f", p_value)))

formatted_diff <- ifelse(formatC(abs(diff_beta), format = "f", digits = 3) == "0.000", 
                         formatC(abs(diff_beta), format = "f", digits = 3), 
                         formatC(diff_beta, format = "f", digits = 3))

label_text <- paste0(formatted_diff, " (", p_label, ")")

y_max <- 0.0015

x_pos <- 1  
x_neg <- 2  

figureS6h <- ggplot(coef_df, aes(x = group, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.6, color = "#C84747") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = sig_cat),
             size = 4, color = "#C84747") +
  scale_shape_manual(
    values = c(
      "p < 0.001" = 16, 
      "p < 0.01" = 15,  
      "p < 0.05"  = 17,  
      "p > 0.05"  = 18   
    )
  ) +
  scale_y_continuous(
    limits = c(-0.004,0.002),
    breaks = c(-0.004, -0.002, 0),
    labels = function(x) ifelse(
      x == 0, "0",
      scales::label_number(accuracy = 0.001)(x)
    )
  ) +
  labs(
    x = "Median household income",
    y = "ln(Numner of visits)",
    shape = NULL,      
    title = NULL,
    tag = "h" 
  ) +
  theme(
    legend.position = "none",  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, size=10),  
    plot.tag = element_text(size = 15, face = "bold"),  
    plot.tag.position = c(0.02, 0.98)  
  )+
  geom_segment(x = x_pos, xend = x_neg, y = y_max, yend = y_max, linewidth = 0.45) + 
  geom_segment(x = x_pos, xend = x_pos, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) + 
  geom_segment(x = x_neg, xend = x_neg, y = y_max, yend = y_max - (y_max * 0.1), linewidth = 0.5) +  
  annotate("text", x = mean(c(x_pos, x_neg)), y = y_max*1.15, 
           label = label_text, size = 3.5) 

figureS6h

figureS6_withoutlegend <- (figureS6a | figureS6b | figureS6c | figureS6d) /
  (figureS6e | figureS6f | figureS6g | figureS6h)

df_shapes <- data.frame(
  x = c(2.8, 3.45, 4.05, 4.65),  
  y = 1,
  p = c("italic(p) < 0.001",
        "italic(p) < 0.01",
        "italic(p) < 0.05",
        "italic(p) > 0.05")
)

p_legend <- ggplot() +
  geom_segment(
    aes(x = 1, xend = 1.3, y = 1, yend = 1),
    color = "#C84747",
    linewidth = 0.6
  ) +
  annotate("text", x = 1.4, y = 1, label = "95% CI", hjust = 0) +
  annotate("text", x = 2, y = 1, label = "Point estimate", hjust = 0) +
  geom_point(
    data = df_shapes,
    aes(x = x, y = y, shape = p),
    color = "#C84747",
    size = 3
  ) +
  geom_text(
    data = df_shapes,
    aes(x = x + 0.12, y = y, label = p),
    hjust = 0,
    parse = TRUE  
  ) +
  scale_shape_manual(values = c(
    "italic(p) < 0.001" = 16,
    "italic(p) < 0.01" = 15,
    "italic(p) < 0.05"  = 17,
    "italic(p) > 0.05"  = 18
  )) +
  coord_cartesian(
    xlim = c(0, 6.1),   
    ylim = c(0.8, 1.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  )

figureS6 <- figureS6_withoutlegend / p_legend +
  plot_layout(heights = c(5, 5, 0.5))   

figureS6