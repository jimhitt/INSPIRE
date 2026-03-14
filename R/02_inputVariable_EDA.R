# This function checks for required packages, installs missing ones,
# updates renv.lock to track them, and loads all packages

check_and_install <- function(packages) {
  installed <- installed.packages()[, "Package"]
  missing <- packages[!packages %in% installed]
  
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    renv::install(missing)
    renv::snapshot()  # Update renv.lock with new packages
  }
  
  # Load all packages
  for (pkg in packages) {
    library(pkg, character.only = TRUE)
  }
}

# Define required packages for this analysis
required_packages <- c("tidyverse", "DBI", "duckdb", "gt", "gtsummary", "scales")

# Check, install if needed, and load
check_and_install(required_packages)

# clean up variables
rm(check_and_install, required_packages)

# Load packages ----------------------------------------------------------
library(tidyverse)
library(DBI)
library(duckdb)
library(gt)
library(gtsummary)
library(scales)

# Connect to data --------------------------------------------------------
# Note: PHI data stored locally, never committed to GitHub

# load data from 01_multiple_operations.R
data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
first_operations = readRDS(file.path(data_dir,"PHI_first_operations.rds"))

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

ops_file = "PHI_labs.csv"
ops_path = file.path(data_dir, ops_file)

# Verify file exists before trying to load
if (!file.exists(ops_path)) {
  stop(sprintf("Data file not found: %s\nPlease ensure PHI data is in the data/ directory", ops_path))
}

# Load operations table
labs <- dbGetQuery(con, sprintf("SELECT * FROM read_csv_auto('%s')", ops_path))

# Close connection (we'll work in memory from here)
dbDisconnect(con, shutdown = TRUE)
rm(con, data_dir, ops_file)

# Prepare output variables ------------------------------------------------

analysis_num = "02"
analysis_name = "preop_labs"
output_dir = file.path("output", sprintf("%s_%s",analysis_num, analysis_name))
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# helper function to create numbered output paths
make_output_path = function(type, num, description, ext) {
  filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
  return(file.path(output_dir, filename))
}
rm(analysis_num, analysis_name)

# 1) Join first_operations on labs ----------------------------------------


# join the first_operations on the labs table, and filter labs that happened within
# 6 months of surgery


# Create model_data with all variables from paper
model_data <- first_operations %>%
  filter(asa < 6) %>%
  select(
    subject_id, 
    op_id, 
    # Outcome variables
    inhosp_death_time, 
    orout_time,
    # Demographics
    age, 
    sex, 
    height, 
    weight, 
    # Surgical characteristics
    emop,           # emergency operation
    department,     # surgical department
    antype,         # anesthesia type
    asa,
    # Timing variables
    anstart_time, 
    anend_time, 
    orin_time
  ) %>%
  mutate(
    # Create outcome: death within 30 days of leaving OR
    inhosp_death_90day = if_else(
      !is.na(inhosp_death_time) & (inhosp_death_time < (orout_time + 90 * 24 * 60)), 
      TRUE, 
      FALSE
    ),
    # Calculate anesthesia duration (in minutes)
    andur = anend_time - anstart_time,
    # Calculate BMI (only where height > 10)
    bmi = if_else(height > 10, weight / (height / 100)^2, NA_real_),
    # Encode sex as numeric (Male = 1, Female = 0) - matching their code
    sex = as.integer(sex == "M")
  )

# Define the 13 labs we need
lab_names <- c('hb', 'platelet', 'wbc', 'aptt', 'ptinr', 'glucose', 
               'bun', 'albumin', 'ast', 'alt', 'creatinine', 'sodium', 'potassium')

# Filter labs to only the ones we need and within 6-month preop window
preop_labs <- first_operations %>%
  select(subject_id, op_id, orin_time, anstart_time, orout_time, inhosp_death_time) %>%
  left_join(
    labs %>% filter(item_name %in% lab_names),
    by = "subject_id",
    relationship = "many-to-many"
  ) %>%
  # Keep only labs from 6 months before to surgery start
  filter(
    chart_time <= orin_time,
    chart_time >= orin_time - 6 * 30 * 24 * 60
  ) %>%
  select(subject_id, op_id, orin_time, item_name, chart_time, value, orout_time, inhosp_death_time)

# Now we have all preop labs (potentially multiple per lab type per operation)
# Can explore distributions, timing, missingness patterns
# Then later we'll collapse to most recent value per lab per operation

cat(sprintf("Preop labs: %d measurements from %d patients (%.1f%% of cohort)\n",
            nrow(preop_labs),
            n_distinct(preop_labs$subject_id),
            100 * n_distinct(preop_labs$subject_id) / n_distinct(first_operations$subject_id)))

# Table 1: Lab coverage by type (before joining)
lab_coverage <- preop_labs %>%
  group_by(item_name) %>%
  summarise(
    n_measurements = n(),
    n_patients = n_distinct(subject_id),
    pct_coverage = 100 * n_patients / n_distinct(first_operations$subject_id),
    measurements_per_patient = n() / n_distinct(subject_id)
  ) %>%
  arrange(desc(pct_coverage))

# Pretty print with gt
table_1_labs = lab_coverage %>%
  gt() %>%
  tab_header(
    title = "Preoperative Lab Coverage",
    subtitle = "6-month lookback window from surgery"
  ) %>%
  cols_label(
    item_name = "Lab Test",
    n_measurements = "Total Measurements",
    n_patients = "Patients with Lab",
    pct_coverage = "Coverage (%)",
    measurements_per_patient = "Measurements per Patient"
  ) %>%
  fmt_number(
    columns = c(n_measurements, n_patients),
    decimals = 0
  ) %>%
  fmt_number(
    columns = pct_coverage,
    decimals = 1
  ) %>%
  fmt_number(
    columns = measurements_per_patient,
    decimals = 2
  )

# save table (for blog)
gtsave(table_1_labs, 
       filename = make_output_path("table",1,"lab_summary","html"),
       inline_css = FALSE)

# Making wide table of preop labs to join on operations -------------------

# Collapse to most recent lab and pivot wide in one step
preop_labs_wide <- preop_labs %>%
  group_by(subject_id, op_id, item_name) %>%
  slice_max(chart_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(subject_id, op_id, item_name, value) %>%
  pivot_wider(
    names_from = item_name,
    values_from = value,
    names_prefix = "preop_"
  )

# left join the widened lab table on the operations model_data table
model_data_labs = model_data %>%
  left_join(preop_labs_wide, by = c('subject_id', 'op_id'))

# Table 2: Descriptive statistics after joining to model_data

# Enhanced lab summary with proper formatting and outlier detection
lab_summary <- tibble(
  lab = str_remove(names(model_data_labs %>% select(starts_with('preop_'))), "^preop_")
) %>%
  mutate(
    n_missing = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~sum(is.na(.))),
    pct_missing = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~100*mean(is.na(.))),
    mean = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~mean(., na.rm=TRUE)),
    sd = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~sd(., na.rm=TRUE)),
    median = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~median(., na.rm=TRUE)),
    q25 = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~quantile(., 0.25, na.rm=TRUE)),
    q75 = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~quantile(., 0.75, na.rm=TRUE)),
    min = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~min(., na.rm=TRUE)),
    max = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~max(., na.rm=TRUE)),
    # Outlier detection using IQR method
    n_outliers = map_dbl(model_data_labs %>% select(starts_with('preop_')), ~{
      q1 <- quantile(., 0.25, na.rm=TRUE)
      q3 <- quantile(., 0.75, na.rm=TRUE)
      iqr <- q3 - q1
      lower_fence <- q1 - 1.5 * iqr
      upper_fence <- q3 + 1.5 * iqr
      sum(. < lower_fence | . > upper_fence, na.rm=TRUE)
    }),
    pct_outliers = 100 * n_outliers / (nrow(model_data_labs) - n_missing)
  ) %>%
  mutate(
    # Create formatted columns for display
    mean_sd = sprintf("%.2f (%.2f)", mean, sd),
    median_iqr = sprintf("%.2f [%.2f, %.2f]", median, q25, q75),
    range = sprintf("%.2f - %.2f", min, max)
  ) %>%
  arrange(pct_missing)

# Create pretty table for blog
table_2_lab_stats <- lab_summary %>%
  select(lab, n_missing, pct_missing, mean_sd, median_iqr, range, n_outliers, pct_outliers) %>%
  gt() %>%
  tab_header(
    title = "Preoperative Laboratory Values - Descriptive Statistics",
    subtitle = sprintf("N = %d operations (first operation per patient, ASA < 6)", nrow(model_data_labs))
  ) %>%
  cols_label(
    lab = "Laboratory Test",
    n_missing = "Missing (N)",
    pct_missing = "Missing (%)",
    mean_sd = "Mean (SD)",
    median_iqr = "Median [Q1, Q3]",
    range = "Range",
    n_outliers = "Outliers (N)",
    pct_outliers = "Outliers (%)"
  ) %>%
  fmt_number(
    columns = c(n_missing, n_outliers),
    decimals = 0
  ) %>%
  fmt_number(
    columns = c(pct_missing, pct_outliers),
    decimals = 1
  ) %>%
  tab_source_note(
    source_note = "Most recent lab value within 6 months prior to surgery. Outliers defined as values beyond Q1 - 1.5×IQR or Q3 + 1.5×IQR."
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFF3CD"),
    locations = cells_body(
      columns = pct_outliers,
      rows = pct_outliers > 5
    )
  )

# Save table (for blog)
gtsave(table_2_lab_stats, 
       filename = make_output_path("table", 2, "lab_descriptive_stats", "html"),
       inline_css = FALSE)


# Lab data plot -----------------------------------------------------------

# Prepare data for plotting
lab_plot_data <- model_data_labs %>%
  select(starts_with('preop_')) %>%
  pivot_longer(everything(), names_to = "lab", values_to = "value") %>%
  mutate(lab = str_remove(lab, "^preop_")) %>%
  # Calculate outlier boundaries for each lab
  group_by(lab) %>%
  mutate(
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_fence = q1 - 1.5 * iqr,
    upper_fence = q3 + 1.5 * iqr,
    is_outlier = value < lower_fence | value > upper_fence
  ) %>%
  ungroup()

# Create faceted distribution plot
# Figure 1: Histograms with smoother density curves
fig_1_lab_distributions <- ggplot(lab_plot_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", alpha = 0.5) +
  geom_density(color = "darkblue", linewidth = 1, adjust = 2) +  # adjust = 2 for smoother curve
  geom_vline(aes(xintercept = q1), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_vline(aes(xintercept = q3), linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~lab, scales = "free", ncol = 3) +
  labs(
    title = "Distribution of Preoperative Laboratory Values",
    subtitle = "Dashed red lines show Q1 and Q3; blue curve shows kernel density estimate",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

ggsave(
  filename = make_output_path("figure", 1, "lab_distributions", "png"),
  plot = fig_1_lab_distributions,
  width = 12,
  height = 10,
  dpi = 300
)


# Figure 2: Normalized (z-score) boxplots to compare outlier patterns
lab_plot_normalized <- model_data_labs %>%
  select(starts_with('preop_')) %>%
  pivot_longer(everything(), names_to = "lab", values_to = "value") %>%
  mutate(lab = str_remove(lab, "^preop_")) %>%
  group_by(lab) %>%
  mutate(
    # Z-score normalization
    z_score = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
    median_val = median(value, na.rm = TRUE)
  ) %>%
  ungroup()

# Add median values as text labels
lab_medians <- lab_plot_normalized %>%
  group_by(lab) %>%
  summarise(median_val = first(median_val)) %>%
  mutate(median_text = sprintf("Median: %.2f", median_val))

# Figure 2: Fixed geom_text with proper x aesthetic
fig_2_lab_boxplots <- ggplot(lab_plot_normalized, aes(x = lab, y = z_score, fill = lab)) +
  geom_boxplot(outlier.color = "red", outlier.alpha = 0.3, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "orange", linewidth = 0.5) +
  geom_text(data = lab_medians, aes(x = lab, label = median_text, y = -5), 
            size = 2.5, hjust = 0) +
  coord_flip(ylim = c(-6, 6)) +
  labs(
    title = "Preoperative Laboratory Values - Normalized Distribution",
    subtitle = "Z-scores (standardized values); red points are outliers; dashed orange lines at ±3 SD",
    x = NULL,
    y = "Z-score (standard deviations from mean)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold")
  )

ggsave(
  filename = make_output_path("figure", 2, "lab_boxplots_normalized", "png"),
  plot = fig_2_lab_boxplots,
  width = 10,
  height = 8,
  dpi = 300
)


# Demographic EDA ---------------------------------------------------------

# Prepare demographic data for summary
demo_data <- model_data_labs %>%
  mutate(
    # Create readable labels for categorical variables
    sex_label = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    emop_label = factor(emop, levels = c(0, 1), labels = c("Elective", "Emergency")),
    asa_label = factor(asa, levels = 1:5, labels = paste("ASA", 1:5)),
    inhosp_death_90day_label = factor(inhosp_death_90day, levels = c(0, 1), 
                                      labels = c("Survived", "Died")),
    # Convert anesthesia duration to hours for interpretability
    andur_hours = andur / 60
  )

# Table 3: Demographic and surgical characteristics
table_3_demographics <- demo_data %>%
  select(age, sex_label, height, weight, bmi, emop_label, andur_hours, 
         department, antype, asa_label, inhosp_death_90day_label) %>%
  tbl_summary(
    label = list(
      age ~ "Age (years)",
      sex_label ~ "Sex",
      height ~ "Height (cm)",
      weight ~ "Weight (kg)",
      bmi ~ "BMI (kg/m²)",
      emop_label ~ "Operation Type",
      andur_hours ~ "Anesthesia Duration (hours)",
      department ~ "Surgical Department",
      antype ~ "Anesthesia Type",
      asa_label ~ "ASA Physical Status",
      inhosp_death_90day_label ~ "90-Day Mortality"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd}), {median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing_text = "Missing",
    type = list(
      age ~ "continuous",
      department ~ "categorical",
      antype ~ "categorical"
    )
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Baseline Characteristics and Surgical Details**") %>%
  bold_labels() %>%
  add_n()

# Convert to gt and save
table_3_demographics_gt <- table_3_demographics %>%
  as_gt() %>%
  tab_source_note(
    source_note = sprintf("N = %d operations (first operation per patient, ASA < 6)", 
                          nrow(model_data_labs))
  )

gtsave(table_3_demographics_gt,
       filename = make_output_path("table", 3, "demographics", "html"),
       inline_css = FALSE)

# Additionally, create a focused summary of just the 5 model input variables
table_4_model_inputs <- demo_data %>%
  select(age, sex_label, emop_label, bmi, andur_hours, asa_label, inhosp_death_90day_label) %>%
  tbl_summary(
    label = list(
      age ~ "Age (years)",
      sex_label ~ "Sex", 
      emop_label ~ "Emergency Operation",
      bmi ~ "BMI (kg/m²)",
      andur_hours ~ "Anesthesia Duration (hours)",
      asa_label ~ "ASA Physical Status",
      inhosp_death_90day_label ~ "90-Day Mortality"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd}), {median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Model Input Variables - Summary Statistics**") %>%
  bold_labels() %>%
  add_n()

table_4_model_inputs_gt <- table_4_model_inputs %>%
  as_gt() %>%
  tab_source_note(
    source_note = "Variables used in INSPIRE mortality prediction model replication"
  )

gtsave(table_4_model_inputs_gt,
       filename = make_output_path("table", 4, "model_inputs", "html"),
       inline_css = FALSE)

# Quick check: Compare our cohort characteristics to the published paper
message("\n=== Quick Validation Against Published Paper ===\n")
message(paste(sprintf("Our cohort size: %d\n", nrow(model_data_labs))),
        sprintf("Published (Table 1): 131,109 operations\n"),
        sprintf("Our mortality rate: %.2f%%\n", 100 * mean(model_data_labs$inhosp_death_90day, na.rm=TRUE)),
        sprintf("Published mortality: 1.21%%\n"),
        sprintf("\nMedian age: %.0f years\n", median(model_data_labs$age, na.rm=TRUE)),
        sprintf("Published median age: 60 years (IQR 45-70)\n"),
        sprintf("\nEmergency operations: %.1f%%\n", 100 * mean(model_data_labs$emop, na.rm=TRUE)),
        sprintf("Published emergency: 9.4%%\n"))

# Save model data with labs ----------------------------------------------------

saveRDS(model_data_labs, 'data/PHI_modelData_v01.rds')


