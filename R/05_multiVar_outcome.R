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

# Table/Figure Helper function --------------------------------------------

# Prepare output variables 
blog_name = "2026-03-05-multiVarOutcomeyModel"
output_dir = file.path('posts', blog_name)
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# helper function to create numbered output paths
make_output_path = function(type, num, description, ext) {
  filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
  return(file.path(output_dir, filename))
}
rm(blog_name)
# --- end prepare output variables


# Connect to data --------------------------------------------------------
# Note: PHI data stored locally, never committed to GitHub

# load data from 01_multiple_operations.R
data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
first_operations = readRDS(file.path(data_dir,"PHI_first_operations.rds"))
# load previous model data (v2)
model_data = readRDS(file.path(data_dir, "PHI_modelData_v02.rds"))

# fix the mortality data
first_operations = first_operations %>%
  mutate(all_death_30d = !is.na(allcause_death_time) & allcause_death_time <= orout_time + (30 * 24 * 60),
         all_death_90d = !is.na(allcause_death_time) & allcause_death_time <= orout_time + (90 * 24 * 60))

model_data = model_data %>%
  left_join(
    first_operations %>% select(subject_id, all_death_30d, all_death_90d),
    by = 'subject_id'
  )
write_rds(model_data, file.path(data_dir, "PHI_modelData_v03.rds"))

# Load Labs ---------------------------------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

ops_file = "PHI_labs.csv"
ops_path = file.path(data_dir, ops_file)
# Verify file exists before trying to load
if (!file.exists(ops_path)) {
  stop(sprintf("Data file not found: %s\nPlease ensure PHI data is in the data/ directory", ops_path))
}
# Load operations table
labs <- dbGetQuery(con, sprintf("SELECT * FROM read_csv_auto('%s')", ops_path))


# Load ward vitals --------------------------------------------------------
ops_file = "PHI_ward_vitals.csv"
ops_path = file.path(data_dir, ops_file)
# Verify file exists before trying to load
if (!file.exists(ops_path)) {
  stop(sprintf("Data file not found: %s\nPlease ensure PHI data is in the data/ directory", ops_path))
}
# Load operations table
ward_vitals <- dbGetQuery(con, sprintf("SELECT * FROM read_csv_auto('%s')", ops_path))

# Close connection (we'll work in memory from here)
dbDisconnect(con, shutdown = TRUE)
rm(con, ops_file)

# Prepare output variables ------------------------------------------------

analysis_num = "05"
analysis_name = "multivar_outcome"
output_dir = file.path("output", sprintf("%s_%s",analysis_num, analysis_name))
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# helper function to create numbered output paths
make_output_path = function(type, num, description, ext) {
  filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
  return(file.path(output_dir, filename))
}
rm(analysis_num, analysis_name)

# 1) Join first_operations on labs ----------------------------------------


# not sure why I am repeating this
if (0) {
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
      allcause_death_time,
      icuin_time,
      icuout_time,
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
      orin_time,
      opstart_time,
      opend_time,
      orout_time,
      anend_time
    ) %>%
    mutate(
      # Create outcome: death within 30 days of leaving OR
      allcause_death_30day = if_else(
        !is.na(allcause_death_time) & (allcause_death_time < (orout_time + 30 * 24 * 60)), 
        TRUE, 
        FALSE
      ),
      # Create outcome: death within 90 days of leaving OR
      allcause_death_90day = if_else(
        !is.na(allcause_death_time) & (allcause_death_time < (orout_time + 90 * 24 * 60)), 
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
  
  
}

# Define the 13 labs we need
lab_names <- c('hb', 'platelet', 'wbc', 'aptt', 'ptinr', 'glucose', 'troponin_i',
               'bun', 'albumin', 'ast', 'alt', 'creatinine', 'sodium', 'potassium')

# Filter labs to only the ones we need and within 6-month preop window
preop_labs <- model_data %>%
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

# make postop labs table
postop_labs <- model_data %>%
  select(subject_id, op_id, orin_time, anstart_time, orout_time, inhosp_death_time) %>%
  left_join(
    labs %>% filter(item_name %in% lab_names),
    by = "subject_id",
    relationship = "many-to-many"
  ) %>%
  # Keep only labs from surgery start to 30 days
  filter(
    chart_time > orin_time,
    chart_time <= orin_time + 30 * 24 * 60
  ) %>%
  select(subject_id, op_id, orin_time, item_name, chart_time, value, orout_time, inhosp_death_time)


# find cases with postop ventilation > 24 hours ---------------------------

postop_vent = model_data %>%
  select(subject_id, op_id, orin_time, anstart_time, orout_time, inhosp_death_time) %>%
  left_join(
    ward_vitals %>% filter(item_name == 'vent'),
    by = "subject_id",
    relationship = "many-to-many"
  ) %>%
  # keep only 7 days of data
  filter(
    chart_time > orout_time,
    chart_time <= orout_time + 7 * 24 * 60
  ) %>%
  select(subject_id, op_id, orin_time, item_name, chart_time, value, orout_time, inhosp_death_time) %>%
  mutate(vent_hour = (chart_time - orout_time) / 60)


postop_vent_subj = postop_vent %>% 
  filter(value == 1) %>%
  group_by(subject_id) %>%
  summarise(min_vent_hr = min(vent_hour),
            max_vent_hr = max(vent_hour),
            n = n()) %>%
  arrange(subject_id) %>%
  mutate(prolonged_vent = max_vent_hr > 24)

# 1.01% of patients in the dataset had prologned ventilation (> 24 hours)
round(nrow(postop_vent_subj %>% filter(prolonged_vent)) / nrow(model_data) * 100,2)

model_data = model_data %>%
  left_join(
    postop_vent_subj %>% select(subject_id, prolonged_vent),
    by = 'subject_id'
  ) %>%
  mutate(prolonged_vent = replace_na(prolonged_vent, FALSE))



# calculate pre_op GFR ----------------------------------------------------

# Most recent preop creatinine within 90 days before surgery
preop_cr <- preop_labs %>%
  mutate(days_before = chart_time * -1 / (24 * 60)) %>%
  filter(item_name == "creatinine",
         days_before >= 0,
         days_before <= 90) %>%
  group_by(subject_id) %>%
  slice_min(days_before, n = 1, with_ties = FALSE) %>%
  select(subject_id, baseline_cr = value, days_before_cr = days_before)

# Join demographics and calculate eGFR
calculate_egfr <- function(cr, age, sex) {
  kappa     <- ifelse(sex == 1, 0.9, 0.7)        # 1=male, 0=female
  alpha     <- ifelse(sex == 1, -0.302, -0.241)
  sex_factor <- ifelse(sex == 1, 1.0, 1.012)
  
  142 * pmin(cr/kappa, 1)^alpha *
    pmax(cr/kappa, 1)^(-1.200) *
    0.9938^age * sex_factor
}

preop_egfr <- model_data %>%
  select(subject_id, age, sex) %>%
  inner_join(preop_cr, by = "subject_id") %>%
  mutate(
    baseline_egfr = calculate_egfr(baseline_cr, age, sex),
    ckd_stage = case_when(
      baseline_egfr >= 90 ~ 1L,
      baseline_egfr >= 60 ~ 2L,
      baseline_egfr >= 45 ~ 3L,
      baseline_egfr >= 30 ~ 4L,
      baseline_egfr >= 15 ~ 5L,
      TRUE               ~ 6L
    )
  )

# join on model_data table
model_data = model_data %>%
  left_join(
    preop_egfr %>% select(-sex, -age),
    by = 'subject_id'
  )


# Calculate AKI -----------------------------------------------------------

# Peak postop creatinine within 7 days
postop_cr <- postop_labs %>%
  filter(item_name == "creatinine",
         chart_time >= 0,
         chart_time <= 7 * 24 * 60) %>%        # 7 days in minutes
  group_by(subject_id) %>%
  summarise(
    peak_cr = max(value),
    n_cr_postop = n(),
    .groups = "drop"
  )

# Join baseline and calculate KDIGO staging
aki_table <- preop_egfr %>%
  select(subject_id, baseline_cr) %>%
  inner_join(postop_cr, by = "subject_id") %>%
  mutate(
    delta_cr   = peak_cr - baseline_cr,
    ratio_cr   = peak_cr / baseline_cr,
    AKI_stage1 = delta_cr >= 0.3 | ratio_cr >= 1.5,
    AKI_stage2 = ratio_cr >= 2.0 & ratio_cr < 3.0,
    AKI_stage3 = ratio_cr >= 3.0 | peak_cr >= 4.0,
    AKI_any = pmax(AKI_stage1, AKI_stage2, AKI_stage3, na.rm = TRUE) == 1,
    AKI_stage  = case_when(
      AKI_stage3 ~ 3L,
      AKI_stage2 ~ 2L,
      AKI_stage1 ~ 1L,
      TRUE       ~ 0L
    )
  )

# join on model_data
model_data = model_data %>%
  left_join(
    aki_table %>% select(subject_id, delta_cr, AKI_stage1, AKI_stage2, AKI_stage3, AKI_any, AKI_stage),
    by = 'subject_id'
  )


# Look for acute cardiac event (troponin) ---------------------------------

# Peak troponin within 7 days postop
troponin_table <- postop_labs %>%
  filter(item_name == "troponin_i",
         chart_time >= 0,
         chart_time <= 10080) %>%
  group_by(subject_id) %>%
  summarise(
    peak_troponin = max(value, na.rm = TRUE),
    n_troponin    = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  filter(n_troponin > 0) %>%        # drop subjects with only NA values
  mutate(
    troponin_measured = TRUE,
    pmi_event = peak_troponin > 0.04
  )

# Left join — untested patients are PMI-negative
pmi_final <- model_data %>%
  select(subject_id) %>%
  left_join(troponin_table, by = "subject_id") %>%
  mutate(
    troponin_measured = replace_na(troponin_measured, FALSE),
    pmi_event         = replace_na(pmi_event, FALSE)
  )
  
model_data = model_data %>%
  left_join(
    pmi_final %>% select(subject_id, troponin_measured, pmi_event),
    by = 'subject_id'
  )


# Check the composite outcomes --------------------------------------------

# Adding an ICU admission column, but the max ICU admission time is 24 hours exactly
# So this does not seem to capture unplanned ICU admissions but only planned ICU admissions
# There are some negative numbers, so some patients in the ICU were taken to surgery
# for their first entry into the surgical outcome table.
# Because of this, ICU admission will not be included as an adverse event
# Many of these cases would have been planned ICU admissions.

first_operations = first_operations %>%
  mutate(icu_admission = !is.na(icuin_time) & (icuin_time - orout_time) >= 0 & (icuin_time - orout_time) < 7 * 24 * 60)

model_data = model_data %>%
  left_join(
    first_operations %>% select(subject_id, icu_admission),
    by = 'subject_id'
  )

model_data = model_data %>%
  mutate(
    any_adverse_event = all_death_30d | 
      AKI_any | pmi_event | prolonged_vent
  )

model_data %>%
  summarise(
    death_30d    = sum(all_death_30d, na.rm = TRUE),
    #icu          = sum(icu_admission, na.rm = TRUE),
    aki          = sum(AKI_any, na.rm = TRUE),
    pmi          = sum(pmi_event),
    prol_vent    = sum(prolonged_vent),
    any_adverse  = sum(any_adverse_event, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "outcome", values_to = "n") %>%
  mutate(pct = round(n / nrow(model_data) * 100, 2))

# some missing subjects with 30 day mortality from first_op to model data
# Patients with ASA class 6 (brain-dead organ donors) and those with missing ASA classification were excluded from the analytic cohort.
# that explains the slight difference in mortality statistics from the first_operations table and model_data


# ----- Outcome Incidence Table -----
outcome_tbl <- tibble(
  Outcome = c(
    "30-day mortality",
    "Acute kidney injury (KDIGO any stage, 7d)",
    "Perioperative myocardial injury (troponin >0.04 ng/mL, 7d)",
    "Prolonged mechanical ventilation (>24h)",
    "Any adverse event (composite)"
  ),
  N_events = c(
    sum(model_data$all_death_30d, na.rm = TRUE),
    sum(model_data$AKI_any, na.rm = TRUE),
    sum(model_data$pmi_event, na.rm = TRUE),
    sum(model_data$prolonged_vent, na.rm = TRUE),
    sum(model_data$any_adverse_event, na.rm = TRUE)
  ),
  N_total = nrow(model_data)
) %>%
  mutate(Pct = round(N_events / N_total * 100, 2))

outcome_gt <- outcome_tbl %>%
  gt() %>%
  tab_header(
    title = "Postoperative Adverse Event Incidence",
    subtitle = "First-operation cohort, INSPIRE dataset"
  ) %>%
  cols_label(
    Outcome  = "Outcome",
    N_events = "Events (N)",
    N_total  = "Cohort (N)",
    Pct      = "Incidence (%)"
  ) %>%
  fmt_integer(columns = c(N_events, N_total)) %>%
  fmt_number(columns = Pct, decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Outcome == "Any adverse event (composite)")
  ) %>%
  tab_footnote(
    footnote = "AKI defined using creatinine-based KDIGO criteria; patients without postoperative creatinine measurement assumed AKI-negative (N=63,829, 65.6%).",
    locations = cells_body(columns = Outcome,
                           rows = Outcome == "Acute kidney injury (KDIGO any stage, 7d)")
  ) %>%
  tab_footnote(
    footnote = "PMI defined as peak troponin I >0.04 ng/mL within 7 days; patients without troponin measurement assumed PMI-negative (N=92,911, 95.5%).",
    locations = cells_body(columns = Outcome,
                           rows = Outcome == "Perioperative myocardial injury (troponin >0.04 ng/mL, 7d)")
  ) %>%
  opt_row_striping() %>%
  opt_table_outline() %>%
  tab_options(
    table.width = pct(100),
    data_row.padding = px(6),
    column_labels.font.weight = "600",
    table.font.size = px(14)
  )

gtsave(outcome_gt, make_output_path("Table", 1, "OutcomeIncidence", "html"))

# remove serious height/weight outliers and create flag columns
model_data <- model_data %>%
  # Rename originals
  rename(
    height_raw = height,
    weight_raw = weight
  ) %>%
  mutate(
    # -------------------------
    # HEIGHT CLEANING
    # -------------------------
    height_flag = case_when(
      is.na(height_raw) ~ NA,  # preserve original NA
      height_raw < 130 ~ TRUE,
      height_raw > 220 ~ TRUE,
      TRUE ~ FALSE
    ),
    height = case_when(
      height_flag == TRUE ~ NA_real_,
      TRUE ~ as.numeric(height_raw)
    ),
    # -------------------------
    # WEIGHT CLEANING
    # -------------------------
    weight_flag = case_when(
      is.na(weight_raw) ~ NA,
      weight_raw == 0 ~ TRUE,
      weight_raw < 35 ~ TRUE,
      weight_raw > 250 ~ TRUE,
      TRUE ~ FALSE
    ),
    weight = case_when(
      weight_flag == TRUE ~ NA_real_,
      TRUE ~ as.numeric(weight_raw)
    ),
    # -------------------------
    # BMI (only when valid)
    # -------------------------
    bmi = ifelse(!is.na(height) & !is.na(weight),
                 weight / (height / 100)^2,
                 NA_real_)
  )

# Save new model data -----------------------------------------------------
data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
saveRDS(model_data, file.path(data_dir,"PHI_modelData_v04.rds"))

model_data = readRDS(file.path(data_dir,"PHI_modelData_v04.rds"))

