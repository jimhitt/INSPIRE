# INSPIRE Dataset - Multiple Operations Analysis
# Purpose: Characterize patients with repeat surgeries and time intervals between operations
# Output: Tables and figures for blog post
# Date: 2026-02-06

# Install packages in renv -----------------------------------------------
# renv = "R environment" - creates project-specific package libraries
# Like Python's venv/conda, it isolates packages per project for reproducibility
# 
# Why use renv?
# - Ensures exact package versions are recorded in renv.lock
# - Anyone can restore your exact environment with renv::restore()
# - Prevents package conflicts between different R projects
# - Essential for reproducible research and collaboration
#
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
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

data_dir = "data"
ops_file = "PHI_operations.csv"
ops_path = file.path(data_dir, ops_file)

# Verify file exists before trying to load
if (!file.exists(ops_path)) {
  stop(sprintf("Data file not found: %s\nPlease ensure PHI data is in the data/ directory", ops_path))
}

# Load operations table
ops <- dbGetQuery(con, sprintf("SELECT * FROM read_csv_auto('%s')", ops_path))

# Close connection (we'll work in memory from here)
dbDisconnect(con, shutdown = TRUE)
rm(con, data_dir, ops_file)


# Prepare output variables ------------------------------------------------

analysis_num = "01"
analysis_name = "multiple_operations"
output_dir = file.path("output", sprintf("%s_%s",analysis_num, analysis_name))
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# helper function to create numbered output paths
make_output_path = function(type, num, description, ext) {
  filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
  return(file.path(output_dir, filename))
}
# Calculate operation sequences ------------------------------------------

# For each subject, number their operations chronologically
ops_sequenced <- ops %>%
  arrange(subject_id, anstart_time) %>%
  group_by(subject_id) %>%
  mutate(
    operation_number = row_number(),
    total_operations = n(),
    is_first_operation = (operation_number == 1),
    multiple_operations = (total_operations > 1)
  ) %>%
  ungroup()

# Calculate time between consecutive operations (in days)
ops_intervals <- ops_sequenced %>%
  filter(total_operations > 1) %>%
  arrange(subject_id, operation_number) %>%
  group_by(subject_id) %>%
  mutate(
    prev_anend_time = lag(anend_time),
    days_since_prev = (anstart_time - prev_anend_time) / (60 * 24)  # Convert minutes to days
  ) %>%
  ungroup() %>%
  filter(!is.na(days_since_prev))  # Remove first operation (no previous)


# Table 1: Distribution of operations per subject -----------------------

# Subject-level summary
subject_summary <- ops_sequenced %>%
  group_by(subject_id) %>%
  summarise(
    total_operations = first(total_operations),
    .groups = "drop"
  )

# Count by number of operations
ops_distribution <- subject_summary %>%
  mutate(
    ops_category = case_when(
      total_operations == 1 ~ "1 operation",
      total_operations == 2 ~ "2 operations",
      total_operations == 3 ~ "3 operations",
      total_operations == 4 ~ "4 operations",
      total_operations >= 5 ~ "5+ operations",
      TRUE ~ "Unknown"
    ),
    ops_category = factor(ops_category, levels = c("1 operation", "2 operations", 
                                                   "3 operations", "4 operations", 
                                                   "5+ operations"))
  ) %>%
  count(ops_category, name = "n_subjects") %>%
  mutate(
    percentage = n_subjects / sum(n_subjects) * 100,
    pct_label = sprintf("%.1f%%", percentage)
  )

# Create publication table
table1_operations_dist <- ops_distribution %>%
  gt() %>%
  cols_label(
    ops_category = "Operations per Subject",
    n_subjects = "Count",
    percentage = "Percentage",
    pct_label = "Formatted"
  ) %>%
  fmt_number(
    columns = n_subjects,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  fmt_number(
    columns = percentage,
    decimals = 1
  ) %>%
  cols_hide(columns = percentage) %>%
  cols_move_to_end(columns = pct_label) %>%
  cols_label(pct_label = "Percentage") %>%
  tab_header(
    title = "Distribution of Operations per Subject",
    subtitle = "INSPIRE Dataset (N=99,886 subjects)"
  ) %>%
  tab_source_note(
    source_note = "78.4% of subjects had a single operation during the study period"
  )

# Save table
gtsave(table1_operations_dist, 
       make_output_path("table",1,"operations_distribution", "html"),
       inline_css = FALSE
       )

# Print summary statistics
message("\n=== OPERATIONS PER SUBJECT ===\n")
print(ops_distribution)

cat("\nTotal subjects with 2+ operations:", 
    sum(ops_distribution$n_subjects[ops_distribution$ops_category != "1 operation"]),
    sprintf("(%.1f%%)", 
            sum(ops_distribution$percentage[ops_distribution$ops_category != "1 operation"])))


# Table 2: Inter-operation time intervals --------------------------------

# Summary statistics for time between operations
interval_summary <- ops_intervals %>%
  summarise(
    n_intervals = n(),
    mean_days = mean(days_since_prev, na.rm = TRUE),
    sd_days = sd(days_since_prev, na.rm = TRUE),
    median_days = median(days_since_prev, na.rm = TRUE),
    p01 = quantile(days_since_prev, 0.01, na.rm = TRUE),
    p05 = quantile(days_since_prev, 0.05, na.rm = TRUE),
    p10 = quantile(days_since_prev, 0.10, na.rm = TRUE),
    p25 = quantile(days_since_prev, 0.25, na.rm = TRUE),
    p75 = quantile(days_since_prev, 0.75, na.rm = TRUE),
    p90 = quantile(days_since_prev, 0.90, na.rm = TRUE),
    p95 = quantile(days_since_prev, 0.95, na.rm = TRUE),
    p99 = quantile(days_since_prev, 0.99, na.rm = TRUE),
    max_days = max(days_since_prev, na.rm = TRUE)
  )

# Create table of interval statistics
# -- DO NOT FIND THIS TABLE HELPFUL -----
table2_intervals <- interval_summary %>%
  pivot_longer(
    cols = everything(),
    names_to = "statistic",
    values_to = "days"
  ) %>%
  mutate(
    statistic_label = case_when(
      statistic == "n_intervals" ~ "Number of interval pairs",
      statistic == "mean_days" ~ "Mean",
      statistic == "sd_days" ~ "Standard deviation",
      statistic == "median_days" ~ "Median",
      statistic == "p01" ~ "1st percentile",
      statistic == "p05" ~ "5th percentile",
      statistic == "p10" ~ "10th percentile",
      statistic == "p25" ~ "25th percentile",
      statistic == "p75" ~ "75th percentile",
      statistic == "p90" ~ "90th percentile",
      statistic == "p95" ~ "95th percentile",
      statistic == "p99" ~ "99th percentile",
      statistic == "max_days" ~ "Maximum",
      TRUE ~ statistic
    )
  ) %>%
  gt() %>%
  cols_label(
    statistic_label = "Statistic",
    days = "Days"
  ) %>%
  fmt_number(
    columns = days,
    decimals = 1,
    use_seps = TRUE
  ) %>%
  tab_header(
    title = "Time Intervals Between Consecutive Operations",
    subtitle = sprintf("Based on %s operation pairs", 
                       format(nrow(ops_intervals), big.mark = ","))
  ) %>%
  tab_source_note(
    source_note = "Median interval: 208 days (IQR: 28-637 days)"
  )

gtsave(table2_intervals, 
       filename = make_output_path("table",2,"interval_stats","html"),
       inline_css = FALSE)

message("\n=== INTER-OPERATION INTERVALS ===\n")
print(interval_summary)


# Table 3: Categorical time intervals -----------------------------------

# Categorize intervals
interval_categories <- ops_intervals %>%
  mutate(
    interval_category = case_when(
      days_since_prev < 1 ~ "Same day (< 1 day)",
      days_since_prev < 7 ~ "1-7 days",
      days_since_prev < 30 ~ "7-30 days",
      days_since_prev < 90 ~ "30-90 days",
      days_since_prev < 180 ~ "90-180 days",
      days_since_prev < 365 ~ "180-365 days",
      days_since_prev >= 365 ~ "> 365 days",
      TRUE ~ "Unknown"
    ),
    interval_category = factor(interval_category, 
                               levels = c("Same day (< 1 day)", "1-7 days", "7-30 days",
                                          "30-90 days", "90-180 days", "180-365 days",
                                          "> 365 days"))
  ) %>%
  count(interval_category, name = "n_intervals") %>%
  mutate(
    percentage = n_intervals / sum(n_intervals) * 100,
    pct_label = sprintf("%.1f%%", percentage)
  )

table3_interval_categories <- interval_categories %>%
  gt() %>%
  cols_label(
    interval_category = "Time Between Operations",
    n_intervals = "Count",
    pct_label = "Percentage"
  ) %>%
  fmt_number(
    columns = n_intervals,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  cols_hide(columns = percentage) %>%
  tab_header(
    title = "Distribution of Inter-Operation Intervals",
    subtitle = "Categorized by clinically relevant time periods"
  ) %>%
  tab_source_note(
    source_note = "5.8% of operations occurred within 7 days of previous surgery"
  )

gtsave(table3_interval_categories, 
       filename = make_output_path("table",3,"interval_categories","html"),
       inline_css = FALSE)

message("\n=== INTERVAL CATEGORIES ===\n")
print(interval_categories)


# Figure 1: Operations per subject (bar chart) --------------------------

fig1_ops_distribution <- ggplot(ops_distribution, 
                                aes(x = ops_category, y = n_subjects, fill = ops_category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(n_subjects)), 
            vjust = -0.5, size = 3.5) +
  geom_text(aes(label = pct_label), 
            vjust = 1.5, size = 3, color = "white", fontface = "bold") +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(
    title = "Distribution of Operations per Subject",
    subtitle = "INSPIRE Dataset (N=99,886 subjects)",
    x = "Number of Operations",
    y = "Number of Subjects",
    caption = "78.4% of subjects had a single operation during 2011-2020"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(make_output_path("fig",1,"operations_distribution","png"),
       plot = fig1_ops_distribution,
       width = 8, height = 6, dpi = 300, bg = "white")

print(fig1_ops_distribution)


# Figure 2: Interval distribution (histogram with density) --------------

fig2_interval_histogram <- ggplot(ops_intervals, aes(x = days_since_prev)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = median(ops_intervals$days_since_prev), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", 
           x = median(ops_intervals$days_since_prev) + 100, 
           y = max(density(ops_intervals$days_since_prev)$y) * 0.9,
           label = sprintf("Median: %d days", round(median(ops_intervals$days_since_prev))),
           color = "red", fontface = "bold", hjust = 0) +
  scale_x_continuous(
    limits = c(0, 1500),
    breaks = seq(0, 1500, 250),
    labels = comma
  ) +
  labs(
    title = "Time Between Consecutive Operations",
    subtitle = "Distribution shows right skew with long tail",
    x = "Days Since Previous Operation",
    y = "Density",
    caption = sprintf("Based on %s operation pairs | Median: %d days | Mean: %d days",
                      format(nrow(ops_intervals), big.mark = ","),
                      round(median(ops_intervals$days_since_prev)),
                      round(mean(ops_intervals$days_since_prev)))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(make_output_path("fig",2,"interval_histogram","png"),
       plot = fig2_interval_histogram,
       width = 10, height = 6, dpi = 300, bg = "white")

print(fig2_interval_histogram)


# Figure 3: Interval categories (bar chart with annotation) -------------

fig3_interval_categories <- ggplot(interval_categories, 
                                   aes(x = interval_category, y = n_intervals, fill = interval_category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(n_intervals)), 
            vjust = -0.5, size = 3) +
  geom_text(aes(label = pct_label), 
            vjust = 1.5, size = 3, color = "white", fontface = "bold") +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
  labs(
    title = "Time Between Operations: Categorical Distribution",
    subtitle = "Clinically relevant time periods",
    x = "Time Category",
    y = "Number of Operation Pairs",
    caption = "Nearly 6% of operations occurred within 7 days of previous surgery"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(make_output_path("fig",3,"interval_categories","png"),
       plot = fig3_interval_categories,
       width = 10, height = 6, dpi = 300, bg = "white")

print(fig3_interval_categories)


# Flag first operations and reoperation outcome ---------------------------
ops_analysis <- ops_sequenced %>%
  mutate(
    # Primary analysis cohort flag
    is_first_operation = (operation_number == 1),
    
    # Patient-level: Did this subject have multiple operations?
    multiple_operations = (total_operations > 1)
  )

# For first operations: Flag if reoperation occurred within 30 days ------
# This is a POSTOPERATIVE OUTCOME, not a preoperative predictor

first_operations <- ops_analysis %>%
  filter(is_first_operation == TRUE)

# Identify subjects who had a second operation within 30 days of first
reoperation_30d_subjects <- ops_intervals %>%
  filter(operation_number == 2, days_since_prev < 30) %>%
  pull(subject_id)

# Add reoperation flag to first operations dataset
first_operations <- first_operations %>%
  mutate(
    reoperation_within_30d = subject_id %in% reoperation_30d_subjects
  )

# Table 4: Summary of Multiple Operations Analysis -----------------------

summary_data <- tibble(
  Metric = c(
    "Total operations in dataset",
    "Total subjects",
    "Subjects with single operation",
    "Subjects with multiple operations",
    "Operations from multi-surgery subjects",
    "First operations (primary analysis cohort)",
    "Subjects with reoperation within 30 days"
  ),
  Count = c(
    nrow(ops_sequenced),
    length(unique(ops_sequenced$subject_id)),
    sum(subject_summary$total_operations == 1),
    sum(subject_summary$total_operations > 1),
    sum(!ops_sequenced$is_first_operation),
    nrow(first_operations),
    sum(first_operations$reoperation_within_30d)
  ),
  Percentage = c(
    NA,  # Total operations (no percentage)
    NA,  # Total subjects (no percentage)
    sum(subject_summary$total_operations == 1) / nrow(subject_summary) * 100,
    sum(subject_summary$total_operations > 1) / nrow(subject_summary) * 100,
    sum(!ops_sequenced$is_first_operation) / nrow(ops_sequenced) * 100,
    nrow(first_operations) / nrow(ops_sequenced) * 100,
    sum(first_operations$reoperation_within_30d) / nrow(first_operations) * 100
  )
)

table4_summary <- summary_data %>%
  gt() %>%
  cols_label(
    Metric = "Metric",
    Count = "N",
    Percentage = "Percentage"
  ) %>%
  fmt_number(
    columns = c(Count),
    decimals = 0,
    use_seps = TRUE
  ) %>%
  fmt_number(
    columns = c(Percentage),
    decimals = 1,
    pattern = "{x}%"
  ) %>%
  sub_missing(
    columns = c(Percentage),
    missing_text = "—"
  ) %>%
  tab_header(
    title = "Multiple Operations Summary",
    subtitle = "INSPIRE Dataset (2011-2020)"
  ) %>%
  tab_source_note(
    source_note = sprintf("Primary analysis cohort: first operations only (N=%s, %.1f%% of all operations)", 
                          format(nrow(first_operations), big.mark = ","),
                          nrow(first_operations) / nrow(ops_sequenced) * 100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(Metric),
      rows = Metric == "First operations (primary analysis cohort)"
    )
  )

# Save table
gtsave(table4_summary, 
       make_output_path("table", 4, "summary_statistics", "html"),
       inline_css = FALSE)

message("\n=== SUMMARY TABLE ===")
print(summary_data)

# Save analysis-ready dataset --------------------------------------------

# Primary analysis cohort: first operations only
saveRDS(first_operations, "data/PHI_first_operations.rds")

# Complete operations dataset with analysis flags
saveRDS(ops_analysis, "data/PHI_operations_analysis.rds")
