# R/03_mortality_modeling.R
# Mortality prediction: From simple baselines to complex models
# Author: Jim Hitt
# Date: 2026-02-10


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

# Required packages for mortality modeling
required_packages <- c(
  # Core tidyverse
  "tidyverse",      # Data manipulation and visualization
  "dplyr",          # Data manipulation (included in tidyverse but explicit)
  "ggplot2",        # Visualization (included in tidyverse but explicit)
  
  # Modeling framework
  "tidymodels",     # Unified modeling interface (includes parsnip, recipes, tune, etc.)
  "parsnip",        # Model specifications (included in tidymodels)
  "recipes",        # Feature engineering (included in tidymodels)
  "rsample",        # Resampling/CV (included in tidymodels)
  "yardstick",      # Model metrics (included in tidymodels)
  "tune",           # Hyperparameter tuning (included in tidymodels)
  "workflows",      # Modeling workflows (included in tidymodels)
  
  # Specific model engines
  "xgboost",        # XGBoost algorithm
  "glmnet",         # Regularized regression (if we use penalized LR)
  
  # Model interpretation
  "vip",            # Variable importance plots
  "SHAPforxgboost", # SHAP values for XGBoost
  
  # Additional utilities
  "probably",       # Calibration plots
  "themis",         # Handling class imbalance if needed
  
  # Tables
  "gt",             # Publication-quality tables
  "gtsummary",       # Summary tables
  
  # Python
  "reticulate"
)

# Check, install if needed, and load
check_and_install(required_packages)

# clean up variables
rm(check_and_install, required_packages)


# Setup and Data Loading --------------------------------------------------

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(tidymodels)
library(duckdb)

# Setup ------------------------------------------------------------------------
set.seed(42)  # For reproducibility

model_data = readRDS('data/PHI_modelData_v01.rds')

#quick data check
glimpse(model_data)
model_data %>% count(inhosp_death_90day) %>% mutate(pct = n/sum(n)*100)


# Data preparation -------------------------------------------------------------
# A 0.44% event rate for 90-day mortality. 
# That's quite low, which will affect our modeling strategy (class imbalance considerations).

# Examine the outcome distribution
message("\n=== Outcome Distribution ===\n")
model_data %>% 
  count(inhosp_death_90day) %>% 
  mutate(
    pct = n/sum(n) * 100,
    pct = round(pct, 2)
  ) %>%
  print()

message(paste0(
  sprintf("\nTotal sample size: %d\n", nrow(model_data)),
  sprintf("Number of deaths: %d\n", sum(model_data$inhosp_death_90day)),
  sprintf("Event rate: %.2f%%\n", mean(model_data$inhosp_death_90day) * 100)
))

# Create stratified cross-validation folds
# Blog comment: Using 5-fold CV, stratified by outcome to preserve event rate in each fold
# Blog comment: Will discuss the importance of this tecnhique for low-incidence events in the blog

# change the outcome variable to a factor for the models
model_data = model_data %>%
  mutate(inhosp_death_90day = factor(inhosp_death_90day, levels = c(FALSE,TRUE)))

cv_folds <- vfold_cv(
  model_data, 
  v = 5, 
  strata = inhosp_death_90day
)

# Verify stratification worked
cat("\nThe stratified CV folds worked! (events per fold):\n")
cv_folds %>%
  mutate(
    n_events = map_int(splits, ~{
      analysis(.x) %>% 
        filter(inhosp_death_90day == TRUE) %>% 
        nrow()
    }),
    n_total = map_int(splits, ~nrow(analysis(.x)))
  ) %>%
  select(id, n_events, n_total) %>%
  mutate(event_rate = round(n_events/n_total * 100, 2)) %>%
  print()

# quick EDA of model variables
model_data %>% select(-subject_id, -op_id) %>% summary()

# Some suspicious values (height max 17,410? weight max 455?) - likely data quality issues we can filter

# Data cleaning ----------------------------------------------------------------
# Blog comment: We could have done this in the previous step, but it better to catch it now before modeling
cat("\n=== Data Cleaning ===\n")

# Set implausible values to NA (don't remove rows)
model_data_clean <- model_data %>%
  mutate(
    # Flag implausible values
    height = if_else(height < 100 | height > 220, NA_real_, height),
    weight = if_else(weight < 30 | weight > 200, NA_real_, weight),
    bmi = if_else(bmi < 10 | bmi > 60, NA_real_, bmi)
  )

message(paste0(
  sprintf("\nData quality summary:\n"),
  sprintf("Total observations: %d\n", nrow(model_data_clean)),
  sprintf("Total events: %d\n\n", sum(model_data_clean$inhosp_death_90day ==TRUE))
))

cat("Missing data by variable:\n")
model_data_clean %>%
  summarise(
    across(
      c(asa, age, sex, height, weight, bmi, emop,
        starts_with("preop_")),
      ~sum(is.na(.)),
      .names = "{.col}_missing"
    )
  ) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = round(n_missing / nrow(model_data_clean) * 100, 2),
    variable = str_remove(variable, "_missing")
  ) %>%
  arrange(desc(pct_missing)) %>%
  print(n = 20)

# recreate stratified cross-validation folds
# Blog comment: Using 5-fold CV, stratified by outcome to preserve event rate in each fold
# Blog comment: We addressed outlier values by converting them to Null/NaN
# Blog comment: Will discuss the importance of this tecnhique for low-incidence events in the blog

# Recreate CV folds with all data (NAs will be handled model-specifically)
cv_folds <- vfold_cv(
  model_data_clean, 
  v = 5, 
  strata = inhosp_death_90day
)

# Verify stratification worked
cat("\nThe stratified CV folds worked! (events per fold):\n")
cv_folds %>%
  mutate(
    n_events = map_int(splits, ~{
      analysis(.x) %>% 
        filter(inhosp_death_90day == TRUE) %>% 
        nrow()
    }),
    n_total = map_int(splits, ~nrow(analysis(.x)))
  ) %>%
  select(id, n_events, n_total) %>%
  mutate(event_rate = round(n_events/n_total * 100, 2)) %>%
  print()

# blog: what if we don't stratify for the low-incidence outcome?
cv_folds_unstratified <- vfold_cv(
  model_data_clean, 
  v = 5
)

cv_folds_unstratified %>%
  mutate(
    n_events = map_int(splits, ~{
      analysis(.x) %>% 
        filter(inhosp_death_90day == TRUE) %>% 
        nrow()
    }),
    n_total = map_int(splits, ~nrow(analysis(.x)))
  ) %>%
  select(id, n_events, n_total) %>%
  mutate(event_rate = round(n_events/n_total * 100, 2)) %>%
  print()

# The unstratfied CV folds are evenly distributed.
# ## When Does Stratification Matter?
#
# With 97,259 observations and 429 events, even unstratified random splitting 
# maintains consistent event rates (340-345 events per fold, range 0.44-0.44%).

# **Stratification becomes critical when:**
#   - Smaller sample sizes (n < 10,000)
#   - Rarer outcomes (< 0.1%)  
#   - Subgroup analyses with small cell counts

# **Our case:** Large enough that stratification provides minimal benefit for 
# event rate balance, but it's still **best practice** and costs nothing.

# Model 1: ASA Classification Only --------------------------------------------
cat("\n=== Model 1: ASA Classification Only ===\n")


# Define model specifications
lr_spec <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

xgb_spec <- boost_tree(trees = 100, tree_depth = 6) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Create recipe with just ASA (as factor) for LR
recipe_asa_lr <- recipe(inhosp_death_90day ~ asa, data = model_data_clean) %>%
  step_mutate(asa = factor(asa))

# Recipe for XGBoost (ASA as numeric)
recipe_asa_xgb <- recipe(inhosp_death_90day ~ asa, data = model_data_clean)

# Create workflows
wf_asa_lr <- workflow() %>%
  add_model(lr_spec) %>%
  add_recipe(recipe_asa_lr)

wf_asa_xgb <- workflow() %>%
  add_model(xgb_spec) %>%
  add_recipe(recipe_asa_xgb)

# Fit logistic regression
cat("\nFitting Logistic Regression (ASA only)...\n")
results_asa_lr <- fit_resamples(
  wf_asa_lr,
  resamples = cv_folds,
  metrics = metric_set(roc_auc, brier_class),
  control = control_resamples(save_pred = TRUE)
)

# Fit XGBoost
cat("Fitting XGBoost (ASA only)...\n")
results_asa_xgb <- fit_resamples(
  wf_asa_xgb,
  resamples = cv_folds,
  metrics = metric_set(roc_auc, brier_class),
  control = control_resamples(save_pred = TRUE)
)

# Collect and compare metrics
metrics_comparison <- bind_rows(
  collect_metrics(results_asa_lr) %>% mutate(model = "Logistic Regression"),
  collect_metrics(results_asa_xgb) %>% mutate(model = "XGBoost")
) %>%
  select(model, .metric, mean, std_err)

cat("\n=== Model Comparison ===\n")
print(metrics_comparison)

# Get AUROC values for both models
auroc_lr <- metrics_comparison %>% 
  filter(model == "Logistic Regression", .metric == "roc_auc") %>% 
  pull(mean)

auroc_xgb <- metrics_comparison %>% 
  filter(model == "XGBoost", .metric == "roc_auc") %>% 
  pull(mean)

message(paste0(
  sprintf("\nResults Summary\nLogistic Regression AUROC: %.3f\n", auroc_lr),
  sprintf("XGBoost AUROC: %.3f\n", auroc_xgb),
  sprintf("Difference: %.3f\n", auroc_xgb - auroc_lr)
))

# Generate figures for blog ------------------------------------------------
# Prepare output variables 
blog_name = "2026-02-11-ASAmortalityModel"
output_dir = file.path('posts', blog_name)
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# helper function to create numbered output paths
make_output_path = function(type, num, description, ext) {
  filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
  return(file.path(output_dir, filename))
}
rm(blog_name)
# --- end prepare output variables


# Fig01 ROC Curve Comparison ----------------------------------------------

preds_asa_lr <- collect_predictions(results_asa_lr)
preds_asa_xgb <- collect_predictions(results_asa_xgb)

# Fix ROC curve - use correct event level
roc_data <- bind_rows(
  roc_curve(preds_asa_lr, inhosp_death_90day, .pred_TRUE, event_level = "second") %>% 
    mutate(model = "Logistic Regression"),
  roc_curve(preds_asa_xgb, inhosp_death_90day, .pred_TRUE, event_level = "second") %>% 
    mutate(model = "XGBoost")
)

p_roc = ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, color = "gray50") +
  labs(
    title = "Model 1: ASA Classification Only",
    subtitle = "Both models achieve AUROC = 0.802",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(make_output_path('Fig',1,'ROCcurves','png'), 
       p_roc, width = 6, height = 5, dpi = 300)

# Blog Entry:
# With a single ordinal predictor, even a complex model like XGBoost can't outperform logistic regression. 
# There's simply no additional signal to extract. This changes dramatically when we add more features and interactions in later models."


# Table01 LR model coefficients -------------------------------------------

# Fit final LR model on full data to get coefficients
final_lr <- fit(wf_asa_lr, data = model_data_clean)

# Extract and format coefficients
coef_table <- tidy(final_lr) %>%
  mutate(
    OR = exp(estimate),
    OR_lower = exp(estimate - 1.96 * std.error),
    OR_upper = exp(estimate + 1.96 * std.error),
    p_value = format.pval(p.value, digits = 3)
  ) %>%
  select(term, estimate, std.error, OR, OR_lower, OR_upper, p_value)

# Create HTML table
library(gt)

gt_table <- coef_table %>%
  gt() %>%
  tab_header(
    title = "Logistic Regression Coefficients",
    subtitle = "ASA Classification predicting 90-day mortality"
  ) %>%
  fmt_number(
    columns = c(estimate, std.error, OR, OR_lower, OR_upper),
    decimals = 3
  ) %>%
  cols_label(
    term = "Term",
    estimate = "Coefficient",
    std.error = "Std Error",
    OR = "Odds Ratio",
    OR_lower = "OR 95% Lower",
    OR_upper = "OR 95% Upper",
    p_value = "P-value"
  )

# Save HTML table
gtsave(gt_table, 
       make_output_path('Table',1,'LRcoefficients','html'))

cat("\nCoefficient table saved to output/SAFE_asa_coefficients.html\n")

# Save clean model data for Python --------------------------------------------
write_csv(model_data_clean, "data/PHI_model_data_v01.csv")
cat("\nSaved model_data_clean.csv for Python\n")

# save the clean model data (with outliers addressed)
write_rds(model_data_clean, 'PHI_modelData_v02.rds')


# Execute Python Code from R ----------------------------------------------
# Blog: I will start running Python code from RStudio. I find RStudio is excellent for
#   data handling, and Python is great for modeling and data analysis, so I will combine
#   them and get the best of both worlds.

library(reticulate)
use_virtualenv("Python/.venv")
system("python Python/01_mortality_asa_xgb.py data/PHI_model_data_v01.csv")
