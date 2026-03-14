# R/04_multiVar_mortality_modeling.R
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
library(gt)

# Setup ------------------------------------------------------------------------
set.seed(42)  # For reproducibility

data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
model_data = readRDS(file.path(data_dir,'PHI_modelData_v01.rds'))


#quick data check
glimpse(model_data)
model_data %>% count(inhosp_death_90day) %>% mutate(pct = n/sum(n)*100)


# Exploring the predictor variables ---------------------------------------
# there are some errors (low and high outliers) in height and weight
# I will truncate those variables to a reasonable range
# There mey be signal in the errors, so I will create a height_ and weight_flag column
# for the LR model (XGBoost can take Null values)
# the lab data range appears physiologic without significant errant outliers.

input_vars = c('age', 'sex', 'height','weight','bmi', 'emop', 'andur', 'preop_hb', 'preop_platelet', 
               'preop_wbc', 'preop_aptt', 'preop_ptinr', 'preop_glucose', 'preop_bun', 'preop_albumin',
               'preop_ast', 'preop_alt', 'preop_creatinine', 'preop_sodium', 'preop_potassium' )

summary_tbl <- model_data %>%
  summarise(across(all_of(input_vars),
                   list(
                     n   = ~sum(!is.na(.)),
                     mean= ~mean(., na.rm = TRUE),
                     sd  = ~sd(., na.rm = TRUE),
                     min = ~min(., na.rm = TRUE),
                     p1  = ~quantile(., 0.01, na.rm = TRUE),
                     p10  = ~quantile(., 0.10, na.rm = TRUE),
                     p25  = ~quantile(., 0.25, na.rm = TRUE),
                     p50  = ~quantile(., 0.50, na.rm = TRUE),
                     p75  = ~quantile(., 0.75, na.rm = TRUE),
                     p90  = ~quantile(., 0.90, na.rm = TRUE),
                     p99 = ~quantile(., 0.99, na.rm = TRUE),
                     max = ~max(., na.rm = TRUE)
                   ),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("variable", "stat"),
               names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(variable)

summary_tbl

# create table for blog
if (1){
  # Prepare output variables 
  blog_name = "2026-02-12-multiVarMortalityyModel"
  output_dir = file.path('posts', blog_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = F)
  
  # helper function to create numbered output paths
  make_output_path = function(type, num, description, ext) {
    filename = sprintf("%s%02d_%s.%s", type, num, description, ext)
    return(file.path(output_dir, filename))
  }
  rm(blog_name)
  # --- end prepare output variables
  
  binary_vars <- c("sex", "emop")
  
  # OPTIONAL: variable display labels (edit as you like)
  var_labels <- c(
    age = "Age (years)",
    sex = "Sex (1=Male)",
    height = "Height (cm)",
    weight = "Weight (kg)",
    bmi = "BMI (kg/m²)",
    emop = "Emergency operation (1=Yes)",
    andur = "Anesthesia duration (min)",
    preop_hb = "Preop Hemoglobin (g/dL)",
    preop_platelet = "Preop Platelets (10^3/µL)",
    preop_wbc = "Preop WBC (10^3/µL)",
    preop_aptt = "Preop aPTT (s)",
    preop_ptinr = "Preop PT-INR",
    preop_glucose = "Preop Glucose (mg/dL)",
    preop_bun = "Preop BUN (mg/dL)",
    preop_albumin = "Preop Albumin (g/dL)",
    preop_ast = "Preop AST (U/L)",
    preop_alt = "Preop ALT (U/L)",
    preop_creatinine = "Preop Creatinine (mg/dL)",
    preop_sodium = "Preop Sodium (mmol/L)",
    preop_potassium = "Preop Potassium (mmol/L)"
  )
  # ----- Table 1: Binary Variables -----
  bin_tbl <- model_data %>%
    summarise(
      across(all_of(binary_vars), list(
        n = ~sum(!is.na(.)),
        n1 = ~sum(. == 1, na.rm = TRUE),
        p1 = ~mean(. == 1, na.rm = TRUE)
      ), .names = "{.col}__{.fn}")
    ) %>%
    pivot_longer(
      everything(),
      names_to = c("variable", "stat"),
      names_sep = "__",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(
      Variable = dplyr::coalesce(unname(var_labels[variable]), variable),
      N = as.integer(round(n)),
      `N (Yes)` = as.integer(round(n1)),
      `% (Yes)` = 100 * p1
    ) %>%
    select(Variable, N, `N (Yes)`, `% (Yes)`)
  
  binary_gt <- bin_tbl %>%
    gt() %>%
    tab_header(
      title = "Binary Predictors",
      subtitle = "Counts and percentages for categorical variables"
    ) %>%
    fmt_integer(columns = c(N, `N (Yes)`)) %>%
    fmt_number(columns = `% (Yes)`, decimals = 1) %>%
    opt_row_striping() %>%
    opt_table_outline() %>%
    tab_options(
      table.width = pct(100),
      data_row.padding = px(6),
      column_labels.font.weight = "600",
      table.font.size = px(14)
    )
  
  gtsave(binary_gt, make_output_path("Table", 1, "BinaryPredictors", "html"))
  
  # ----- Table 2: Continuous Variables -----
  cont_tbl <- summary_tbl %>%
    filter(!variable %in% binary_vars) %>%
    mutate(
      Variable = dplyr::coalesce(unname(var_labels[variable]), variable),
      N = as.integer(round(n))
    ) %>%
    transmute(
      Variable, N,
      Mean = mean, SD = sd,
      Min = min, P25 = p25, Median = p50, P75 = p75, Max = max
    )
  
  continuous_gt <- cont_tbl %>%
    gt() %>%
    tab_header(
      title = "Continuous Predictors",
      subtitle = "Distribution summaries for numeric variables"
    ) %>%
    fmt_integer(columns = N) %>%
    fmt_number(columns = c(Mean, SD, Min, P25, Median, P75, Max),
               decimals = 2, use_seps = TRUE) %>%
    opt_row_striping() %>%
    opt_table_outline() %>%
    tab_options(
      table.width = pct(100),
      data_row.padding = px(6),
      column_labels.font.weight = "600",
      table.font.size = px(14)
    )
  
  gtsave(continuous_gt, make_output_path("Table", 2, "ContinuousPredictors", "html")) 
 
}

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

# save model data with corrected height, weight, BMI
data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
saveRDS(model_data,file.path(data_dir, "PHI_modelData_v02.rds"))

# Modeling ----------------------------------------------------------------
library(tidyverse)
library(tidymodels)
tidymodels_prefer()
set.seed(2026)

# typing function
type_model_data_clean <- function(df) {
  df %>%
    mutate(
      # binary predictors: treat as categorical, not continuous
      sex  = if ("sex"  %in% names(.)) factor(sex,  levels = c(0, 1)) else sex,
      emop = if ("emop" %in% names(.)) factor(emop, levels = c(0, 1)) else emop,
      
      # ordinal predictor
      asa  = if ("asa"  %in% names(.)) factor(asa, levels = 1:5, ordered = TRUE) else asa,
      
      # optional: force logical flags to 0/1 numeric (recommended if you model them)
      height_flag = if ("height_flag" %in% names(.)) as.integer(height_flag) else height_flag,
      weight_flag = if ("weight_flag" %in% names(.)) as.integer(weight_flag) else weight_flag
    )
}

model_data_clean = type_model_data_clean(model_data)

# ---------------------------
# Shared helpers for modeling
# ---------------------------
set.seed(2026)

mort_metrics <- metric_set(roc_auc, pr_auc, brier_class)
ctrl <- control_resamples(save_pred = TRUE)

prep_outcome <- function(df, outcome = "inhosp_death_90day") {
  df %>%
    mutate(
      !!outcome := factor(if_else(.data[[outcome]], "yes", "no"),
                          levels = c("no", "yes"))
    )
}

build_formula <- function(data, outcome, input_vars) {
  input_vars <- intersect(input_vars, names(data))
  as.formula(paste(outcome, "~", paste(input_vars, collapse = " + ")))
}

make_recipe_lr <- function(formula, data) {
  recipe(formula, data = data) %>%
    step_indicate_na(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_zv(all_predictors())
}

make_recipe_xgb <- function(formula, data) {
  recipe(formula, data = data) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_zv(all_predictors())
}

run_lr_resamples <- function(data, input_vars, cv_folds,
                             outcome = "inhosp_death_90day",
                             engine = c("glmnet", "glm"),
                             penalty = 0.01, mixture = 1,
                             metrics = mort_metrics, control = ctrl) {
  
  engine <- match.arg(engine)
  
  f <- build_formula(data, outcome, input_vars)
  rec <- make_recipe_lr(f, data)
  
  spec <- if (engine == "glmnet") {
    logistic_reg(penalty = penalty, mixture = mixture) %>%
      set_engine("glmnet") %>%
      set_mode("classification")
  } else {
    logistic_reg() %>%
      set_engine("glm") %>%
      set_mode("classification")
  }
  
  wf <- workflow() %>% add_model(spec) %>% add_recipe(rec)
  res <- fit_resamples(wf, resamples = cv_folds, metrics = metrics, control = control)
  
  list(
    model = paste0("LR (", engine, ")"),
    workflow = wf,
    resamples = res,
    metrics = collect_metrics(res),
    preds = collect_predictions(res)
  )
}

run_xgb_resamples <- function(data, input_vars, cv_folds,
                              outcome = "inhosp_death_90day",
                              trees = 800, tree_depth = 6, learn_rate = 0.05,
                              min_n = 10, loss_reduction = 0,
                              sample_size = 0.8, mtry_prop = 0.8,
                              metrics = mort_metrics, control = ctrl) {
  
  # class imbalance for scale_pos_weight
  pos <- sum(data[[outcome]] == "yes", na.rm = TRUE)
  neg <- sum(data[[outcome]] == "no",  na.rm = TRUE)
  spw <- if (pos > 0) neg / pos else 1
  
  f <- build_formula(data, outcome, input_vars)
  used_vars <- all.vars(f)[-1]
  mtry_val <- max(1, floor(length(used_vars) * mtry_prop))
  
  rec <- make_recipe_xgb(f, data)
  
  spec <- boost_tree(
    trees = trees,
    tree_depth = tree_depth,
    learn_rate = learn_rate,
    min_n = min_n,
    loss_reduction = loss_reduction,
    sample_size = sample_size,
    mtry = mtry_val
  ) %>%
    set_engine("xgboost", eval_metric = "logloss", scale_pos_weight = spw) %>%
    set_mode("classification")
  
  wf <- workflow() %>% add_model(spec) %>% add_recipe(rec)
  res <- fit_resamples(wf, resamples = cv_folds, metrics = metrics, control = control)
  
  list(
    model = "XGB",
    settings = list(scale_pos_weight = spw, mtry = mtry_val),
    workflow = wf,
    resamples = res,
    metrics = collect_metrics(res),
    preds = collect_predictions(res)
  )
}

compare_model_metrics <- function(...) {
  runs <- list(...)
  bind_rows(lapply(runs, function(x) {
    x$metrics %>% mutate(model = x$model) %>% select(model, .metric, mean, std_err)
  })) %>%
    arrange(.metric, desc(mean))
}

# -------------------------
# One-time: outcome + typing
# -------------------------
model_data_clean2 <- model_data %>%
  prep_outcome("inhosp_death_90day") %>%
  type_model_data_clean()  # <- from earlier snippet

# 2) recreate folds on the finalized data
set.seed(2026)
cv_folds2 <- vfold_cv(model_data_clean2, v = 5, strata = inhosp_death_90day)


# small feature set -------------------------------------------------------

input_vars <- c(
  "age", "sex", "bmi"
)

# fit with the new folds
lr_run  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                            engine = "glm", penalty = 0.01, mixture = 1)

xgb_run <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

metrics_comparison01 <- compare_model_metrics(lr_run, xgb_run)
print(metrics_comparison01)

# Collect predictions with labels:
lr_run01  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                              engine = "glm", penalty = 0.01, mixture = 1)
xgb_run01 <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

# Collect predictions with metadata
preds_01 <- bind_rows(
  lr_run01$preds %>% mutate(model = "LR", feature_set = "Demographics"),
  xgb_run01$preds %>% mutate(model = "XGB", feature_set = "Demographics")
)

# medium feature set ------------------------------------------------------

input_vars <- c(
  "age", "sex", "bmi", "asa", 
  "emop", "andur"
)

# fit with the new folds
lr_run  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                            engine = "glm", penalty = 0.01, mixture = 1)

xgb_run <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

metrics_comparison02 <- compare_model_metrics(lr_run, xgb_run)
print(metrics_comparison02)

# Collect predictions with labels:
lr_run02  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                              engine = "glm", penalty = 0.01, mixture = 1)
xgb_run02 <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

# Collect predictions with metadata
preds_02 <- bind_rows(
  lr_run02$preds %>% mutate(model = "LR", feature_set = "ClinicalSeverity"),
  xgb_run02$preds %>% mutate(model = "XGB", feature_set = "ClinicalSeverity")
)

# labs and age ------------------------------------------------------------

input_vars <- c(
  "age", 
  "preop_hb", "preop_platelet", "preop_wbc",
  "preop_aptt", "preop_ptinr", "preop_glucose",
  "preop_bun", "preop_albumin", "preop_ast", "preop_alt",
  "preop_creatinine", "preop_sodium", "preop_potassium"
)

# fit with the new folds
lr_run  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                            engine = "glm", penalty = 0.01, mixture = 1)

xgb_run <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

metrics_comparison03 <- compare_model_metrics(lr_run, xgb_run)
print(metrics_comparison03)

# Collect predictions with labels:
lr_run03  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                              engine = "glm", penalty = 0.01, mixture = 1)
xgb_run03 <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

# Collect predictions with metadata
preds_03 <- bind_rows(
  lr_run03$preds %>% mutate(model = "LR", feature_set = "Labs"),
  xgb_run03$preds %>% mutate(model = "XGB", feature_set = "Labs")
)

# full feature set --------------------------------------------------------

input_vars <- c(
  "age", "sex", "bmi", "asa", 
  "emop", "andur",
  "height_flag", "weight_flag",
  "preop_hb", "preop_platelet", "preop_wbc",
  "preop_aptt", "preop_ptinr", "preop_glucose",
  "preop_bun", "preop_albumin", "preop_ast", "preop_alt",
  "preop_creatinine", "preop_sodium", "preop_potassium"
)

# fit with the new folds
lr_run  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                            engine = "glm", penalty = 0.01, mixture = 1)

xgb_run <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

metrics_comparison04 <- compare_model_metrics(lr_run, xgb_run)
print(metrics_comparison04)

# Collect predictions with labels:
lr_run04  <- run_lr_resamples(model_data_clean2, input_vars, cv_folds2,
                              engine = "glm", penalty = 0.01, mixture = 1)
xgb_run04 <- run_xgb_resamples(model_data_clean2, input_vars, cv_folds2)

# Collect predictions with metadata
preds_04 <- bind_rows(
  lr_run04$preds %>% mutate(model = "LR", feature_set = "FullFeatures"),
  xgb_run04$preds %>% mutate(model = "XGB", feature_set = "FullFeatures")
)

print(metrics_comparison01)
print(metrics_comparison02)
print(metrics_comparison03)
print(metrics_comparison04)

# Combine all metrics with feature set labels ===================
metrics_combined <- bind_rows(
  metrics_comparison01 %>% mutate(feature_set = "Demographics"),
  metrics_comparison02 %>% mutate(feature_set = "Clinical Severity"),
  metrics_comparison03 %>% mutate(feature_set = "Laboratories"),
  metrics_comparison04 %>% mutate(feature_set = "Full Model")
) %>%
  mutate(
    feature_set = factor(feature_set, levels = c("Demographics", "Clinical Severity", 
                                                 "Laboratories", "Full Model"))
  )

# Reshape for better table display
metrics_wide <- metrics_combined %>%
  select(feature_set, model, .metric, mean, std_err) %>%
  pivot_wider(
    names_from = .metric,
    values_from = c(mean, std_err),
    names_glue = "{.metric}_{.value}"
  ) %>%
  select(feature_set, model, 
         roc_auc_mean, roc_auc_std_err,
         pr_auc_mean, pr_auc_std_err,
         brier_class_mean, brier_class_std_err) %>%
  arrange(feature_set, desc(roc_auc_mean))

# Create gt table
gt_metrics <- metrics_wide %>%
  gt(groupname_col = "feature_set") %>%
  tab_header(
    title = "Model Performance Across Feature Sets",
    subtitle = "5-fold cross-validation results for mortality prediction"
  ) %>%
  cols_label(
    model = "Model",
    roc_auc_mean = "AUC",
    roc_auc_std_err = "SE",
    pr_auc_mean = "PR-AUC",
    pr_auc_std_err = "SE",
    brier_class_mean = "Brier",
    brier_class_std_err = "SE"
  ) %>%
  tab_spanner(
    label = "ROC AUC",
    columns = c(roc_auc_mean, roc_auc_std_err)
  ) %>%
  tab_spanner(
    label = "Precision-Recall AUC",
    columns = c(pr_auc_mean, pr_auc_std_err)
  ) %>%
  tab_spanner(
    label = "Brier Score",
    columns = c(brier_class_mean, brier_class_std_err)
  ) %>%
  fmt_number(
    columns = c(roc_auc_mean, pr_auc_mean),
    decimals = 3
  ) %>%
  fmt_number(
    columns = c(roc_auc_std_err, pr_auc_std_err),
    decimals = 4
  ) %>%
  fmt_number(
    columns = c(brier_class_mean, brier_class_std_err),
    decimals = 5
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  opt_row_striping() %>%
  opt_table_outline() %>%
  tab_options(
    table.width = pct(100),
    data_row.padding = px(6),
    column_labels.font.weight = "600",
    table.font.size = px(14)
  ) %>%
  tab_footnote(
    footnote = "Lower Brier scores indicate better calibration (closer to 0 is better)",
    locations = cells_column_spanners(spanners = "Brier Score")
  )

# Save
gtsave(gt_metrics, make_output_path("Table", 3, "ModelPerformance", "html"))

metrics_all <- bind_rows(
  metrics_comparison01 %>% mutate(feature_set = "Demographics"),
  metrics_comparison02 %>% mutate(feature_set = "Clinical Severity"),
  metrics_comparison03 %>% mutate(feature_set = "Laboratories"),
  metrics_comparison04 %>% mutate(feature_set = "Full Model")
)


metrics_all = metrics_all %>%
  mutate(feature_set = factor(
    feature_set,
    levels = c(
      "Demographics",
      "Clinical Severity",
      "Laboratories",
      "Full Model"
    )))

metrics_all %>%
  filter(.metric == "roc_auc") %>%
  ggplot(aes(x = feature_set, y = mean, color = model, group = model)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mortality Discrimination Across Feature Groups",
    x = "Feature Group",
    y = "ROC AUC"
  )


roc_progression <- metrics_all %>%
  filter(.metric == "roc_auc") %>%
  select(feature_set, model, mean) %>%
  arrange(model, feature_set) %>%
  group_by(model) %>%
  mutate(delta_auc = mean - dplyr::lag(mean))

roc_progression
# > roc_progression
# # A tibble: 8 × 4
# # Groups:   model [2]
# feature_set       model     mean delta_auc
# <fct>             <chr>    <dbl>     <dbl>
#   1 Demographics      LR (glm) 0.763  NA      
# 2 Clinical Severity LR (glm) 0.901   0.138  
# 3 Laboratories      LR (glm) 0.905   0.00365
# 4 Full Model        LR (glm) 0.934   0.0288 
# 5 Demographics      XGB      0.752  NA      
# 6 Clinical Severity XGB      0.902   0.150  
# 7 Laboratories      XGB      0.896  -0.00583
# 8 Full Model        XGB      0.928   0.0317 

# combine all predictions
all_preds <- bind_rows(preds_01, preds_02, preds_03, preds_04)

# roc_curve() returns sensitivity and specificity correctly
roc_data <- all_preds %>%
  group_by(model, feature_set) %>%
  roc_curve(truth = inhosp_death_90day, .pred_no) %>%  # Use .pred_no
  ungroup() %>%
  mutate(fpr = 1 - specificity)

fig1 = ggplot(roc_data, aes(x = fpr, y = sensitivity, color = model)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray40", linewidth = 0.5) +
  facet_wrap(~feature_set, nrow = 2) +
  scale_color_manual(
    values = c("LR" = "#E64B35", "XGB" = "#4DBBD5"),
    labels = c("LR" = "Logistic Regression", "XGB" = "XGBoost")
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  labs(
    title = "ROC Curves: Mortality Prediction Across Feature Sets",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

ggsave(make_output_path('Fig',1,'ROCcurves','png'), 
       fig1, width = 6, height = 5, dpi = 300)
