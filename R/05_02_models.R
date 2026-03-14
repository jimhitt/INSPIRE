library(tidyverse)
library(tidymodels)
library(tidyverse)
library(DBI)
library(duckdb)
library(gt)
library(gtsummary)
library(scales)
tidymodels_prefer()
set.seed(2026)

data_dir = '/Volumes/ResearchDATA_v1/INSPIRE/data'
model_data = readRDS(file.path(data_dir,"PHI_modelData_v04.rds"))

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

model_data_clean2 = model_data_clean2 %>%
  mutate(any_adverse_event = replace_na(any_adverse_event, "no"))

# ---------------------------
# Shared helpers for modeling
# ---------------------------
set.seed(2026)

mort_metrics <- metric_set(roc_auc, pr_auc, brier_class)
ctrl <- control_resamples(save_pred = TRUE)

prep_outcome <- function(df, outcome = "any_adverse_event") {
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
                             outcome = "any_adverse_event",
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
                              outcome = "any_adverse_event",
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
  prep_outcome("any_adverse_event") %>%
  type_model_data_clean()  # <- from earlier snippet

# 2) recreate folds on the finalized data
set.seed(2026)
cv_folds2 <- vfold_cv(model_data_clean2, v = 5, strata = any_adverse_event)

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
  "preop_creatinine", "preop_sodium", "preop_potassium",
  "baseline_egfr"
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
  "preop_creatinine", "preop_sodium", "preop_potassium",
  "baseline_egfr"
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

# Resutls



# Make ROC plot -----------------------------------------------------------
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

# combine all predictions
all_preds <- bind_rows(preds_01, preds_02, preds_03, preds_04)

# roc_curve() returns sensitivity and specificity correctly
roc_data <- all_preds %>%
  group_by(model, feature_set) %>%
  roc_curve(truth = any_adverse_event, .pred_no) %>%  # Use .pred_no
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
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

ggsave(make_output_path('Fig',1,'ROCcurves','png'), 
       fig1, width = 6, height = 5, dpi = 300)


# Confusion Matrix --------------------------------------------------------

library(probably)

# Pull the full-feature predictions
preds_04_clean <- preds_04 %>%
  select(model, .pred_yes, any_adverse_event)

# Youden's J threshold per model
youden_thresholds <- preds_04_clean %>%
  group_by(model) %>%
  threshold_perf(any_adverse_event, .pred_yes, 
                 thresholds = seq(0.01, 0.99, by = 0.01),
                 metrics = metric_set(j_index)) %>%
  group_by(model) %>%
  slice_max(.estimate, n = 1) %>%
  select(model, youden_thresh = .threshold)

# Function to compute confusion matrix stats at a given threshold
confusion_stats <- function(df, threshold, threshold_name) {
  df %>%
    mutate(
      .pred_class = factor(if_else(.pred_yes >= threshold, "yes", "no"),
                           levels = c("no", "yes"))
    ) %>%
    group_by(model) %>%
    summarise(
      TP = sum(.pred_class == "yes" & any_adverse_event == "yes"),
      TN = sum(.pred_class == "no"  & any_adverse_event == "no"),
      FP = sum(.pred_class == "yes" & any_adverse_event == "no"),
      FN = sum(.pred_class == "no"  & any_adverse_event == "yes"),
      .groups = "drop"
    ) %>%
    mutate(
      threshold_name = threshold_name,
      threshold      = threshold,
      sensitivity    = TP / (TP + FN),
      specificity    = TN / (TN + FP),
      ppv            = TP / (TP + FP),
      npv            = TN / (TN + FN),
      fpr            = FP / (FP + TN)
    )
}

# Apply all three thresholds
results_prevalence <- confusion_stats(preds_04_clean, 0.05, "Prevalence (0.05)")

youden_lr  <- youden_thresholds %>% filter(model == "LR")  %>% pull(youden_thresh)
youden_xgb <- youden_thresholds %>% filter(model == "XGB") %>% pull(youden_thresh)

results_youden <- bind_rows(
  confusion_stats(preds_04_clean %>% filter(model == "LR"),  youden_lr,  "Youden's J"),
  confusion_stats(preds_04_clean %>% filter(model == "XGB"), youden_xgb, "Youden's J")
)

results_sensitivity <- confusion_stats(preds_04_clean, 
                                       # threshold that gets ~80% sensitivity for each model
                                       quantile(preds_04_clean$.pred_yes[
                                         preds_04_clean$any_adverse_event == "yes"], 0.20),
                                       "80% Sensitivity")

confusion_summary <- bind_rows(results_prevalence, results_youden, results_sensitivity) %>%
  arrange(model, threshold_name)

print(confusion_summary)
# This is one of the most important teaching moments in clinical prediction modeling, and your data illustrates it perfectly.
# The Youden's J result is telling. At the optimal threshold by that criterion, both models essentially refuse to predict any 
# events — LR calls 5 patients positive out of 97K, XGBoost calls 8. The "optimal" threshold mathematically is nearly 1.0 
# because with 5% prevalence, the penalty for false positives is so asymmetric that the model learns to be extremely 
# conservative. Youden's J is poorly suited to imbalanced outcomes.
# The core tension is PPV vs sensitivity. At 80% sensitivity (LR), you're flagging 71,000 patients to catch 5,000 
# real events — a PPV of 7%. In a clinical workflow that means 13 alerts for every 1 true event. That's textbook 
# alert fatigue territory. XGBoost handles this better at the same threshold (PPV 36%, FPR 6%) because its probability 
# calibration is different, but sensitivity drops to 61%.
# The prevalence threshold is probably the most honest operating point for a clinical audience. LR at 0.05 catches 77% of 
# events with a PPV of 20% — meaning 1 in 5 flagged patients truly had an adverse outcome. Whether that's useful depends 
# entirely on what the intervention is and what it costs to act on a false positive.
# This is actually not a failure of the model — it's a fundamental property of predicting rare events. A ROC-AUC 
# of 0.905 means the model ranks patients well, but ranking well and having actionable thresholds are different 
# things. The blog post should make this distinction explicit, because it's where a lot of clinical AI tools 
# oversell their utility.


# Create Table For Blog ---------------------------------------------------

# ============================================================
# Confusion Matrix Analysis — Threshold Selection
# ============================================================
#
# The XGBoost model achieves an AUROC of 0.905, which sounds
# impressive. But AUROC measures ranking ability — how well the
# model orders patients by risk — not clinical utility. The
# confusion matrices below illustrate why a high AUROC does not
# guarantee a useful clinical tool.
#
# The fundamental challenge is threshold selection. A predicted
# probability is continuous; a clinical decision is binary. The
# threshold that converts one to the other determines the entire
# character of the tool:
#
# THRESHOLD A — High Sensitivity (0.010):
#   Catches 61% of true events, but generates 5,520 false
#   positives (PPV 36%). For every patient correctly flagged,
#   roughly 2 are flagged unnecessarily. In a busy surgical
#   service, this volume of false alerts would quickly erode
#   clinician trust — a phenomenon known as alert fatigue.
#
# THRESHOLD B — Prevalence-Based (0.050):
#   A principled middle ground: flag patients whose predicted
#   risk exceeds the population event rate of 5.2%. Sensitivity
#   drops to 39% — we are now missing most events — but the
#   false positive burden falls to 1,550 (PPV 56%). Whether
#   this tradeoff is acceptable depends entirely on the cost
#   and feasibility of the intended intervention.
#
# THRESHOLD C — Youden's J (0.990):
#   Mathematically "optimal" in that it maximizes the sum of
#   sensitivity and specificity. In practice, it almost never
#   predicts an event (sensitivity 0.2%, PPV 100%). At 5%
#   prevalence, the algorithm learns that the safest strategy
#   is near-silence. High specificity dominates because true
#   negatives vastly outnumber true positives. Youden's J is
#   poorly suited to imbalanced outcomes.
#
# The core lesson: no threshold is universally correct. The
# right operating point depends on the clinical context —
# specifically, the relative costs of missing a true event
# (false negative) versus acting on a false alarm (false
# positive). A screening tool tolerates false positives;
# a tool that triggers an invasive intervention cannot.
#
# AUROC tells you the model can rank patients. Threshold
# analysis tells you whether that ranking translates into
# something actionable.
# ============================================================

library(htmltools)

# XGBoost rows only, in clinical order
xgb_confusion <- confusion_summary %>%
  filter(model == "XGB") %>%
  arrange(factor(threshold_name, levels = c("80% Sensitivity", 
                                            "Prevalence (0.05)", 
                                            "Youden's J")))

# Build a single annotated confusion matrix as a gt table
make_confusion_gt <- function(row_data, title, subtitle) {
  
  TP <- row_data$TP;  TN <- row_data$TN
  FP <- row_data$FP;  FN <- row_data$FN
  
  sens <- sprintf("%.1f%%", row_data$sensitivity * 100)
  spec <- sprintf("%.1f%%", row_data$specificity * 100)
  ppv  <- sprintf("%.1f%%", row_data$ppv * 100)
  npv  <- sprintf("%.1f%%", row_data$npv * 100)
  
  total_actual_pos  <- TP + FN
  total_actual_neg  <- TN + FP
  total_pred_neg    <- TN + FN
  total_pred_pos    <- TP + FP
  total_n           <- TP + TN + FP + FN
  
  fmt <- function(x) format(x, big.mark = ",")
  
  tibble(
    `Actual Outcome`     = c("Event (Yes)",  "No Event (No)",  "Total Predicted", "Metric"),
    `Predicted: No Event` = c(fmt(FN),        fmt(TN),          fmt(total_pred_neg), paste0("NPV = ", npv)),
    `Predicted: Event`    = c(fmt(TP),        fmt(FP),          fmt(total_pred_pos), paste0("PPV = ", ppv)),
    `Total Actual`        = c(fmt(total_actual_pos), fmt(total_actual_neg), fmt(total_n), ""),
    `Metric`              = c(paste0("Sensitivity = ", sens), paste0("Specificity = ", spec), "", "")
  ) %>%
    gt(rowname_col = "Actual Outcome") %>%
    tab_header(title = title, subtitle = subtitle) %>%
    # Color the 4 cells of the confusion matrix
    tab_style(
      style = cell_fill(color = "#d4edda"),   # TP: green
      locations = cells_body(columns = `Predicted: Event`, rows = 1)
    ) %>%
    tab_style(
      style = cell_fill(color = "#f8d7da"),   # FN: red (missed events)
      locations = cells_body(columns = `Predicted: No Event`, rows = 1)
    ) %>%
    tab_style(
      style = cell_fill(color = "#f8d7da"),   # FP: red (false alarms)
      locations = cells_body(columns = `Predicted: Event`, rows = 2)
    ) %>%
    tab_style(
      style = cell_fill(color = "#d4edda"),   # TN: green
      locations = cells_body(columns = `Predicted: No Event`, rows = 2)
    ) %>%
    # Style margin rows
    tab_style(
      style = list(cell_fill(color = "#f2f2f2"), cell_text(weight = "bold")),
      locations = cells_body(rows = c(3, 4))
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_stub()
    ) %>%
    # Spanners to clarify column groups
    tab_spanner(
      label    = "Predicted Outcome",
      columns  = c(`Predicted: No Event`, `Predicted: Event`)
    ) %>%
    tab_spanner(
      label    = "Marginal Statistics",
      columns  = c(`Total Actual`, `Metric`)
    ) %>%
    opt_table_outline() %>%
    tab_options(
      table.width                  = pct(100),
      data_row.padding             = px(8),
      column_labels.font.weight    = "600",
      table.font.size              = px(14),
      stub.font.weight             = "bold"
    )
}

# Define titles and subtitles for each threshold
thresh_specs <- list(
  list(
    name     = "80% Sensitivity",
    title    = "Threshold A: High Sensitivity (threshold = 0.010)",
    subtitle = "Catches 61% of events — but flags 5,520 patients who will not have a complication"
  ),
  list(
    name     = "Prevalence (0.05)",
    title    = "Threshold B: Prevalence-Based (threshold = 0.050)",
    subtitle = "Flag patients whose predicted risk exceeds the population event rate of 5.2%"
  ),
  list(
    name     = "Youden's J",
    title    = "Threshold C: Youden's J — Mathematically Optimal (threshold = 0.990)",
    subtitle = "Maximizes sensitivity + specificity simultaneously — but nearly never predicts an event"
  )
)

# Build gt tables
gt_tables <- lapply(thresh_specs, function(x) {
  row <- xgb_confusion %>% filter(threshold_name == x$name)
  make_confusion_gt(row, x$title, x$subtitle)
})

# Combine into a single HTML file using htmltools
combined_html <- htmltools::tagList(
  htmltools::tags$h2("XGBoost Model — Full Feature Set (AUROC = 0.905)"),
  htmltools::tags$p(
    "Three operating thresholds illustrating the tension between sensitivity, 
     specificity, and clinical utility. Green cells indicate correct classifications; 
     red cells indicate errors. PPV = positive predictive value; NPV = negative predictive value."
  ),
  htmltools::HTML(as_raw_html(gt_tables[[1]])),
  htmltools::tags$br(),
  htmltools::HTML(as_raw_html(gt_tables[[2]])),
  htmltools::tags$br(),
  htmltools::HTML(as_raw_html(gt_tables[[3]]))
)

htmltools::save_html(
  combined_html,
  make_output_path("Table", 2, "ConfusionMatrices", "html")
)
