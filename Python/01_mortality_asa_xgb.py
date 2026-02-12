"""
Replicate R XGBoost model: ASA classification predicting 90-day mortality
Demonstrates R/Python equivalence for mortality prediction
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
import os
import sys

if len(sys.argv) > 1:
    data_path = sys.argv[1]
else:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(os.path.dirname(script_dir), 'data', 'PHI_model_data_v01.csv')

# Load data
print("Loading data...")
data = pd.read_csv(data_path)

# Prepare outcome and predictor
X = data[['asa']].values
y = (data['inhosp_death_90day'] == True).astype(int)

print(f"\nData loaded: {len(data):,} observations")
print(f"Events: {y.sum():,} ({y.mean()*100:.2f}%)")

# 5-fold stratified cross-validation
print("\nRunning 5-fold cross-validation...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

aurocs = []
fold_num = 1

for train_idx, test_idx in cv.split(X, y):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # Train XGBoost (match R parameters)
    model = XGBClassifier(
        n_estimators=100,
        max_depth=6,
        random_state=42,
        eval_metric='logloss'
    )
    model.fit(X_train, y_train)
    
    # Predict and calculate AUROC
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    auroc = roc_auc_score(y_test, y_pred_proba)
    aurocs.append(auroc)
    
    print(f"Fold {fold_num}: AUROC = {auroc:.4f}")
    fold_num += 1

# Summary statistics
print(f"\n{'='*50}")
print(f"Mean AUROC: {np.mean(aurocs):.4f}")
print(f"Std Dev:    {np.std(aurocs):.4f}")
print(f"Min:        {np.min(aurocs):.4f}")
print(f"Max:        {np.max(aurocs):.4f}")
print(f"{'='*50}")

# Compare to R result
r_auroc = 0.802
print(f"\nR XGBoost AUROC:      {r_auroc:.3f}")
print(f"Python XGBoost AUROC: {np.mean(aurocs):.3f}")
print(f"Difference:           {abs(np.mean(aurocs) - r_auroc):.3f}")

# Save results
results_df = pd.DataFrame({
    'fold': range(1, 6),
    'auroc': aurocs
})
output_path = os.path.join(os.path.dirname(os.path.dirname(data_path)), "output", "SAFE_python_xgb_results.csv")
results_df.to_csv(output_path, index=False)
print("\nResults saved to output/SAFE_python_xgb_results.csv")