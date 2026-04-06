# %%
import pandas as pd
import numpy as np

import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, root_mean_squared_error
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt

import joblib

import shap

import GEOparse





# %% [markdown]
# # DNAm Clock - GSE125105

# %%
file_path = r"C:\Users\serly\Documents\Bioinformatics\Learning\Jupyter\DNA Methylation\GSE125105_matrix_normalized.txt.gz"
df = pd.read_csv(file_path, sep="\t", compression="gzip", comment='!', low_memory=False)
# print(a, b, sep="\t") prints a and b separated by a tab space rather than the default single space.

# %%
pval_cols = [ c for c in df.columns if "DetectionPval" in c]
beta_cols = [ c for c in df.columns if "DetectionPval" not in c]

pval_samples = set(c.replace("_DetectionPval", "") for c in pval_cols)
beta_samples = set(beta_cols)

print("In beta but NOT in pval:", beta_samples - pval_samples)
print("In pval but NOT in beta:", pval_samples - beta_samples)

# {'ID_REF'} is the extra columns


# Spliting the data

pval_cols = [ c for c in df.columns if "DetectionPval" in c]
beta_cols = [ c for c in df.columns if "DetectionPval" not in c and c!= "ID_REF"]
# c!= "ID_REF" to filter out the ID_REF

# Use ID_REF as the index
beta_df = df.set_index("ID_REF")[beta_cols].copy()
pval_df = df.set_index("ID_REF")[pval_cols].copy()
pval_df.columns = beta_df.columns  # now 699 == 699 

print(f"beta_df: {beta_df.shape}")   # (474239, 699)
print(f"pval_df: {pval_df.shape}")   # (474239, 699)



array = beta_df.values.copy() 
#to change from pandas df to numpy
# copy() -> to create a copy so modification will not affect the original df

probe_means = np.nanmean(array, axis=1) # 1D array

row_idx, col_idx = np.where(np.isnan(array))
array[row_idx, col_idx] = probe_means[row_idx]

beta_df = pd.DataFrame(array, index=beta_df.index, columns=beta_df.columns)

print(f"NaNs remaining: {np.isnan(array).sum()}")

X = beta_df.T
print(f"Final X shape: {X.shape}")

# %%
# Getting the y (target) from another supplementary data

# Load the data from GEOparse
gse = GEOparse.get_GEO("GSE125105", destdir="./", silent=True)

pheno = gse.phenotype_data


# Extract the age

# ── Extract age and sample name from pheno ────────────────────────────────────
pheno["age"] = pd.to_numeric(pheno["characteristics_ch1.1.age"], errors="coerce")

# Extract sample name from title e.g. "genomic DNA from sample291" → "sample291"
pheno["sample_name"] = pheno["title"].str.extract(r"(sample\d+)")

print(pheno[["sample_name", "age"]].head(10))
print(f"\nSamples with valid age: {pheno['age'].notna().sum()}")
print(pheno["age"].describe())

# ── Align with X ──────────────────────────────────────────────────────────────
pheno_indexed = pheno.set_index("sample_name")

# Find common samples -> this is like the intersection from data in X and y
common = X.index.intersection(pheno_indexed.index)
print(f"\nSamples in both X and pheno: {len(common)}")

X_aligned = X.loc[common]
y = pheno_indexed.loc[common, "age"]


# %%
# ML Elastic Net

x = X_aligned
y = y

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# training
model = ElasticNet(alpha=0.1, l1_ratio=0.9)

model_EN = model.fit(x_train, y_train)

y_pred_train = model_EN.predict(x_train)

MAE_train = mean_absolute_error(y_train, y_pred_train)
MSE_train = mean_squared_error(y_train, y_pred_train)
RMSE_train = root_mean_squared_error(y_train, y_pred_train)
R2_train = r2_score(y_train, y_pred_train)

print("MAE_train:{:.3f}".format(MAE_train))
print("MSE_train:{:.3f}".format(MSE_train))
print("RMSE_train:{:.3f}".format(RMSE_train))
print("R2_train:{:.3f}".format(R2_train))

# Test on test data

y_pred_test = model_EN.predict(x_test)

MAE_test = mean_absolute_error(y_test, y_pred_test)
MSE_test = mean_squared_error(y_test, y_pred_test)
RMSE_test = root_mean_squared_error(y_test, y_pred_test)
R2_test = r2_score(y_test, y_pred_test)

print("MAE_test:{:.3f}".format(MAE_test))
print("MSE_test:{:.3f}".format(MSE_test))
print("RMSE_test:{:.3f}".format(RMSE_test))
print("R2_test:{:.3f}".format(R2_test))



