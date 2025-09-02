import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.multivariate.manova import MANOVA

# Sample data for both pipelines (average pairwise identity and average error rate)
pipeline1_data = np.array([[99.5, 0.5]] * 96)
pipeline2_data = np.array([[98.6, 1.46]] * 96)

# Create a DataFrame with the data and group labels
data = np.vstack((pipeline1_data, pipeline2_data))
df = pd.DataFrame(data, columns=["Pairwise Identity", "Error Rate"])
df["Pipeline"] = np.array(["Pipeline 1"] * 96 + ["Pipeline 2"] * 96)

# Perform MANOVA
manova = MANOVA.from_formula('Q("Pairwise Identity") + Q("Error Rate") ~ Pipeline', data=df)

# Print the results
print("Multivariate Analysis of Variance (MANOVA):")
print(manova.mv_test())
