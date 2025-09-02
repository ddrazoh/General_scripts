
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
import numpy as np

# Read the CSV file with explicit column names
df = pd.read_csv('/Users/drake/Desktop/tf_trial/analysis/analysis2/combined_sequences.csv', delimiter='\t', names=['Sample_ID', 'MinION1_ErrorRate', 'MinION2_ErrorRate', 'MinION1_SNPs', 'MinION1_Indels', 'MinION2_SNPs', 'MinION2_Indels', 'MinION1_PairwiseIdentity', 'MinION2_PairwiseIdentity'])

# Calculate descriptive statistics
minion1_error_mean = df['MinION1_ErrorRate'].mean()
minion2_error_mean = df['MinION2_ErrorRate'].mean()
minion1_error_std = df['MinION1_ErrorRate'].std()
minion2_error_std = df['MinION2_ErrorRate'].std()

# Perform paired t-test for MinION1_ErrorRate and MinION2_ErrorRate
t_stat, p_value = ttest_rel(df['MinION1_ErrorRate'], df['MinION2_ErrorRate'])

# Create box plots for MinION1_ErrorRate and MinION2_ErrorRate
plt.figure(figsize=(8, 6))
sns.boxplot(data=df[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.xlabel('MinION')
plt.ylabel('Error Rate')
plt.title('MinION1 vs MinION2 Error Rate')
plt.show()

# Create violin plots for MinION1_ErrorRate and MinION2_ErrorRate
plt.figure(figsize=(8, 6))
sns.violinplot(data=df[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.xlabel('MinION')
plt.ylabel('Error Rate')
plt.title('MinION1 vs MinION2 Error Rate')
plt.show()

# Create a heatmap for pairwise correlations
correlation_matrix = df[['MinION1_ErrorRate', 'MinION2_ErrorRate', 'MinION1_SNPs', 'MinION1_Indels', 'MinION2_SNPs', 'MinION2_Indels', 'MinION1_PairwiseIdentity', 'MinION2_PairwiseIdentity']].corr()
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', square=True)
plt.title('Correlation Heatmap')
plt.show()

# Print the calculated statistics and test results
print('Descriptive Statistics:')
print(f"MinION1 Error Rate: Mean={minion1_error_mean:.4f}, Std={minion1_error_std:.4f}")
print(f"MinION2 Error Rate: Mean={minion2_error_mean:.4f}, Std={minion2_error_std:.4f}\n")
print('Paired t-test:')
print(f"t-statistic: {t_stat:.4f}")
print(f"p-value: {p_value:.4f}")
