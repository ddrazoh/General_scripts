import pandas as pd
from Bio import pairwise2
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Read the data from the CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/5p_transformed_data.csv')

# Extract the sample IDs, Sanger sequences, and MinION sequences from the dataset
sample_ids = data['Sample_id']
sanger_sequences = data['Sanger']
minion_sequences = data['MinION']

# Initialize lists to store error rates and pairwise identities
minion_error_rates = []
pairwise_identities = []

# Step 2: Calculate Error Rates and Pairwise Identities using pairwise alignment
for i, (sanger_seq, minion_seq) in enumerate(zip(sanger_sequences, minion_sequences)):
    alignment = pairwise2.align.globalms(sanger_seq, minion_seq, 2, -1, -2, -2)  # Perform pairwise alignment
    aligned_sanger_seq, aligned_minion_seq, _, _, _ = alignment[0]

    mismatches = 0
    indels = 0
    seq_length = len(aligned_sanger_seq)

    for j in range(seq_length):
        if aligned_sanger_seq[j] != aligned_minion_seq[j]:
            if aligned_sanger_seq[j] == '-' or aligned_minion_seq[j] == '-':
                indels += 1
            else:
                mismatches += 1

    error_rate = (mismatches + indels) / seq_length * 100
    pairwise_identity = (seq_length - mismatches - indels) / seq_length * 100

    minion_error_rates.append(error_rate)
    pairwise_identities.append(pairwise_identity)

    # Print the error rate and pairwise identity for each sample
    print(f"Sample ID: {sample_ids[i]}")
    print(f"MinION Error Rate: {error_rate:.2f}%")
    print(f"Pairwise Identity: {pairwise_identity:.2f}%")
    print()

# Calculate the overall mean error rate, standard deviation, and range
mean_error_rate = np.mean(minion_error_rates)
std_error_rate = np.std(minion_error_rates)
error_rate_range = np.max(minion_error_rates) - np.min(minion_error_rates)

# Calculate the overall mean pairwise identity
mean_pairwise_identity = np.mean(pairwise_identities)

# Step 3: Statistical Analysis
t_statistic, p_value = stats.ttest_ind(minion_error_rates, np.ones_like(minion_error_rates) * mean_error_rate)

# Print the statistical results for error rate
print("Overall Mean MinION Error Rate:", f"{mean_error_rate:.2f}%")
print("MinION Error Rate Standard Deviation:", f"{std_error_rate:.2f}%")
print("MinION Error Rate Range:", f"{error_rate_range:.2f}%")
print("T-test p-value (Error Rate):", p_value)

# Create a figure and axes
fig, ax = plt.subplots(figsize=(10, 6))

# Plot MinION Error Rate histogram with KDE overlay
sns.histplot(minion_error_rates, bins=10, kde=True, ax=ax)
ax.set_xlabel('MinION Error Rate (%)')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of MinION Error Rates with KDE Overlay')

# Save the figure as an image
plt.savefig('histogram.png')

# Create a violin plot
fig, ax = plt.subplots(figsize=(10, 6))
sns.violinplot(data=minion_error_rates, ax=ax)
ax.set_ylabel('MinION Error Rate (%)')
ax.set_title('Violin Plot of MinION Error Rates')

# Save the figure as an image
plt.savefig('violin_plot.png')

# Create a box plot
fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=minion_error_rates, ax=ax)
ax.set_ylabel('MinION Error Rate (%)')
ax.set_title('Box Plot of MinION Error Rates')

# Save the figure as an image
plt.savefig('box_plot.png')

# Show the plots
plt.show()
