import pandas as pd
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Step 1: Read the data from the combined CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/combined_sequences.csv')

# Extract the sample IDs, Sanger sequences, MinION sequences 1, and MinION sequences 2 from the dataset
sample_ids = data['Sample_id']
sanger_sequences = data['Sanger']
minion_sequences1 = data['MinION1']
minion_sequences2 = data['MinION2']

# Initialize lists to store error rates, pairwise identities, SNPs/mismatches, and indels for each MinION dataset
minion_error_rates1 = []
minion_error_rates2 = []
snps_counts1 = []
snps_counts2 = []
indels_counts1 = []
indels_counts2 = []
pairwise_identities1 = []
pairwise_identities2 = []

for i, (sanger_seq, minion_seq1, minion_seq2) in enumerate(zip(sanger_sequences, minion_sequences1, minion_sequences2)):
    sanger_length = len(sanger_seq)
    minion_length1 = len(minion_seq1)
    minion_length2 = len(minion_seq2)
    seq_length = min(sanger_length, minion_length1, minion_length2)

    alignment1 = pairwise2.align.localms(sanger_seq, minion_seq1, 2, -1, -2, -2)
    aligned_sanger_seq1, aligned_minion_seq1, _, _, _ = alignment1[0]

    alignment2 = pairwise2.align.localms(sanger_seq, minion_seq2, 2, -1, -2, -2)
    aligned_sanger_seq2, aligned_minion_seq2, _, _, _ = alignment2[0]

    snps1 = 0
    indels1 = 0

    snps2 = 0
    indels2 = 0

    for j in range(seq_length):
        if aligned_sanger_seq1[j] != aligned_minion_seq1[j]:
            if aligned_sanger_seq1[j] == '-' or aligned_minion_seq1[j] == '-':
                indels1 += 1
            else:
                snps1 += 1

        if aligned_sanger_seq2[j] != aligned_minion_seq2[j]:
            if aligned_sanger_seq2[j] == '-' or aligned_minion_seq2[j] == '-':
                indels2 += 1
            else:
                snps2 += 1

    indels_counts1.append(indels1)
    indels_counts2.append(indels2)
    snps_counts1.append(snps1)
    snps_counts2.append(snps2)

    error_rate1 = (snps1 + indels1) / seq_length * 100
    pairwise_identity1 = (seq_length - snps1 - indels1) / seq_length * 100

    error_rate2 = (snps2 + indels2) / seq_length * 100
    pairwise_identity2 = (seq_length - snps2 - indels2) / seq_length * 100

    minion_error_rates1.append(error_rate1)
    minion_error_rates2.append(error_rate2)
    pairwise_identities1.append(pairwise_identity1)
    pairwise_identities2.append(pairwise_identity2)

    print(f"Sample ID: {sample_ids[i]}")
    print(f"MinION 1 SNPs/Mismatches: {snps1}")
    print(f"MinION 1 Indels: {indels1}")
    print(f"MinION 2 SNPs/Mismatches: {snps2}")
    print(f"MinION 2 Indels: {indels2}")
    print()

# Create a DataFrame to store the results
results = pd.DataFrame({
    'Sample ID': sample_ids,
    'MinION 1 SNPs/Mismatches': snps_counts1,
    'MinION 1 Indels': indels_counts1,
    'MinION 2 SNPs/Mismatches': snps_counts2,
    'MinION 2 Indels': indels_counts2,
    'MinION 1 Error Rate': minion_error_rates1,
    'MinION 2 Error Rate': minion_error_rates2,
    'MinION 1 Pairwise Identity': pairwise_identities1,
    'MinION 2 Pairwise Identity': pairwise_identities2
})

# Perform paired t-test between MinION 1 and MinION 2 error rates
t_statistic, p_value = stats.ttest_rel(minion_error_rates1, minion_error_rates2)
results['MinION 1 vs MinION 2 t-statistic'] = t_statistic
results['MinION 1 vs MinION 2 p-value'] = p_value

# Perform paired Wilcoxon signed-rank test between MinION 1 and MinION 2 error rates
statistic, p_value = stats.wilcoxon(minion_error_rates1, minion_error_rates2)
results['MinION 1 vs MinION 2 Wilcoxon Test Statistic'] = statistic
results['MinION 1 vs MinION 2 p-value (Wilcoxon)'] = p_value

# Save the results to a CSV file
results.to_csv('results.csv', index=False)

# Visualize the data

# Violin Plot for error rates
plt.figure(figsize=(8, 6))
sns.violinplot(data=results[['MinION 1 Error Rate', 'MinION 2 Error Rate']], palette='Set2')
plt.xlabel('MinION Dataset')
plt.ylabel('Error Rate')
plt.title('Error Rates: MinION 1 vs MinION 2')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.savefig('error_rate_violinplot.png')
plt.show()

# Violin Plot for SNPs/mismatches
plt.figure(figsize=(8, 6))
sns.violinplot(data=results[['MinION 1 SNPs/Mismatches', 'MinION 2 SNPs/Mismatches']], palette='Set2')
plt.xlabel('MinION Dataset')
plt.ylabel('SNPs/Mismatches')
plt.title('SNPs/Mismatches: MinION 1 vs MinION 2')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.savefig('snps_violinplot.png')
plt.show()

# Violin Plot for indels
plt.figure(figsize=(8, 6))
sns.violinplot(data=results[['MinION 1 Indels', 'MinION 2 Indels']], palette='Set2')
plt.xlabel('MinION Dataset')
plt.ylabel('Indels')
plt.title('Indels: MinION 1 vs MinION 2')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.savefig('indels_violinplot.png')
plt.show()

# Visualize the data

# Grouped Bar Chart for error rates and SNPs/mismatches
x = np.arange(len(results))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))
error_rates = ax.bar(x - width/2, results['MinION 1 Error Rate'], width, label='MinION 1 Error Rate')
snps_counts = ax.bar(x + width/2, results['MinION 1 SNPs/Mismatches'], width, label='MinION 1 SNPs/Mismatches')

error_rates2 = ax.bar(x - width/2, results['MinION 2 Error Rate'], width, label='MinION 2 Error Rate')
snps_counts2 = ax.bar(x + width/2, results['MinION 2 SNPs/Mismatches'], width, label='MinION 2 SNPs/Mismatches')

ax.set_xlabel('Sample ID')
ax.set_ylabel('Counts / Error Rate')
ax.set_title('Error Rates and SNPs/Mismatches: MinION 1 vs MinION 2')
ax.set_xticks(x)
ax.set_xticklabels(results['Sample ID'], rotation=90)
ax.legend()

plt.tight_layout()
plt.savefig('error_rate_snps_barchart.png')
plt.show()

# Create a color gradient based on error rates or SNPs/mismatches
color_palette = sns.color_palette("viridis", len(results))

# Plot the grouped bar chart with color gradient
plt.figure(figsize=(10, 6))
sns.barplot(data=results[['MinION 1 Error Rate', 'MinION 2 Error Rate', 'MinION 1 SNPs/Mismatches', 'MinION 2 SNPs/Mismatches']],
            palette=color_palette)
plt.xlabel('Sample ID')
plt.ylabel('Count')
plt.title('Error Rates and SNPs/Mismatches: MinION 1 vs MinION 2')
plt.xticks(rotation=90)
plt.legend(['MinION 1 Error Rate', 'MinION 2 Error Rate', 'MinION 1 SNPs/Mismatches', 'MinION 2 SNPs/Mismatches'])
plt.show()

# Create a heatmap for error rates and SNPs/mismatches
plt.figure(figsize=(10, 8))
sns.heatmap(data=results[['MinION 1 Error Rate', 'MinION 2 Error Rate', 'MinION 1 SNPs/Mismatches', 'MinION 2 SNPs/Mismatches']],
            cmap='YlOrRd', annot=True, fmt=".1f", linewidths=0.5)
plt.xlabel('Metrics')
plt.ylabel('Sample ID')
plt.title('Error Rates and SNPs/Mismatches: MinION 1 vs MinION 2')
plt.show()


# Create a bubble chart for error rates and SNPs/mismatches
plt.figure(figsize=(10, 8))
sns.scatterplot(x='MinION 1 Error Rate', y='MinION 1 SNPs/Mismatches', size='Sample ID', data=results, alpha=0.8)
sns.scatterplot(x='MinION 2 Error Rate', y='MinION 2 SNPs/Mismatches', size='Sample ID', data=results, alpha=0.8)
plt.xlabel('Error Rate')
plt.ylabel('SNPs/Mismatches')
plt.title('Error Rates vs SNPs/Mismatches: MinION 1 vs MinION 2')
plt.legend(['MinION 1', 'MinION 2'])
plt.show()
