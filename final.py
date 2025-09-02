import pandas as pd
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, wilcoxon, chi2_contingency
import scipy.stats



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

    alignment1 = pairwise2.align.globalms(sanger_seq, minion_seq1, 2, -1, -1, -0.5)
    aligned_sanger_seq1, aligned_minion_seq1, _, _, _ = alignment1[0]

    alignment2 = pairwise2.align.globalms(sanger_seq, minion_seq2, 2, -1, -1, -0.5)
    aligned_sanger_seq2, aligned_minion_seq2, _, _, _ = alignment2[0]

    snps1 = 0
    indels1 = 0

    snps2 = 0
    indels2 = 0

    if len(aligned_sanger_seq2) == len(aligned_minion_seq2):
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
    print(f"MinION 1 Error Rate: {error_rate1}%")
    print(f"MinION 2 Error Rate: {error_rate2}%")
    print(f"MinION 1 Pairwise Identity: {pairwise_identity1}%")
    print(f"MinION 2 Pairwise Identity: {pairwise_identity2}%")
    print()

# Step 2: Create a new DataFrame to store the results
results = pd.DataFrame({'Sample_ID': sample_ids,
                        'MinION1_ErrorRate': minion_error_rates1,
                        'MinION2_ErrorRate': minion_error_rates2,
                        'MinION1_SNPs': snps_counts1,
                        'MinION1_Indels': indels_counts1,
                        'MinION2_SNPs': snps_counts2,
                        'MinION2_Indels': indels_counts2,
                        'MinION1_PairwiseIdentity': pairwise_identities1,
                        'MinION2_PairwiseIdentity': pairwise_identities2})

# Step 3: Save the results to a CSV file
results.to_csv('error_snp_mismatch_results.csv', index=False)

# Step 3: Perform statistical tests
# Wilcoxon signed-rank test for error rates
wilcoxon_error_rates = wilcoxon(results['MinION1_ErrorRate'], results['MinION2_ErrorRate'])
print("Wilcoxon Signed-Rank Test Results:")
print(f"Error Rates - p-value: {wilcoxon_error_rates.pvalue}")

ttest_error_rates = ttest_ind(results['MinION1_ErrorRate'], results['MinION2_ErrorRate'])
print("T-test Results:")
print(f"Error Rates - t-statistic: {ttest_error_rates.statistic}")
print(f"Error Rates - p-value: {ttest_error_rates.pvalue}")

# McNemar's test for SNPs
contingency_table = pd.crosstab(results['MinION1_SNPs'], results['MinION2_SNPs'])
mcnemar_snps = chi2_contingency(contingency_table)
print("McNemar's Test Results (SNPs):")
print(f"Chi-square statistic: {mcnemar_snps[0]}")
print(f"p-value: {mcnemar_snps[1]}")

# McNemar's test for indels
contingency_table = pd.crosstab(results['MinION1_Indels'], results['MinION2_Indels'])
mcnemar_indels = chi2_contingency(contingency_table)
print("McNemar's Test Results (Indels):")
print(f"Chi-square statistic: {mcnemar_indels[0]}")
print(f"p-value: {mcnemar_indels[1]}")

minion1_error_rates = results['MinION1_ErrorRate']
minion2_error_rates = results['MinION2_ErrorRate']

# Standard deviation
minion1_std = np.std(minion1_error_rates)
minion2_std = np.std(minion2_error_rates)

# Variance
minion1_var = np.var(minion1_error_rates)
minion2_var = np.var(minion2_error_rates)

print("Measures of Dispersion:")
print(f"MinION 1 Standard Deviation: {minion1_std}")
print(f"MinION 2 Standard Deviation: {minion2_std}")
print(f"MinION 1 Variance: {minion1_var}")
print(f"MinION 2 Variance: {minion2_var}")

# Measures of Central Tendency
minion1_mean = np.mean(minion1_error_rates)
minion2_mean = np.mean(minion2_error_rates)
minion1_median = np.median(minion1_error_rates)
minion2_median = np.median(minion2_error_rates)

print("Measures of Central Tendency:")
print(f"MinION 1 Mean: {minion1_mean}")
print(f"MinION 2 Mean: {minion2_mean}")
print(f"MinION 1 Median: {minion1_median}")
print(f"MinION 2 Median: {minion2_median}")

# Confidence Interval
confidence_level = 0.95

# MinION 1 Confidence Interval
minion1_error_rates_ci = scipy.stats.t.interval(confidence_level, len(minion1_error_rates)-1,
                                                loc=np.mean(minion1_error_rates),
                                                scale=scipy.stats.sem(minion1_error_rates))
minion1_error_rates_ci_lower, minion1_error_rates_ci_upper = minion1_error_rates_ci

print(f"MinION 1 Confidence Interval ({confidence_level * 100}%):")
print(f"Lower Bound: {minion1_error_rates_ci_lower}")
print(f"Upper Bound: {minion1_error_rates_ci_upper}")

# MinION 2 Confidence Interval
minion2_error_rates_ci = scipy.stats.t.interval(confidence_level, len(minion2_error_rates)-1,
                                                loc=np.mean(minion2_error_rates),
                                                scale=scipy.stats.sem(minion2_error_rates))
minion2_error_rates_ci_lower, minion2_error_rates_ci_upper = minion2_error_rates_ci

print(f"MinION 2 Confidence Interval ({confidence_level * 100}%):")
print(f"Lower Bound: {minion2_error_rates_ci_lower}")
print(f"Upper Bound: {minion2_error_rates_ci_upper}")

# Measures of Dispersion
minion1_error_rates = results['MinION1_ErrorRate']
minion2_error_rates = results['MinION2_ErrorRate']

# Standard deviation
minion1_std = np.std(minion1_error_rates)
minion2_std = np.std(minion2_error_rates)

# Variance
minion1_var = np.var(minion1_error_rates)
minion2_var = np.var(minion2_error_rates)

print("Measures of Dispersion:")
print(f"MinION 1 Standard Deviation: {minion1_std}")
print(f"MinION 2 Standard Deviation: {minion2_std}")
print(f"MinION 1 Variance: {minion1_var}")
print(f"MinION 2 Variance: {minion2_var}")

# Measures of Central Tendency
minion1_mean = np.mean(minion1_error_rates)
minion2_mean = np.mean(minion2_error_rates)
minion1_median = np.median(minion1_error_rates)
minion2_median = np.median(minion2_error_rates)

print("Measures of Central Tendency:")
print(f"MinION 1 Mean: {minion1_mean}")
print(f"MinION 2 Mean: {minion2_mean}")
print(f"MinION 1 Median: {minion1_median}")
print(f"MinION 2 Median: {minion2_median}")

# Confidence Interval
confidence_level = 0.95

# MinION 1 Confidence Interval
minion1_error_rates_ci = scipy.stats.t.interval(confidence_level, len(minion1_error_rates)-1,
                                                loc=np.mean(minion1_error_rates),
                                                scale=scipy.stats.sem(minion1_error_rates))
minion1_error_rates_ci_lower, minion1_error_rates_ci_upper = minion1_error_rates_ci

print(f"MinION 1 Confidence Interval ({confidence_level * 100}%):")
print(f"Lower Bound: {minion1_error_rates_ci_lower}")
print(f"Upper Bound: {minion1_error_rates_ci_upper}")

# MinION 2 Confidence Interval
minion2_error_rates_ci = scipy.stats.t.interval(confidence_level, len(minion2_error_rates)-1,
                                                loc=np.mean(minion2_error_rates),
                                                scale=scipy.stats.sem(minion2_error_rates))
minion2_error_rates_ci_lower, minion2_error_rates_ci_upper = minion2_error_rates_ci

print(f"MinION 2 Confidence Interval ({confidence_level * 100}%):")
print(f"Lower Bound: {minion2_error_rates_ci_lower}")
print(f"Upper Bound: {minion2_error_rates_ci_upper}")

# Step 4: Plot the distributions of error rates, SNPs, and indels
plt.figure(figsize=(10, 6))
sns.histplot(data=results, x='MinION1_ErrorRate', kde=True, label='MinION 1')
sns.histplot(data=results, x='MinION2_ErrorRate', kde=True, label='MinION 2')
plt.xlabel('Error Rate (%)')
plt.ylabel('Count')
plt.title('Distribution of MinION Error Rates')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
sns.histplot(data=results, x='MinION1_SNPs', kde=True, label='MinION 1')
sns.histplot(data=results, x='MinION2_SNPs', kde=True, label='MinION 2')
plt.xlabel('SNPs/Mismatches')
plt.ylabel('Count')
plt.title('Distribution of MinION SNPs/Mismatches')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
sns.histplot(data=results, x='MinION1_Indels', kde=True, label='MinION 1')
sns.histplot(data=results, x='MinION2_Indels', kde=True, label='MinION 2')
plt.xlabel('Indels')
plt.ylabel('Count')
plt.title('Distribution of MinION Indels')
plt.legend()
plt.show()

# Additional plots: Violin plots and Heatmap
plt.figure(figsize=(10, 6))
sns.violinplot(data=results[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.xlabel('MinION Dataset')
plt.ylabel('Error Rate (%)')
plt.title('Error Rate Distribution - MinION 1 vs MinION 2')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.show()

plt.figure(figsize=(10, 6))
sns.heatmap(data=results[['MinION1_PairwiseIdentity', 'MinION2_PairwiseIdentity']], cmap='YlOrRd')
plt.xlabel('MinION Dataset')
plt.ylabel('Sample ID')
plt.title('Pairwise Identity Comparison - MinION 1 vs MinION 2')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.show()

# Heat map for error rates
error_rates_heatmap = results.pivot_table(values='MinION1_ErrorRate', index='Sample_ID', columns='MinION2_ErrorRate')
plt.figure(figsize=(10, 8))
sns.heatmap(error_rates_heatmap, cmap='YlOrRd', annot=True, fmt=".1f", cbar=True)
plt.xlabel('MinION 2 Error Rate (%)')
plt.ylabel('MinION 1 Error Rate (%)')
plt.title('Error Rate Heat Map')
plt.show()

# Scatter plot for error rates and pairwise identities
plt.figure(figsize=(10, 6))
plt.scatter(results['MinION1_ErrorRate'], results['MinION1_PairwiseIdentity'], label='MinION 1')
plt.scatter(results['MinION2_ErrorRate'], results['MinION2_PairwiseIdentity'], label='MinION 2')
plt.xlabel('Error Rate (%)')
plt.ylabel('Pairwise Identity (%)')
plt.title('Error Rates vs Pairwise Identities')
plt.legend()
plt.show()

# Line graph for error rates and pairwise identities
plt.figure(figsize=(10, 6))
plt.plot(results.index, results['MinION1_ErrorRate'], label='MinION 1', marker='o')
plt.plot(results.index, results['MinION2_ErrorRate'], label='MinION 2', marker='o')
plt.plot(results.index, results['MinION1_PairwiseIdentity'], label='MinION 1 Pairwise Identity', linestyle='--', marker='o')
plt.plot(results.index, results['MinION2_PairwiseIdentity'], label='MinION 2 Pairwise Identity', linestyle='--', marker='o')
plt.xlabel('Sample')
plt.ylabel('Error Rate / Pairwise Identity')
plt.title('Error Rates and Pairwise Identities')
plt.xticks(results.index, results['Sample_ID'], rotation=45)
plt.legend()
plt.show()

# Grouped violin plot for Indels and SNPs
plt.figure(figsize=(10, 6))
sns.violinplot(data=results[['MinION1_Indels', 'MinION2_Indels', 'MinION1_SNPs', 'MinION2_SNPs']], inner='quartile', cut=0)
plt.xlabel('Pipeline')
plt.ylabel('Counts')
plt.title('Distribution of Indels and SNPs by Pipeline')
plt.xticks([0, 1, 2, 3], ['MinION 1 Indels', 'MinION 2 Indels', 'MinION 1 SNPs', 'MinION 2 SNPs'])
plt.show()

plt.figure(figsize=(10, 6))
sns.violinplot(data=results[['MinION1_Indels', 'MinION1_SNPs']], inner='quartile', cut=0, split=True, linewidth=0.8)
sns.violinplot(data=results[['MinION2_Indels', 'MinION2_SNPs']], inner='quartile', cut=0, split=True, linewidth=0.8)
plt.xlabel('Pipeline')
plt.ylabel('Counts')
plt.title('Distribution of Indels and SNPs by Pipeline')
plt.xticks([0, 1], ['MinION 1', 'MinION 2'])
plt.legend(['Indels', 'SNPs'])
plt.show()


# Data for the split violins
minion1_indels = results['MinION1_Indels']
minion1_snps = results['MinION1_SNPs']
minion2_indels = results['MinION2_Indels']
minion2_snps = results['MinION2_SNPs']

# Position of the violins
positions = np.array([0, 0.5])

# Create the figure and axes
plt.figure(figsize=(10, 6))
ax = plt.gca()

# Plot the violins
ax.violinplot([minion1_indels, minion2_indels], positions=positions-0.1, widths=0.2, showmedians=True)
ax.violinplot([minion1_snps, minion2_snps], positions=positions+0.1, widths=0.2, showmedians=True)

# Set the labels and title
plt.xlabel('Pipeline')
plt.ylabel('Counts')
plt.title('Distribution of Indels and SNPs by Pipeline')

# Set the x-axis ticks and labels
ax.set_xticks(positions)
ax.set_xticklabels(['MinION 1', 'MinION 2'])

# Set the legend
indels_patch = plt.Line2D([0], [0], color='blue', lw=4)
snps_patch = plt.Line2D([0], [0], color='orange', lw=4)
ax.legend([indels_patch, snps_patch], ['Indels', 'SNPs'])

# Show the plot
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

# Set the figure size
plt.figure(figsize=(10, 6))

# Create the split violin plot for MinION 1
sns.violinplot(data=results[['MinION1_Indels', 'MinION1_SNPs']], inner='quartile', cut=0, split=True,
               linewidth=0.8, palette=['blue', 'red'])

# Create the split violin plot for MinION 2
sns.violinplot(data=results[['MinION2_Indels', 'MinION2_SNPs']], inner='quartile', cut=0, split=True,
               linewidth=0.8, palette=['green', 'orange'])

# Set the labels and title
plt.xlabel('Pipeline')
plt.ylabel('Counts')
plt.title('Distribution of Indels and SNPs by Pipeline')

# Create the custom legend
indels_patch = plt.Line2D([0], [0], color='blue', linewidth=8)
snps_patch = plt.Line2D([0], [0], color='red', linewidth=8)
legend = plt.legend([indels_patch, snps_patch], ['Indels', 'SNPs'], loc='upper right')

# Show the plot
plt.show()


plt.figure(figsize=(10, 6))

# Create the split violin plot for MinION 1
sns.violinplot(data=results[['MinION1_Indels', 'MinION1_SNPs']].stack().reset_index(level=1), hue='level_1', y=0,
               inner='quartile', cut=0, split=True, linewidth=0.8, palette=['purple', 'orange'])

# Create the split violin plot for MinION 2
sns.violinplot(data=results[['MinION2_Indels', 'MinION2_SNPs']].stack().reset_index(level=1), hue='level_1', y=0,
               inner='quartile', cut=0, split=True, linewidth=0.8, palette=['blue', 'green'])

# Set the labels and title
plt.xlabel('Pipeline')
plt.ylabel('Counts')
plt.title('Distribution of Indels and SNPs by Pipeline')

# Create the custom legend
indels_patch = plt.Line2D([0], [0], color='purple', linewidth=8)
snps_patch = plt.Line2D([0], [0], color='orange', linewidth=8)
legend = plt.legend([indels_patch, snps_patch], ['Indels', 'SNPs'], loc='upper right')

# Show the plot
plt.show()





