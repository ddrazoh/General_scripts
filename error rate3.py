import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Step 1: Load the data
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/5p_transformed_data.csv')

# Step 2: Process the data
sample_ids = data['Sample_id']
sanger_seqs = data['Sanger']
minion_seqs1 = data['MinION1']
minion_seqs2 = data['MinION2']

seq_length = len(sanger_seqs[0])  # Assuming all sequences have the same length

minion_error_rates1 = []
minion_error_rates2 = []
pairwise_identities1 = []
pairwise_identities2 = []
mismatch_counts1 = []
mismatch_counts2 = []
indel_counts1 = []
indel_counts2 = []

total_sequences = len(sample_ids)
progress_interval = max(total_sequences // 10, 1)  # Print progress every 10% or for small datasets, every sequence

for i, (sanger_seq, minion_seq1, minion_seq2) in enumerate(zip(sanger_seqs, minion_seqs1, minion_seqs2), 1):
    aligned_sanger_seq1 = sanger_seq.replace('-', '')
    aligned_minion_seq1 = minion_seq1.replace('-', '')
    aligned_sanger_seq2 = sanger_seq.replace('-', '')
    aligned_minion_seq2 = minion_seq2.replace('-', '')

    # Calculate the mismatch and indel counts for MinION 1
    mismatches1 = sum(a != b for a, b in zip(aligned_sanger_seq1, aligned_minion_seq1))
    indels1 = aligned_minion_seq1.count("-")
    mismatch_counts1.append(mismatches1)
    indel_counts1.append(indels1)

    # Calculate the mismatch and indel counts for MinION 2
    mismatches2 = sum(a != b for a, b in zip(aligned_sanger_seq2, aligned_minion_seq2))
    indels2 = aligned_minion_seq2.count("-")
    mismatch_counts2.append(mismatches2)
    indel_counts2.append(indels2)

    # Calculate the error rate and pairwise identity for MinION 1
    error_rate1 = (mismatches1 + indels1) / seq_length * 100
    pairwise_identity1 = (seq_length - mismatches1 - indels1) / seq_length * 100

    # Calculate the error rate and pairwise identity for MinION 2
    error_rate2 = (mismatches2 + indels2) / seq_length * 100
    pairwise_identity2 = (seq_length - mismatches2 - indels2) / seq_length * 100

    # Append the error rates and pairwise identities to the respective lists
    minion_error_rates1.append(error_rate1)
    minion_error_rates2.append(error_rate2)
    pairwise_identities1.append(pairwise_identity1)
    pairwise_identities2.append(pairwise_identity2)

    # Print progress
    if i % progress_interval == 0 or i == total_sequences:
        print(f"Processing sequence {i}/{total_sequences}")

# Step 3: Calculate statistics and save the results to a CSV file
result_data = pd.DataFrame({
    'Sample ID': sample_ids,
    'MinION Error Rate 1': minion_error_rates1,
    'MinION Error Rate 2': minion_error_rates2,
    'Pairwise Identity 1': pairwise_identities1,
    'Pairwise Identity 2': pairwise_identities2,
    'MinION Mismatch Count 1': mismatch_counts1,
    'MinION Mismatch Count 2': mismatch_counts2,
    'MinION Indel Count 1': indel_counts1,
    'MinION Indel Count 2': indel_counts2,
})

result_data['MinION Error Rate Mean 1'] = result_data['MinION Error Rate 1'].mean()
result_data['MinION Error Rate Mean 2'] = result_data['MinION Error Rate 2'].mean()
result_data['MinION Error Rate STD 1'] = result_data['MinION Error Rate 1'].std()
result_data['MinION Error Rate STD 2'] = result_data['MinION Error Rate 2'].std()

result_data['MinION Mismatch Count Mean 1'] = result_data['MinION Mismatch Count 1'].mean()
result_data['MinION Mismatch Count Mean 2'] = result_data['MinION Mismatch Count 2'].mean()
result_data['MinION Mismatch Count STD 1'] = result_data['MinION Mismatch Count 1'].std()
result_data['MinION Mismatch Count STD 2'] = result_data['MinION Mismatch Count 2'].std()

result_data['MinION Indel Count Mean 1'] = result_data['MinION Indel Count 1'].mean()
result_data['MinION Indel Count Mean 2'] = result_data['MinION Indel Count 2'].mean()
result_data['MinION Indel Count STD 1'] = result_data['MinION Indel Count 1'].std()
result_data['MinION Indel Count STD 2'] = result_data['MinION Indel Count 2'].std()

result_data.to_csv('result_data.csv', index=False)

# Step 4: Visualize the data
sns.set(style="whitegrid")
plt.figure(figsize=(12, 6))

# Boxplot
plt.subplot(1, 2, 1)
sns.boxplot(data=[minion_error_rates1, minion_error_rates2], palette=["#1f77b4", "#ff7f0e"])
plt.xticks([0, 1], ["MinION 1", "MinION 2"])
plt.xlabel("MinION Dataset")
plt.ylabel("Error Rate (%)")
plt.title("Error Rates Comparison")

# Profiled errors
plt.subplot(1, 2, 2)
plt.bar(np.arange(total_sequences), mismatch_counts1, color="#1f77b4", label="MinION 1 Mismatches")
plt.bar(np.arange(total_sequences), indel_counts1, bottom=mismatch_counts1, color="#ff7f0e", label="MinION 1 Indels")
plt.bar(np.arange(total_sequences), mismatch_counts2, color="#1f77b4", alpha=0.5, label="MinION 2 Mismatches")
plt.bar(np.arange(total_sequences), indel_counts2, bottom=mismatch_counts2, color="#ff7f0e", alpha=0.5, label="MinION 2 Indels")
plt.xlabel("Sequence")
plt.ylabel("Count")
plt.title("Profiled Errors")
plt.legend()

# Save the figures
plt.savefig('visualization.png')

plt.show()
