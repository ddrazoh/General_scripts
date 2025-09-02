import pandas as pd
from Bio import pairwise2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Read the data from the combined CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/combined_sequences.csv')

# Extract the sample IDs, Sanger sequences, MinION sequences 1, and MinION sequences 2 from the dataset
sample_ids = data['Sample_id']
sanger_sequences = data['Sanger']
minion_sequences1 = data['MinION1']
minion_sequences2 = data['MinION2']

# Initialize lists to store error rates and pairwise identities for each MinION dataset
minion_error_rates1 = []
minion_error_rates2 = []
pairwise_identities1 = []
pairwise_identities2 = []

# Step 2: Calculate Error Rates and Pairwise Identities using pairwise alignment for each MinION dataset
for i, (sanger_seq, minion_seq1, minion_seq2) in enumerate(zip(sanger_sequences, minion_sequences1, minion_sequences2)):
    sanger_length = len(sanger_seq)
    minion_length1 = len(minion_seq1)
    minion_length2 = len(minion_seq2)
    seq_length = min(sanger_length, minion_length1, minion_length2)

    alignment1 = pairwise2.align.globalxx(sanger_seq, minion_seq1, one_alignment_only=True)
    aligned_sanger_seq1, aligned_minion_seq1, _, _, _ = alignment1[0]

    alignment2 = pairwise2.align.globalxx(sanger_seq, minion_seq2, one_alignment_only=True)
    aligned_sanger_seq2, aligned_minion_seq2, _, _, _ = alignment2[0]

    mismatches1 = 0
    indels1 = 0

    mismatches2 = 0
    indels2 = 0

    for j in range(seq_length):
        if aligned_sanger_seq1[j] != aligned_minion_seq1[j]:
            if aligned_sanger_seq1[j] == '-' or aligned_minion_seq1[j] == '-':
                indels1 += 1
            else:
                mismatches1 += 1

        if aligned_sanger_seq2[j] != aligned_minion_seq2[j]:
            if aligned_sanger_seq2[j] == '-' or aligned_minion_seq2[j] == '-':
                indels2 += 1
            else:
                mismatches2 += 1

    error_rate1 = (mismatches1 + indels1) / seq_length * 100
    pairwise_identity1 = (seq_length - mismatches1 - indels1) / seq_length * 100

    error_rate2 = (mismatches2 + indels2) / seq_length * 100
    pairwise_identity2 = (seq_length - mismatches2 - indels2) / seq_length * 100

    minion_error_rates1.append(error_rate1)
    minion_error_rates2.append(error_rate2)
    pairwise_identities1.append(pairwise_identity1)
    pairwise_identities2.append(pairwise_identity2)

    # Print the error rates and pairwise identities for each sample
    print(f"Sample ID: {sample_ids[i]}")
    print(f"MinION 1 Error Rate: {error_rate1:.2f}%")
    print(f"MinION 2 Error Rate: {error_rate2:.2f}%")
    print(f"MinION 1 Pairwise Identity: {pairwise_identity1:.2f}%")
    print(f"MinION 2 Pairwise Identity: {pairwise_identity2:.2f}%")
    print()

# Create a DataFrame to store the results
results = pd.DataFrame({
    'Sample ID': sample_ids,
    'MinION 1 Error Rate': minion_error_rates1,
    'MinION 2 Error Rate': minion_error_rates2,
    'MinION 1 Pairwise Identity': pairwise_identities1,
    'MinION 2 Pairwise Identity': pairwise_identities2
})

# Save the results to a CSV file
results.to_csv('error_identity_results.csv', index=False)

# Plot the error rates for MinION 1 and MinION 2
plt.figure(figsize=(10, 6))
sns.boxplot(data=results[['MinION 1 Error Rate', 'MinION 2 Error Rate']])
plt.xlabel('MinION Dataset')
plt.ylabel('Error Rate (%)')
plt.title('Box Plot of MinION Error Rates')
plt.savefig('error_rates_plot.png')
plt.show()

# Plot the pairwise identities for MinION 1 and MinION 2
plt.figure(figsize=(10, 6))
sns.boxplot(data=results[['MinION 1 Pairwise Identity', 'MinION 2 Pairwise Identity']])
plt.xlabel('MinION Dataset')
plt.ylabel('Pairwise Identity (%)')
plt.title('Box Plot of MinION Pairwise Identities')
plt.savefig('pairwise_identities_plot.png')
plt.show()
