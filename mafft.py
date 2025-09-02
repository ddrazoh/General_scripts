import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
import scipy.stats as stats

# Step 1: Read the data from the CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/combined_sequences.csv')

# Step 2: Perform sequence alignment
def align_sequences(sanger_seq, minion_seq):
    # Check for missing values
    if pd.isnull(sanger_seq) or pd.isnull(minion_seq):
        return None, None

    sanger_seq = str(sanger_seq)
    minion_seq = str(minion_seq)

    alignments = pairwise2.align.globalxx(sanger_seq, minion_seq, gap_char="-")
    best_alignment = max(alignments, key=lambda x: x.score)
    aligned_sanger_seq, aligned_minion_seq = best_alignment[:2]
    return aligned_sanger_seq, aligned_minion_seq

# Step 3: Calculate error rates and pairwise identities
def calculate_error_rate(sanger_seq, minion_seq):
    # Check for missing values
    if pd.isnull(sanger_seq) or pd.isnull(minion_seq):
        return np.nan, np.nan

    seq_length = len(sanger_seq)
    mismatches = sum(a != b for a, b in zip(sanger_seq, minion_seq))
    indels = minion_seq.count('-')
    error_rate = (mismatches + indels) / seq_length * 100
    pairwise_identity = (seq_length - mismatches - indels) / seq_length * 100
    return error_rate, pairwise_identity

sample_ids = data['Sample_id']
sanger_sequences = data['Sanger']
minion_sequences1 = data['MinION1']
minion_sequences2 = data['MinION2']

minion_error_rates1 = []
minion_error_rates2 = []
pairwise_identities1 = []
pairwise_identities2 = []

for sanger_seq, minion_seq1, minion_seq2 in zip(sanger_sequences, minion_sequences1, minion_sequences2):
    aligned_sanger_seq1, aligned_minion_seq1 = align_sequences(sanger_seq, minion_seq1)
    aligned_sanger_seq2, aligned_minion_seq2 = align_sequences(sanger_seq, minion_seq2)

    if aligned_sanger_seq1 is None or aligned_minion_seq1 is None:
        continue

    if aligned_sanger_seq2 is None or aligned_minion_seq2 is None:
        continue
    
    error_rate1, pairwise_identity1 = calculate_error_rate(aligned_sanger_seq1, aligned_minion_seq1)
    error_rate2, pairwise_identity2 = calculate_error_rate(aligned_sanger_seq2, aligned_minion_seq2)
    
    minion_error_rates1.append(error_rate1)
    minion_error_rates2.append(error_rate2)
    pairwise_identities1.append(pairwise_identity1)
    pairwise_identities2.append(pairwise_identity2)

# Step 4: Analyze and compare the results
mean_error_rate1 = np.mean(minion_error_rates1)
std_error_rate1 = np.std(minion_error_rates1)
error_rate_range1 = np.max(minion_error_rates1) - np.min(minion_error_rates1)

mean_error_rate2 = np.mean(minion_error_rates2)
std_error_rate2 = np.std(minion_error_rates2)
error_rate_range2 = np.max(minion_error_rates2) - np.min(minion_error_rates2)

# Perform statistical tests (Mann-Whitney U test)
statistic, p_value = stats.mannwhitneyu(minion_error_rates1, minion_error_rates2)

# Step 5: Save the results to a CSV file
result_df = pd.DataFrame({
    'Sample_id': sample_ids,
    'MinION1 Error Rate': minion_error_rates1,
    'MinION2 Error Rate': minion_error_rates2,
    'MinION1 Pairwise Identity': pairwise_identities1,
    'MinION2 Pairwise Identity': pairwise_identities2
})

result_df.to_csv('validation_results.csv', index=False)

statistics_df = pd.DataFrame({
    'MinION Dataset': ['MinION1', 'MinION2'],
    'Mean Error Rate': [mean_error_rate1, mean_error_rate2],
    'Standard Deviation': [std_error_rate1, std_error_rate2],
    'Error Rate Range': [error_rate_range1, error_rate_range2]
})

statistics_df.to_csv('validation_statistics.csv', index=False)

print(f"Statistical test (Mann-Whitney U test):\nStatistic: {statistic}\np-value: {p_value}")
