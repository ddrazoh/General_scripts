import pandas as pd
import numpy as np
from Bio import pairwise2, SeqIO
from Bio.Align.Applications import MafftCommandline
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Step 1: Read the data from the combined CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/5p_transformed_data.csv')

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
for i, (sample_id, sanger_seq, minion_seq1, minion_seq2) in enumerate(
    zip(sample_ids, sanger_sequences, minion_sequences1, minion_sequences2)
):
    # Check the orientation of the sequences and adjust if needed
    if sanger_seq != minion_seq1:
        minion_seq1 = str(Seq(minion_seq1).reverse_complement())

    if sanger_seq != minion_seq2:
        minion_seq2 = str(Seq(minion_seq2).reverse_complement())

    # Perform pairwise alignment using MAFFT with auto-adjustment
    mafft_cline = MafftCommandline()
    input_seqs = [SeqRecord(Seq(sanger_seq), id="Sanger"), SeqRecord(Seq(minion_seq1), id="MinION1")]
    mafft_input = SeqIO.write(input_seqs, "alignment.fasta", "fasta")
    mafft_cline.input = "alignment.fasta"
    mafft_cline.auto = True
    stdout, stderr = mafft_cline()

    aligned_seqs = stdout.split("\n")
    aligned_sanger_seq1 = aligned_seqs[1].upper()
    aligned_minion_seq1 = aligned_seqs[3].upper()

    # Perform pairwise alignment for MinION 2
    mafft_cline = MafftCommandline()
    input_seqs = [SeqRecord(Seq(sanger_seq), id="Sanger"), SeqRecord(Seq(minion_seq2), id="MinION2")]
    mafft_input = SeqIO.write(input_seqs, "alignment.fasta", "fasta")
    mafft_cline.input = "alignment.fasta"
    mafft_cline.auto = True
    stdout, stderr = mafft_cline()

    aligned_seqs = stdout.split("\n")
    aligned_sanger_seq2 = aligned_seqs[1].upper()
    aligned_minion_seq2 = aligned_seqs[3].upper()

    # Calculate the error rate and pairwise identity for MinION 1
    seq_length = len(sanger_seq)
    mismatches1 = sum(a != b for a, b in zip(aligned_sanger_seq1, aligned_minion_seq1))
    indels1 = aligned_minion_seq1.count("-")
    error_rate1 = (mismatches1 + indels1) / seq_length * 100
    pairwise_identity1 = (seq_length - mismatches1 - indels1) / seq_length * 100

    # Calculate the error rate and pairwise identity for MinION 2
    mismatches2 = sum(a != b for a, b in zip(aligned_sanger_seq2, aligned_minion_seq2))
    indels2 = aligned_minion_seq2.count("-")
    error_rate2 = (mismatches2 + indels2) / seq_length * 100
    pairwise_identity2 = (seq_length - mismatches2 - indels2) / seq_length * 100

    minion_error_rates1.append(error_rate1)
    minion_error_rates2.append(error_rate2)
    pairwise_identities1.append(pairwise_identity1)
    pairwise_identities2.append(pairwise_identity2)

    # Print the error rates and pairwise identities for each sample
    print(f"Sample ID: {sample_id}")
    print(f"MinION 1 Error Rate: {error_rate1:.2f}%")
    print(f"MinION 2 Error Rate: {error_rate2:.2f}%")
    print(f"MinION 1 Pairwise Identity: {pairwise_identity1:.2f}%")
    print(f"MinION 2 Pairwise Identity: {pairwise_identity2:.2f}%")
    print()

# Calculate the overall mean error rates, standard deviations, and ranges for each MinION dataset
mean_error_rate1 = np.mean(minion_error_rates1)
std_error_rate1 = np.std(minion_error_rates1)
error_rate_range1 = np.max(minion_error_rates1) - np.min(minion_error_rates1)

mean_error_rate2 = np.mean(minion_error_rates2)
std_error_rate2 = np.std(minion_error_rates2)
