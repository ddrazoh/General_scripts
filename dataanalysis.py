import pandas as pd
from Bio import pairwise2
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from io import BytesIO
from pptx import Presentation
from pptx.util import Inches
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# Step 1: Read the data from the CSV file
data = pd.read_csv('/Users/drake/Desktop/tf_trial/analysis/3p_transformed_data.csv')

# Extract the sample IDs, Sanger sequences, and MinION sequences from the dataset
sample_ids = data['Sample_id']
sanger_sequences = data['Sanger']
minion_sequences = data['MinION']

# Step 2: Calculate Error Rates using pairwise alignment
minion_error_rates = []
pairwise_identities = []

for i, (sanger_seq, minion_seq) in enumerate(zip(sanger_sequences, minion_sequences)):
    alignment = pairwise2.align.globalms(sanger_seq, minion_seq, 2, -1, -2, -2)  # Perform pairwise alignment
    aligned_sanger_seq, aligned_minion_seq, _, _, _ = alignment[0]

    matches = sum(a == b for a, b in zip(aligned_sanger_seq, aligned_minion_seq))
    seq_length = len(aligned_sanger_seq)

    mismatches = seq_length - matches
    indels = aligned_minion_seq.count('-') + aligned_sanger_seq.count('-')

    error_rate = (mismatches + indels) / seq_length * 100
    pairwise_identity = matches / seq_length * 100

    minion_error_rates.append(error_rate)
    pairwise_identities.append(pairwise_identity)

    # Print the error rate and pairwise identity for each sample
    print(f"Sample ID: {sample_ids[i]}")
    print(f"MinION Error Rate: {error_rate:.2f}%")
    print(f"Pairwise Identity: {pairwise_identity:.2f}%")
    print()

# Calculate statistics
mean_error_rate = np.mean(minion_error_rates)
median_error_rate = np.median(minion_error_rates)
std_error_rate = np.std(minion_error_rates)
q1_error_rate = np.percentile(minion_error_rates, 25)
q3_error_rate = np.percentile(minion_error_rates, 75)
error_rate_range = np.ptp(minion_error_rates)

# Step 3: Statistical Analysis
t_statistic, p_value = stats.ttest_1samp(minion_error_rates, mean_error_rate)

# Confidence Interval Calculation
confidence_interval = stats.t.interval(0.95, len(minion_error_rates) - 1, loc=mean_error_rate, scale=stats.sem(minion_error_rates))

# Print the statistical results
print("Overall Mean MinION Error Rate:", f"{mean_error_rate:.2f}%")
print("Median MinION Error Rate:", f"{median_error_rate:.2f}%")
print("Standard Deviation:", f"{std_error_rate:.2f}")
print("1st Quartile (Q1):", f"{q1_error_rate:.2f}%")
print("3rd Quartile (Q3):", f"{q3_error_rate:.2f}%")
print("MinION Error Rate Range:", f"{error_rate_range:.2f}%")
print("T-test p-value:", p_value)
print("Confidence Interval:", confidence_interval)

# Create a DataFrame for the statistics
table_data = {'Sample ID': sample_ids,
              'MinION Error Rate (%)': minion_error_rates}
statistics
