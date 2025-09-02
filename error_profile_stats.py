from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import matplotlib.patches as mpatches

# === Input folder and files ===
folder = "/Users/drake/Documents/MRC_UVRI/MINION/RUNS/Analysis/All_3_prime/TF_alignments/*.fasta"
file_list = glob.glob(folder)

# === Colors for bases and error types ===
base_colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red', '-': 'grey'}
error_colors = {'match': 'lightgrey', 'SNP': 'red', 'insertion': 'blue', 'deletion': 'purple'}

stats_list = []
error_matrix = []

# === Parse sequences, compute stats ===
for filepath in file_list:
    sample_name = os.path.basename(filepath).split('.')[0]
    seqs = {rec.id: str(rec.seq).replace("\n", "").upper() for rec in SeqIO.parse(filepath, "fasta")}

    ont_seq = next((s for name, s in seqs.items() if "ont" in name.lower()), None)
    sanger_seq = next((s for name, s in seqs.items() if "sanger" in name.lower()), None)

    if ont_seq is None or sanger_seq is None:
        print(f"Skipping {sample_name} — missing ONT or Sanger sequence")
        continue

    seq_length = min(len(ont_seq), len(sanger_seq))
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    sample_error_row = []

    for o, s in zip(ont_seq[:seq_length], sanger_seq[:seq_length]):
        if o == s:
            matches += 1
            sample_error_row.append(0)  # 0 = match
        elif o == '-':
            deletions += 1
            sample_error_row.append(3)  # 3 = deletion
        elif s == '-':
            insertions += 1
            sample_error_row.append(2)  # 2 = insertion
        else:
            mismatches += 1
            sample_error_row.append(1)  # 1 = SNP

    total_ref = seq_length
    error_rate = (mismatches + insertions + deletions) / total_ref * 100
    similarity = matches / total_ref * 100

    stats_list.append({
        "Sample": sample_name,
        "Length": total_ref,
        "Matches": matches,
        "Mismatches(SNPs)": mismatches,
        "Insertions": insertions,
        "Deletions": deletions,
        "ErrorRate(%)": round(error_rate, 2),
        "Similarity(%)": round(similarity, 2)
    })
    error_matrix.append(sample_error_row)

# === Convert to DataFrame ===
stats_df = pd.DataFrame(stats_list).sort_values(by="Similarity(%)", ascending=False)
error_matrix = np.array(error_matrix)

# === Highlighter Plot ===
fig, ax = plt.subplots(figsize=(15, len(stats_df)*0.6))
for i, row in enumerate(error_matrix):
    y = len(stats_df)-1 - i
    for pos, val in enumerate(row):
        if val == 1:
            ax.plot(pos, y, 's', color=error_colors['SNP'], markersize=4)
        elif val == 2:
            ax.plot(pos, y, 's', color=error_colors['insertion'], markersize=4)
        elif val == 3:
            ax.plot(pos, y, 's', color=error_colors['deletion'], markersize=4)
ax.set_yticks(range(len(stats_df)))
ax.set_yticklabels(stats_df['Sample'])
ax.set_xlabel("Sequence Position")
ax.set_title("Highlighter Plot: ONT vs Sanger (SNPs, Insertions, Deletions)")
patches = [mpatches.Patch(color=color, label=key) for key, color in error_colors.items()]
ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# === Error rate bar plot ===
plt.figure(figsize=(12,6))
plt.bar(stats_df['Sample'], stats_df['ErrorRate(%)'], color='orange')
plt.xticks(rotation=90)
plt.ylabel('Error Rate (%)')
plt.title('ONT Sequencing Error Rate vs Sanger Reference')
plt.tight_layout()
plt.show()

# === Stacked base differences plot ===
plt.figure(figsize=(12,6))
plt.bar(stats_df['Sample'], stats_df['Matches'], label='Matches', color='green')
plt.bar(stats_df['Sample'], stats_df['Mismatches(SNPs)'], bottom=stats_df['Matches'], label='SNPs', color='red')
plt.bar(stats_df['Sample'], stats_df['Insertions'], bottom=stats_df['Matches']+stats_df['Mismatches(SNPs)'], label='Insertions', color='blue')
plt.bar(stats_df['Sample'], stats_df['Deletions'], bottom=stats_df['Matches']+stats_df['Mismatches(SNPs)']+stats_df['Insertions'], label='Deletions', color='purple')
plt.xticks(rotation=90)
plt.ylabel('Number of Bases')
plt.title('ONT vs Sanger Base Differences')
plt.legend()
plt.tight_layout()
plt.show()

# === Save stats to CSV ===
stats_df.to_csv("ONT_vs_Sanger_stats_sorted.csv", index=False)

print(stats_df.head())


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Convert stats/error matrix to long-format DataFrame for plotting
long_data = []

for i, row in enumerate(error_matrix):
    sample_name = stats_df.iloc[i]['Sample']
    for val in row:
        if val == 1:  # SNP
            long_data.append({'Sample': sample_name, 'ErrorType': 'SNP'})
        elif val == 3:  # Deletion
            long_data.append({'Sample': sample_name, 'ErrorType': 'Deletion'})

plot_df = pd.DataFrame(long_data)

# Violin plot
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Convert stats_df into long format
plot_df = stats_df.melt(
    id_vars=["Sample"],
    value_vars=["Mismatches(SNPs)", "Deletions"],
    var_name="ErrorType",
    value_name="Count"
)

plt.figure(figsize=(12,6))
sns.violinplot(
    x="Sample",
    y="Count",
    hue="ErrorType",
    data=plot_df,
    split=True,
    palette={"Mismatches(SNPs)": "red", "Deletions": "purple"}
)

plt.xticks(rotation=90)
plt.ylabel("Counts per sample")
plt.title("Distribution of SNPs and Deletions across Samples")
plt.tight_layout()

avg_similarity = stats_df["Similarity(%)"].mean()
min_similarity = stats_df["Similarity(%)"].min()
max_similarity = stats_df["Similarity(%)"].max()

print(f"Average similarity: {avg_similarity:.2f}%")
print(f"Similarity range: {min_similarity:.2f}% – {max_similarity:.2f}%")

plt.show()


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# --- Similarity stats ---
similarities = stats_df["Similarity(%)"]

avg_similarity = similarities.mean()
min_similarity = similarities.min()
max_similarity = similarities.max()

# 95% confidence interval using t-distribution
conf_int = stats.t.interval(
    0.95,                      # 95% CI
    len(similarities)-1,       # degrees of freedom
    loc=avg_similarity,        # mean
    scale=stats.sem(similarities)  # standard error of mean
)

print(f"Average similarity: {avg_similarity:.2f}%")
print(f"Similarity range: {min_similarity:.2f}% – {max_similarity:.2f}%")
print(f"95% CI for mean similarity: {conf_int[0]:.2f}% – {conf_int[1]:.2f}%")

import numpy as np
import scipy.stats as st

# Extract similarities
similarities = stats_df["Similarity(%)"].astype(float)

# Mean
mean_sim = np.mean(similarities)

# Range
min_sim = np.min(similarities)
max_sim = np.max(similarities)

# Sample size
n = len(similarities)

# Standard error of the mean
sem = st.sem(similarities)

# 95% CI using t-distribution
ci_low, ci_high = st.t.interval(0.95, df=n-1, loc=mean_sim, scale=sem)

print(f"Mean similarity: {mean_sim:.2f}%")
print(f"Range: {min_sim:.2f}% – {max_sim:.2f}%")
print(f"95% CI: {ci_low:.2f}% – {ci_high:.2f}%")


# --- Plot histogram + KDE ---
plt.figure(figsize=(10,6))
sns.histplot(similarities, bins=20, kde=True, color="green")
plt.axvline(avg_similarity, color="red", linestyle="--", label=f"Mean = {avg_similarity:.2f}%")
plt.axvspan(conf_int[0], conf_int[1], color="red", alpha=0.2, label="95% CI")
plt.xlabel("Similarity (%)")
plt.ylabel("Count")
plt.title("Distribution of ONT vs Sanger T/F Pairwise Sequence Similarities")
plt.legend()
plt.tight_layout()
plt.show()
