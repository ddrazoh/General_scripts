from Bio import SeqIO
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.patches as mpatches
import pandas as pd

# === Input folder and files ===
folder = "/Users/drake/Documents/MRC_UVRI/MINION/RUNS/Analysis/All_3_prime/TF_alignments/*.fasta"
file_list = glob.glob(folder)

# === Colors for bases ===
base_colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red', '-': 'grey'}

# === First, collect sequences and compute similarity for all samples ===
samples_data = []

for filepath in file_list:
    sample_name = os.path.basename(filepath).split('.')[0]
    seqs = {rec.id: str(rec.seq).replace("\n", "").upper() for rec in SeqIO.parse(filepath, "fasta")}

    ont_seq = None
    sanger_seq = None
    for name, seq in seqs.items():
        if "ont" in name.lower():
            ont_seq = seq
        elif "sanger" in name.lower():
            sanger_seq = seq

    if ont_seq is None or sanger_seq is None:
        print(f"Skipping {sample_name} — missing ONT or Sanger sequence")
        continue

    seq_length = min(len(ont_seq), len(sanger_seq))
    matches = sum(1 for o, s in zip(ont_seq[:seq_length], sanger_seq[:seq_length]) if o == s)
    similarity = round(matches / seq_length * 100, 2)

    samples_data.append({
        "sample_name": sample_name,
        "ont_seq": ont_seq[:seq_length],
        "sanger_seq": sanger_seq[:seq_length],
        "seq_length": seq_length,
        "similarity": similarity
    })

# === Sort samples by similarity (highest first) ===
samples_data.sort(key=lambda x: x['similarity'], reverse=True)

# === Create figure ===
num_samples = len(samples_data)
fig, ax = plt.subplots(figsize=(15, 1.5 * num_samples), facecolor="white")
ax.set_facecolor("white")

tick_size = 0.15
y_pos = 0
similarity_scores = []
sample_labels = []

# === Plotting loop ===
for info in samples_data:
    ont_seq = info['ont_seq']
    sanger_seq = info['sanger_seq']
    seq_length = info['seq_length']
    sample_name = info['sample_name']

    # Grey backbone lines
    ax.hlines(y=y_pos + 0.2, xmin=0, xmax=seq_length-1, color='lightgrey', linewidth=2)
    ax.hlines(y=y_pos - 0.2, xmin=0, xmax=seq_length-1, color='lightgrey', linewidth=2)

    # ONT mismatches
    for i, base in enumerate(ont_seq):
        if base != sanger_seq[i]:
            ax.vlines(x=i,
                      ymin=y_pos + 0.2 - tick_size,
                      ymax=y_pos + 0.2 + tick_size,
                      color=base_colors.get(base.upper(), 'black'),
                      linewidth=1.2)

    # Sanger mismatches
    for i, base in enumerate(sanger_seq):
        if base != ont_seq[i]:
            ax.vlines(x=i,
                      ymin=y_pos - 0.2 - tick_size,
                      ymax=y_pos - 0.2 + tick_size,
                      color=base_colors.get(base.upper(), 'black'),
                      linewidth=1.2)

    # Labels
    ax.text(-150, y_pos, sample_name, ha='right', va='center', fontsize=7)
    ax.text(-100, y_pos + 0.2, "ONT", va='center', fontsize=6)
    ax.text(-150, y_pos - 0.2, "Sanger", va='center', fontsize=6)

    similarity_scores.append(info['similarity'])
    sample_labels.append(y_pos)
    y_pos -= 1

# === Similarity column ===
x_fixed_data = samples_data[0]['seq_length'] + 110
for sim, ypos in zip(similarity_scores, sample_labels):
    ax.text(x_fixed_data, ypos, f"{sim:.2f}%", va='center', ha='center', fontsize=7, color='black')

# Vertical label
x_label = samples_data[0]['seq_length'] + 230
y_center = (min(sample_labels) + max(sample_labels)) / 2
ax.text(x_label, y_center, "Similarity %", va='center', ha='center', rotation=90, fontsize=10)

# Formatting
ax.set_xlabel("Alignment Position", labelpad=10)
ax.set_ylabel("Samples", labelpad=75)
ax.set_title("Highlighter plot: ONT vs Sanger derived T/F sequences")
ax.set_yticks([])
ax.set_xlim(-10, samples_data[0]['seq_length'] + 10)

# Legend
legend_patches = [mpatches.Patch(color=color, label=base) for base, color in base_colors.items()]
ax.legend(handles=legend_patches,
          title="Mismatch Base",
          title_fontsize=8,
          bbox_to_anchor=(0.1, -0.05),
          loc='upper center',
          ncol=len(base_colors),
          handlelength=0.1,
          handleheight=0.5,
          borderaxespad=0.5,
          fontsize=8)

plt.subplots_adjust(bottom=0.18, right=0.85)
plt.savefig("highlighter_plot_sorted.png", dpi=1000, bbox_inches="tight", pad_inches=0.3, facecolor="white")
plt.savefig("highlighter_plot_sorted.pdf", bbox_inches="tight", pad_inches=0.3, facecolor="white")
plt.show()


stats_list = []

for info in samples_data:
    ont_seq = info['ont_seq']
    sanger_seq = info['sanger_seq']
    sample_name = info['sample_name']

    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0

    for o, s in zip(ont_seq, sanger_seq):
        if o == s:
            matches += 1
        elif o == '-':
            deletions += 1
        elif s == '-':
            insertions += 1
        else:
            mismatches += 1

    total_ref = len(sanger_seq)
    error_rate = (mismatches + insertions + deletions) / total_ref * 100

    stats_list.append({
        "Sample": sample_name,
        "Length": total_ref,
        "Matches": matches,
        "Mismatches(SNPs)": mismatches,
        "Insertions": insertions,
        "Deletions": deletions,
        "ErrorRate(%)": round(error_rate, 2),
        "Similarity(%)": info['similarity']
    })

# Convert to DataFrame
stats_df = pd.DataFrame(stats_list)
stats_df.to_csv("ONT_vs_Sanger_stats.csv", index=False)
print(stats_df.head())


plt.figure(figsize=(10,6))
plt.bar(stats_df['Sample'], stats_df['ErrorRate(%)'], color='orange')
plt.xticks(rotation=90)
plt.ylabel('Error Rate (%)')
plt.title('ONT Sequencing Error Rate vs Sanger Reference')
plt.tight_layout()
plt.show()

plt.figure(figsize=(12,6))
plt.bar(stats_df['Sample'], stats_df['Matches'], label='Matches', color='green')
plt.bar(stats_df['Sample'], stats_df['Mismatches(SNPs)'], bottom=stats_df['Matches'], label='Mismatches', color='red')
plt.bar(stats_df['Sample'], stats_df['Insertions'], bottom=stats_df['Matches'] + stats_df['Mismatches(SNPs)'], label='Insertions', color='blue')
plt.bar(stats_df['Sample'], stats_df['Deletions'], bottom=stats_df['Matches'] + stats_df['Mismatches(SNPs)'] + stats_df['Insertions'], label='Deletions', color='grey')

plt.xticks(rotation=90)
plt.ylabel('Number of bases')
plt.title('ONT vs Sanger Base Differences')
plt.legend()
plt.tight_layout()
plt.show()