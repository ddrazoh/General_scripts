import random
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np

def simulate_sequence(length=300):
    """Generate a random DNA sequence."""
    return ''.join(random.choices("ACGT", k=length))

def introduce_errors(seq, error_rate=0.02):
    """Introduce random sequencing errors (substitutions, insertions, deletions)."""
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            mutation_type = random.choice(['sub', 'ins', 'del'])
            if mutation_type == 'sub':
                seq_list[i] = random.choice("ACGT".replace(seq_list[i], ''))
            elif mutation_type == 'ins':
                seq_list.insert(i, random.choice("ACGT"))
            elif mutation_type == 'del':
                seq_list[i] = ''
    return ''.join(seq_list)

def pairwise_comparison(seq1, seq2):
    """Perform global alignment between two sequences."""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    score = alignments[0].score / max(len(seq1), len(seq2))  # Normalize by length
    return score

# Simulate sequences for 100 TFs
num_samples = 100
sanger_sequences = [simulate_sequence() for _ in range(num_samples)]
ont_sequences = [introduce_errors(seq) for seq in sanger_sequences]

# Perform pairwise comparison
scores = [pairwise_comparison(sanger_sequences[i], ont_sequences[i]) for i in range(num_samples)]

# Plot sequence similarity scores
plt.hist(scores, bins=10, alpha=0.7, color='blue', edgecolor='black')
plt.xlabel("Similarity Score (ONT vs Sanger)")
plt.ylabel("Frequency")
plt.title("Pairwise Sequence Similarity")
plt.show()
