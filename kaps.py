import argparse
from Bio import pairwise2, SeqIO

def parse_fasta(fasta_file):
    """Reads two sequences from a FASTA file."""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) != 2:
        raise ValueError("FASTA file must contain exactly two sequences.")
    return str(sequences[0].seq), str(sequences[1].seq)

def align_sequences(seq1, seq2):
    """Align two sequences using global alignment."""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]  # Get the best alignment (highest score)
    aligned_seq1, aligned_seq2 = best_alignment[:2]
    return aligned_seq1, aligned_seq2

def analyze_alignment(aligned_seq1, aligned_seq2):
    """Calculate the number of gaps and mismatches in the alignment."""
    gaps = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-')
    mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
    return gaps, mismatches

def main():
    parser = argparse.ArgumentParser(description="Align two sequences in a FASTA file and calculate gaps and mismatches.")
    parser.add_argument("input_file", help="Input FASTA file containing two sequences.")
    args = parser.parse_args()

    # Parse sequences from the FASTA file
    seq1, seq2 = parse_fasta(args.input_file)

    # Perform sequence alignment
    aligned_seq1, aligned_seq2 = align_sequences(seq1, seq2)

    # Analyze the alignment
    gaps, mismatches = analyze_alignment(aligned_seq1, aligned_seq2)

    # Output the results
    print("Alignment Results:")
    print(f"Aligned Sequence 1: {aligned_seq1}")
    print(f"Aligned Sequence 2: {aligned_seq2}")
    print(f"Gaps: {gaps}")
    print(f"Mismatches: {mismatches}")

if __name__ == "__main__":
    main()
