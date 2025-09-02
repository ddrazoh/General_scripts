import sys
import os
import logging
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re

logging.basicConfig(level=logging.INFO)

def read_primer_file(primer_file):
    if not os.path.isfile(primer_file):
        logging.error(f"Primer file '{primer_file}' not found.")
        sys.exit(1)

    with open(primer_file, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]  # Filter out empty lines
        logging.info(f"Read lines from primer file: {lines}")  # Debug line
        logging.info(f"Number of lines read: {len(lines)}")  # Debug line
        if len(lines) != 3:
            logging.error("Primer file must contain exactly three lines (5' primer, 3' primer, mismatch tolerance).")
            sys.exit(1)

        primer_5_prime, primer_3_prime, mismatch_tolerance = lines
        return primer_5_prime, primer_3_prime, int(mismatch_tolerance)

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def extract_amplicon(sequence, pattern_5_prime, pattern_3_prime):
    match_5_prime = re.search(pattern_5_prime, sequence)
    match_3_prime = re.search(pattern_3_prime, sequence)

    if match_5_prime and match_3_prime:
        start = match_5_prime.end()
        end = match_3_prime.start()
        if start < end:
            return sequence[start:end]
    return None

# Get input arguments
if len(sys.argv) != 3:
    logging.error("Usage: python3 extract_amplicons.py <input_fasta_file> <primers_file>")
    sys.exit(1)

input_fasta = sys.argv[1]
primer_file = sys.argv[2]

# Read primers and mismatch tolerance
primer_5_prime, primer_3_prime, mismatch_tolerance = read_primer_file(primer_file)

# Prepare regex patterns
pattern_5_prime = f"({primer_5_prime}){{e<={mismatch_tolerance}}}"
pattern_3_prime = f"({reverse_complement(primer_3_prime)}){{e<={mismatch_tolerance}}}"

output_fasta = "amplicon_sequences.fasta"
bin_fasta = "sequences_without_primers.fasta"

with open(output_fasta, "w") as out_handle, open(bin_fasta, "w") as bin_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        logging.info(f"Processing sequence {record.id}")

        # Ensure the sequence is in the 5' to 3' direction
        sequence_str = str(record.seq)
        oriented_seq = sequence_str
        
        # Use the reverse complement if the sequence doesn't match the expected 5' to 3' orientation
        if oriented_seq != reverse_complement(oriented_seq):
            oriented_seq = reverse_complement(sequence_str)

        amplicon_seq = extract_amplicon(oriented_seq, pattern_5_prime, pattern_3_prime)

        if amplicon_seq:
            record.seq = Seq(amplicon_seq)
            SeqIO.write(record, out_handle, "fasta")
        else:
            SeqIO.write(record, bin_handle, "fasta")
