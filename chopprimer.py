import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
from concurrent.futures import ThreadPoolExecutor

# Function to read primers from a .txt file
def read_primer_file(primer_file):
    with open(primer_file, 'r') as file:
        primer_5_prime = file.readline().strip()  # 5' primer
        primer_3_prime = file.readline().strip()  # 3' primer
        mismatch_tolerance = int(file.readline().strip())
    return primer_5_prime, primer_3_prime, mismatch_tolerance

# Function to get the reverse complement of a DNA sequence
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# Function to process each sequence and extract amplicons
def process_sequence(record, pattern_5_prime, pattern_3_prime):
    sequence = str(record.seq)
    match_5_prime = pattern_5_prime.search(sequence)
    match_3_prime = pattern_3_prime.search(sequence)

    if match_5_prime and match_3_prime:
        amplicon_start = match_5_prime.end()
        amplicon_end = match_3_prime.start()
        if amplicon_start < amplicon_end:
            return record.id, sequence[amplicon_start:amplicon_end], True
    return record.id, sequence, False

# Get input arguments
if len(sys.argv) != 3:
    print("Usage: python3 extract_amplicons.py <input_fasta_file> <primers_file>")
    sys.exit(1)

input_fasta = sys.argv[1]
primer_file = sys.argv[2]

# Read primers and prepare regex patterns
primer_5_prime, primer_3_prime, mismatch_tolerance = read_primer_file(primer_file)
pattern_5_prime = re.compile(f"({primer_5_prime}){{e<={mismatch_tolerance}}}")
pattern_3_prime = re.compile(f"({reverse_complement(primer_3_prime)}){{e<={mismatch_tolerance}}}")

# Output file paths
output_fasta = "amplicon_sequences.fasta"
bin_fasta = "sequences_without_primers.fasta"

# Open output files and process records with threading for speedup
with open(output_fasta, "w") as out_handle, open(bin_fasta, "w") as bin_handle:
    with ThreadPoolExecutor() as executor:
        results = executor.map(
            lambda record: process_sequence(record, pattern_5_prime, pattern_3_prime),
            SeqIO.parse(input_fasta, "fasta")
        )
        for record_id, seq, is_amplicon in results:
            if is_amplicon:
                out_handle.write(f">{record_id}\n{seq}\n")
            else:
                bin_handle.write(f">{record_id}\n{seq}\n")
