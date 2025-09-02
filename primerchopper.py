import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re

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

# Get input arguments
if len(sys.argv) != 3:
    print("Usage: python3 extract_amplicons.py <input_fasta_file> <primers_file>")
    sys.exit(1)

input_fasta = sys.argv[1]
primer_file = sys.argv[2]

# Read primers from the specified .txt file
primer_5_prime, primer_3_prime, mismatch_tolerance = read_primer_file(primer_file)

# Prepare regex patterns with allowed mismatches for each primer
pattern_5_prime = f"({primer_5_prime}){{e<={mismatch_tolerance}}}"
pattern_3_prime = f"({reverse_complement(primer_3_prime)}){{e<={mismatch_tolerance}}}"

# Output file paths
output_fasta = "amplicon_sequences.fasta"
bin_fasta = "sequences_without_primers.fasta"  # File for sequences missing primers

# Debugging: Print primer patterns and mismatch tolerance
print("5' Primer Pattern:", pattern_5_prime)
print("3' Primer Pattern:", pattern_3_prime)
print("Mismatch Tolerance:", mismatch_tolerance)

# Open output files
with open(output_fasta, "w") as out_handle, open(bin_fasta, "w") as bin_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq)  # Keep the original sequence
        print(f"Processing sequence ID: {record.id}, Length: {len(sequence)}")
        
        # Search for the 5' primer (start of amplicon)
        match_5_prime = re.search(pattern_5_prime, sequence)
        print(f"5' Primer Match: {match_5_prime}")

        # Search for the 3' primer (end of amplicon) using reverse complement
        match_3_prime = re.search(pattern_3_prime, sequence)
        print(f"3' Primer Match: {match_3_prime}")

        if match_5_prime and match_3_prime:
            amplicon_start = match_5_prime.end()  # End of 5' primer
            amplicon_end = match_3_prime.start()  # Start of 3' primer
            
            # Debugging output for amplicon boundaries
            print(f"Amplicon Start: {amplicon_start}, Amplicon End: {amplicon_end}")

            # Ensure amplicon boundaries are valid
            if amplicon_start < amplicon_end:
                amplicon_seq = sequence[amplicon_start:amplicon_end]
                
                # Check if the amplicon sequence is valid (not empty)
                if amplicon_seq:  # This should hold true based on your debug output
                    record.seq = Seq(amplicon_seq)  # Update the record's sequence
                    print(f"Writing record ID: {record.id} with sequence: {record.seq}")  # Debugging line
                    SeqIO.write(record, out_handle, "fasta")
                else:
                    print(f"No valid amplicon found in {record.id}. Amplicon sequence is empty.")
            else:
                print(f"Invalid amplicon boundaries for {record.id}: start {amplicon_start} is not less than end {amplicon_end}.")
        else:
            # Write the original sequence to the bin file if either primer is missing
            if sequence:  # Ensure the sequence is not empty
                print(f"Writing original sequence to bin file for {record.id}.")
                record.seq = Seq(sequence)  # Keep the original sequence
                SeqIO.write(record, bin_handle, "fasta")
            else:
                print(f"Skipping {record.id}: invalid sequence.")
