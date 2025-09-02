import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from tqdm import tqdm

def kmer_overlap(seq1, seq2, k, min_overlap, max_mismatches):
    kmer_dict = {}
    
    # Build k-mer dictionary for seq1
    for i in range(len(seq1) - k + 1):
        kmer = str(seq1[i:i+k])
        kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    
    # Check for overlaps in seq2
    for i in range(len(seq2) - k + 1):
        kmer = str(seq2[i:i+k])
        if kmer in kmer_dict:
            # Calculate overlap size
            overlap_start = i
            overlap_size = 0
            mismatches = 0
            
            while overlap_start < len(seq2) and overlap_size < min_overlap:
                if seq1[-(overlap_size + 1)] == seq2[overlap_start]:
                    overlap_size += 1
                else:
                    mismatches += 1
                overlap_start += 1
                
                if mismatches > max_mismatches:
                    break
            
            if overlap_size >= min_overlap and mismatches <= max_mismatches:
                return str(seq1)[:-overlap_size] + str(seq2), overlap_size, mismatches
            
    return None, 0, 0

def process_read(read, reads_file2, used_reads1, used_reads2, lock, min_overlap, max_mismatches):
    for read2 in reads_file2:
        with lock:
            if read.id in used_reads1 or read2.id in used_reads2:
                continue
        
        joined_seq, overlap_size, mismatches = kmer_overlap(read.seq, read2.seq, k=20, min_overlap=min_overlap, max_mismatches=max_mismatches)
        
        if joined_seq:
            with lock:
                used_reads1.add(read.id)
                used_reads2.add(read2.id)
            print(f"Joined: {read.id} and {read2.id} | Overlap Size: {overlap_size} | Mismatches: {mismatches}")
            return read.id, read2.id, joined_seq
            
    return None

# Command-line arguments
parser = argparse.ArgumentParser(description="Join paired reads from two FASTA files at overlapping regions.")
parser.add_argument("input_fasta1", help="First input FASTA file")
parser.add_argument("input_fasta2", help="Second input FASTA file")
parser.add_argument("-o", "--output", default="joined_sequences.fasta", help="Output FASTA file")
parser.add_argument("-m", "--min_overlap", type=int, default=45, help="Minimum overlap length")
parser.add_argument("-e", "--mismatches", type=int, default=2, help="Allowed mismatches in overlap")
parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")

args = parser.parse_args()

# Load reads from both files
reads_file1 = list(SeqIO.parse(args.input_fasta1, "fasta"))
reads_file2 = list(SeqIO.parse(args.input_fasta2, "fasta"))

joined_sequences = []
used_reads1 = set()  # Set to keep track of used reads from the first file
used_reads2 = set()  # Set to keep track of used reads from the second file

print(f"Processing {len(reads_file1)} reads from {args.input_fasta1} and {len(reads_file2)} reads from {args.input_fasta2}...")

# Create a lock for managing access to used reads
lock = Lock()

# Process each read in file1 against all reads in file2 using threading
with ThreadPoolExecutor(max_workers=args.threads) as executor:
    futures = [
        executor.submit(process_read, read, reads_file2, used_reads1, used_reads2, lock, args.min_overlap, args.mismatches)
        for read in reads_file1
        if read.id not in used_reads1  # Only process if read1 is not already used
    ]
    
    # Use tqdm to show progress for the futures
    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Joining reads"):
        result = future.result()
        if result:
            read_id1, read_id2, joined_seq = result
            joined_record = SeqIO.SeqRecord(Seq(joined_seq), id=f"{read_id1}_{read_id2}_joined", description="Joined sequence")
            joined_sequences.append(joined_record)

# Write joined sequences to the output file
if joined_sequences:
    with open(args.output, "w") as output_handle:
        SeqIO.write(joined_sequences, output_handle, "fasta")
    print(f"Joined sequences written to {args.output}")
else:
    print("No sequences were joined.")
