import sys

def convert_to_fasta(input_path, output_path):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        outfile.write(">contig_1\n")  # Write the single header
        for line in infile:
            if not line.startswith(">"):  # Skip headers in the input file
                outfile.write(line.strip())  # Write sequence as a continuous line

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_fasta.py <input.txt> <output.fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_to_fasta(input_file, output_file)

