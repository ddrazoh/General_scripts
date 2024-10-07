#!/bin/bash

# Define the input and output directories
input_folder="~/output"
output_folder="~/renamed_output"

# Ensure the output directory exists
mkdir -p "$output_folder"

# Iterate over each fasta file in the input folder
for fasta_file in "$input_folder"/*.fasta; do
    # Extract the base name of the file without extension
    base_name=$(basename "$fasta_file" .fasta)
    # Initialize a counter for the sequence number
    seq_count=0
    # Read the fasta file line by line
    while IFS= read -r line; do
        # Check if the line starts with '>'
        if [[ $line == \>* ]]; then
            # Increment the sequence counter
            ((seq_count++))
            # Construct the new sequence name
            new_header=">${base_name}_${seq_count}"
            # Print the new sequence name
            echo "$new_header" >> "$output_folder/$base_name.fasta"
        else
            # Print the line with the new sequence name
            echo "$line" >> "$output_folder/$base_name.fasta"
        fi
    done < "$fasta_file"
done
