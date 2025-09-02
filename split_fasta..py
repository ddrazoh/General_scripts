def create_fasta_files(input_file):
    with open(input_file, 'r') as f:
        content = f.read().split('>')[1:]

    for entry in content:
        lines = entry.split('\n')
        sequence_name = lines[0].split()[0]
        sequence = ''.join(lines[1:])
        output_filename = f'{sequence_name}.fasta'

        with open(output_filename, 'w') as output_file:
            output_file.write(f'>{sequence_name}\n{sequence}\n')


# Usage
input_fasta = '/Users/drake/Desktop/Alfred/sorted_kleb.fasta'
create_fasta_files(input_fasta)
