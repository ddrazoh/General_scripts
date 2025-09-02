import csv

# Define the input and output file paths
input_file = '/Users/drake/Desktop/tf_trial/analysis/5_prime.csv'
output_file = '5p_transformed_data.csv'

# Read the original CSV file
with open(input_file, 'r') as file:
    reader = csv.reader(file)
    rows = list(reader)

# Extract the headers and data rows
headers = rows[0]
data_rows = rows[1:]

# Create a new list of dictionaries with the desired format
transformed_rows = []
for i in range(0, len(data_rows), 2):
    sanger_name = data_rows[i][0]
    sanger_sequence = data_rows[i][1]
    minion_sequence = data_rows[i+1][1] if i+1 < len(data_rows) else ""

    transformed_row = {
        headers[0]: sanger_name,
        'sanger': sanger_sequence,
        'minion': minion_sequence
    }
    transformed_rows.append(transformed_row)

# Write the transformed data to a new CSV file
fieldnames = [headers[0], 'sanger', 'minion']
with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(transformed_rows)

print("Data transformation complete. The transformed data is saved to:", output_file)
