import pandas as pd

# Step 1: Read the data from the CSV file
data = pd.read_csv('~/Desktop/tf_trial/analysis/analysis2/5p_transformed_data.csv')

# Step 2: Select the "Sample_ID" and "MinION2" columns
selected_data = data[['Sample_id', 'MinION2']]

# Step 3: Save the selected data as a FASTA file
with open('output.fasta', 'w') as f:
    for _, row in selected_data.iterrows():
        sample_id = row['Sample_id']
        sequence = row['MinION2']
        f.write(f'>{sample_id}\n{sequence}\n')
