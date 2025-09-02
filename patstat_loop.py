import pandas as pd
import numpy as np
import glob
import os

# Folder containing your CSV files
input_folder = '/Users/drake/Documents/USYD_HIVPHYLO/2025_Mathew/Phyloscanner2.0/to_nadua/patstats_0'
output_folder = '/Users/drake/Documents/USYD_HIVPHYLO/2025_Mathew/Phyloscanner2.0/to_nadua/577/patstats_0/expanded'

os.makedirs(output_folder, exist_ok=True)


for csv_file in glob.glob(os.path.join(input_folder, "*.csv")):
    print(f"Processing {csv_file}...")
    
    df = pd.read_csv(csv_file, sep=",")
    df.columns = df.columns.str.strip()
    
    expanded_rows = []
    
    for _, row in df.iterrows():
        start, end = map(int, row["tree.id"].split("_to_"))
        new_xcoords = np.linspace(start, end, num=100)
        
        for x in new_xcoords:
            new_row = row.copy()
            new_row["xcoord"] = x
            expanded_rows.append(new_row)
    
    expanded_df = pd.DataFrame(expanded_rows)
    
   
    base_name = os.path.basename(csv_file)
    output_file = os.path.join(output_folder, base_name.replace(".csv", "_expanded.csv"))
    expanded_df.to_csv(output_file, index=False)
    
    print(f"Expanded file written: {output_file}")
