import pandas as pd
import numpy as np

# Load your patstats file (comma-delimited, not tab-delimited)
df = pd.read_csv('/Users/drake/Documents/USYD_HIVPHYLO/2025_Mathew/Phyloscanner2.0/to_nadua/577/E_577_patStats.csv', sep=",")

df.columns = df.columns.str.strip()
print("Columns found:", df.columns.tolist())

expanded_rows = []

for _, row in df.iterrows():
    start, end = map(int, row["tree.id"].split("_to_"))
    
    # Generate 10 evenly spaced xcoords
    new_xcoords = np.linspace(start, end, num=100)
    
    for x in new_xcoords:
        new_row = row.copy()
        new_row["xcoord"] = x
        expanded_rows.append(new_row)

expanded_df = pd.DataFrame(expanded_rows)
expanded_df.to_csv("patstats_expanded.csv", index=False)
print("Expanded file written: patstats_expanded2.csv")
