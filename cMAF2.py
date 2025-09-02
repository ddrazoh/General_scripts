#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
from glob import glob

def cmaf(row):
    try:
        return (row[2:6].sum() - row[2:6].max()) / row[2:6].sum()
    except ZeroDivisionError:
        return np.nan

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    if not os.path.isdir(input_folder):
        print(f"Error: '{input_folder}' is not a valid directory.")
        sys.exit(1)

    # Recursive search for files
    flist = glob(os.path.join(input_folder, '**', '*BaseFreqs.csv'), recursive=True)
    if not flist:
        print(f"No '*BaseFreqs.csv' files found in {input_folder}")
        sys.exit(1)

    dfs = []
    for f in flist:
        try:
            df = pd.read_csv(f)
        except pd.errors.EmptyDataError:
            print(f"Warning: skipping empty file {f}")
            continue
        except Exception as e:
            print(f"Warning: skipping file {f} due to error: {e}")
            continue

        df['sampleid'] = os.path.basename(f).split('_')[0]
        df['cmaf'] = df.apply(cmaf, axis=1)
        dfs.append(df)

    if not dfs:
        print("No valid data to process.")
        sys.exit(1)

    ddf = pd.concat(dfs, ignore_index=True)
    output_file = f"cMAF_{pd.Timestamp.today().date()}.csv"
    ddf.groupby(['sampleid', 'position'])['cmaf'].max().unstack().to_csv(output_file)

    print(f"cMAF dataset saved to {output_file}")

if __name__ == "__main__":
    main()
