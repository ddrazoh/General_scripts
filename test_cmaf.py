import numpy as np
import pandas as pd
from glob import glob

flist = glob('*/*BaseFreqs.csv')

def cmaf(x):
    try:
        return (x[2:6].sum() - x[2:6].max()) / x[2:6].sum()
    except ZeroDivisionError:
        return np.nan

dfs = []
for f in flist:
    df = pd.read_csv(f)
    df.columns = ['position', 'recall', 'a', 'c', 'g', 't', 'gap', 'n']
    df['sampleid'] = f.split('/')[-1].split('_')[0]
    df['cmaf'] = df.apply(cmaf, axis=1)
    dfs.append(df)

ddf = pd.concat(dfs)

grouped = ddf.groupby(['sampleid', 'position']).cmaf.max().unstack()
grouped.to_csv('cMAF_2024-10-03.csv')





##use HXB2 CODINATES

import numpy as np
import pandas as pd
from glob import glob

flist = glob('*/*HXB2.csv')

def cmaf(x):
    try:
        return (x[2:6].sum() - x[2:6].max()) / x[2:6].sum()
    except ZeroDivisionError:
        return np.nan

dfs = []
for f in flist:
    df = pd.read_csv(f)
    df = df.iloc[:, [0, 2, 3, 4, 5, 6, 7, 8]]

    df.columns = ['position', 'recall', 'a', 'c', 'g', 't', 'gap', 'n']
    df['sampleid'] = f.split('/')[-1].split('_')[0]
    df['cmaf'] = df.apply(cmaf, axis=1)
    dfs.append(df)
df 
ddf = pd.concat(dfs)

ddf = ddf[ddf['position'] != '-']

# Convert position to numeric so it sorts as numbers, not strings
ddf['position'] = pd.to_numeric(ddf['position'])


grouped = ddf.groupby(['sampleid', 'position']).cmaf.max().unstack()
grouped.to_csv('cMAF_2027-10-03.csv')
