
###WORKS WELL
for d in */; do
  (
    cd "$d" || exit
    #samtools bam2fq "$(basename "$d").bam" > "$(basename "$d")_reads.fastq"

    minimap2 -t 7 -ax map-ont --secondary=no *ref.fasta *reads* > "$(basename "$d").sam"

    samtools view -b -q 7 -F 4 "$(basename "$d").sam" > "$(basename "$d")_.bam"
    samtools sort "$(basename "$d")_.bam" -o "$(basename "$d").sorted.bam"
    samtools index "$(basename "$d").sorted.bam"

    samtools faidx *ref.fasta

    samtools mpileup --no-BAQ --min-BQ 5 --max-depth 1000000 --reference *ref.fasta "$(basename "$d").sorted.bam" > "$(basename "$d").pileup"

    python3 /Users/drake/opt/anaconda3/envs/minion/bin/tools/AnalysePileup.py "$(basename "$d").pileup" *ref.fasta > "$(basename "$d")_BaseFreqs.csv"

    cat *ref.fasta /Users/drake/opt/anaconda3/pkgs/shiver-1.7.3-hdfd78af_0/data/external/B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta | sed -e 's/?/N/g' -e 's/\(.\)\(>\)/\1\n\2/g' > "HXB2_$(basename "$d").fasta"


    mafft --adjustdirection "HXB2_$(basename "$d").fasta" > "msa_HXB2_$(basename "$d").fasta"

    /Users/drake/opt/anaconda3/envs/minion/bin/tools/MergeBaseFreqsAndCoords.py --pairwise-aln "msa_HXB2_$(basename "$d").fasta" "$(basename "$d")_BaseFreqs.csv" > "$(basename "$d")_BaseFreqs_HXB2.csv"
  )
done



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
grouped.to_csv('cMAF_2025-08-10.csv')




##testmafs.py

/Users/drake/Downloads/HIV-phyloTSI-main/HIVPhyloTSI_maftest.py -m /Users/drake/Documents/USYD_HIVPHYLO/cMAF_2024-11-09.csv -d /Users/drake/Downloads/HIV-phyloTSI-main/Model -o testmaf_out.csv





import numpy as np
import pandas as pd
from glob import glob

# Get a list of all CSV files matching the pattern */*HXB2.csv
flist = glob('*/*HXB2.csv')

# Function to calculate cMAF (complementary minor allele frequency)
def cmaf(x):
    try:
        # Take sum of counts for A, C, G, T (columns 2–5 in zero-based indexing)
        # Subtract the maximum count (major allele)
        # Divide by the total counts for A, C, G, T to get minor allele frequency
        return (x[2:6].sum() - x[2:6].max()) / x[2:6].sum()
    except ZeroDivisionError:
        # If all counts are zero, return NaN to avoid division error
        return np.nan

dfs = []
for f in flist:
    # Read each CSV file
    df = pd.read_csv(f)
    # Keep only position, recall, A, C, G, T, gap, and N columns
    df = df.iloc[:, [0, 2, 3, 4, 5, 6, 7, 8]]

    # Rename columns for clarity
    df.columns = ['position', 'recall', 'a', 'c', 'g', 't', 'gap', 'n']
    # Extract sample ID from filename (before the first underscore)
    df['sampleid'] = f.split('/')[-1].split('_')[0]
    # Apply cmaf() to each row
    df['cmaf'] = df.apply(cmaf, axis=1)
    # Store processed DataFrame
    dfs.append(df)

# Concatenate all DataFrames into one
ddf = pd.concat(dfs)

# Remove rows where position is missing ("-")
ddf = ddf[ddf['position'] != '-']

# Convert position column to numeric type for correct sorting
ddf['position'] = pd.to_numeric(ddf['position'])

# Group by sample and position, take maximum cMAF per position, and reshape table
grouped = ddf.groupby(['sampleid', 'position']).cmaf.max().unstack()

# Save grouped data to CSV
grouped.to_csv('cMAF_2025-08-10.csv')

#concactenate_patstats_keeping one header

awk 'FNR==1 && NR!=1 { next } { print }' *.csv > merged.csv

##run_hiv_phylotsi

/Users/drake/Downloads/HIV-phyloTSI-main/HIVPhyloTSI.py -d /Users/drake/Downloads/HIV-phyloTSI-main/Model  -p nadua1_patstats.csv -m ../dup2_nadua_cMAF_2025-08-10.csv -o ../nadua_duration_estimates.csv


###run phyloscanner 
phyloscanner_make_trees InputFileList.csv --windows 475,5068,4956,9636 --no-trees -P -A /mnt/myvolume/drake/recency/hiv-db_refs.fasta -2 B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XR B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XC '823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886' --min-read-count 1 --merging-threshold-a 0

