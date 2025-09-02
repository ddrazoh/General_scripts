#!/usr/bin/env python
import argparse
import pandas as pd

def load_reference_data(modeldir):
    hxb2 = pd.read_csv(f'{modeldir}/HXB2_refdata.csv')
    hxb2['position'] = hxb2['HXB2 base position']
    hxb2.set_index('position', inplace=True)
    rf3_3cp = hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].max()
    rf2_3cp = hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].max()
    rf1_3cp = hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].max()
    t1 = set(rf1_3cp.reset_index()['HXB2 base position'])
    t2 = set(rf2_3cp.reset_index()['HXB2 base position'])
    t3 = set(rf3_3cp.reset_index()['HXB2 base position'])
    rf1_12 = set(hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf2_12 = set(hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf3_12 = set(hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    first_second_codon_pos = set(rf1_12 | rf2_12 | rf3_12)
    third_codon_pos = set(t1 | t2 | t3) - first_second_codon_pos
    gag = set(hxb2[hxb2['RF1 protein']=='gag'].index)
    pol = set(hxb2[hxb2['RF3 protein']=='pol'].index)
    gp120 = set(hxb2[hxb2['RF3 protein']=='gp120'].index)
    gp41 = set(hxb2[hxb2['RF3 protein']=='gp41'].index)
    return first_second_codon_pos, third_codon_pos, gag, pol, gp120, gp41

def load_maf(fpath):
    df = pd.read_csv(fpath, index_col=0)
    df.columns = df.columns.astype(int)
    return df

def compute_maf_aggregates(Xmaf, modeldir):
    first_second_codon_pos, third_codon_pos, gag, pol, gp120, gp41 = load_reference_data(modeldir)
    results = pd.DataFrame(index=Xmaf.index)
    results['genome_maf12c'] = Xmaf.loc[:, list(first_second_codon_pos & set(Xmaf.columns))].mean(axis=1)
    results['genome_maf3c'] = Xmaf.loc[:, list(third_codon_pos & set(Xmaf.columns))].mean(axis=1)
    for gene, gene_set in zip(['gag', 'pol', 'gp120', 'gp41'], [gag, pol, gp120, gp41]):
        pos12 = first_second_codon_pos & gene_set & set(Xmaf.columns)
        pos3 = third_codon_pos & gene_set & set(Xmaf.columns)
        results[f'{gene}_maf12c'] = Xmaf.loc[:, list(pos12)].mean(axis=1)
        results[f'{gene}_maf3c'] = Xmaf.loc[:, list(pos3)].mean(axis=1)
    results.reset_index(inplace=True)
    results.rename(columns={'index': 'sampleid'}, inplace=True)
    return results

def main():
    import sys
    parser = argparse.ArgumentParser(description="Aggregate MAFs by codon position and gene")
    parser.add_argument('-m', '--mafpath', required=True, help='Input CMAF CSV file')
    parser.add_argument('-d', '--modeldir', required=True, help='Directory with HXB2_refdata.csv')
    parser.add_argument('-o', '--outpath', required=True, help='Output CSV file path')
    args = parser.parse_args()
    
    Xmaf = load_maf(args.mafpath)
    maf_summary = compute_maf_aggregates(Xmaf, args.modeldir)
    maf_summary.to_csv(args.outpath, index=False)
    print(f"MAF summary saved to {args.outpath}")

if __name__ == '__main__':
    main()
