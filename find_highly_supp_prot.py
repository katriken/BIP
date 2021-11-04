#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes a .tsv file with peptides from proteins supported by >= 2 
protein-specific peptides from one tissue and a .gtf file with these proteins.
"""

import pandas as pd
import os
import re
import sys


def find_highly_supp_prot(dirpath):
    '''
    Input: directory with mapped peptides (has a .tsv file for each tissue).
    Output: df with protein-specific peptides from proteins with >= 2 
    protein-specific peptides. It has following columns: peptide sequence, 
    protein id, tissues.
    '''
    # Opening mapping results
    print('Opening mapping results...')
    df = pd.DataFrame()
    for filename in os.listdir(dirpath):
        filepath = os.path.join(dirpath, filename)
        df_i = pd.read_csv(filepath, sep='\t', index_col=False, names=[
                           'sequence', 'protein', 'unique'])
        print(filename+':', len(df_i))
        df = pd.concat([df, df_i], ignore_index=True)
    # Selecting only protein-specific peptides
    print('Selecting only protein-specific peptides...')
    print('Total number of peptides:', len(
        set(df['sequence'].values.tolist())))
    df = df[df['unique'] == '+']
    df = df[['sequence', 'protein']]
    df = df.drop_duplicates(subset=['sequence', 'protein'], keep='first')
    pept = df['sequence'].values.tolist()
    print('Number of protein-specific peptides:', len(pept))
    # Finding tissue information for the peptides
    print('Finding tissue information for the peptides...')
    df2 = tissue_info(dirpath, pept)
    df = pd.merge(df, df2, how='inner', on='sequence')
    # Finding proteins supported by >=2 protein-specific peptides from one tissue
    print('Finding proteins supported by >=2 protein-specific peptides from one tissue...')
    prot = list(set(df['protein'].values.tolist()))
    hsp = []
    for p in prot:
        df_i = df[df['protein'] == p]
        tissues = df_i['tissue'].values.tolist()
        tissues = [x.split(',') for x in tissues]
        tissues = [y for x in tissues for y in x]
        if len(set(tissues)) < len(tissues):
            hsp.append(p)
    df = df[df['protein'].isin(hsp)]
    print('Number of protein-specific peptides from proteins with >= 2 such peptides:',
          len(df['sequence'].values.tolist()))
    print('Number of highly supported proteins (>=2 protein-specific peptides):', len(hsp))
    return df


def tissue_info(dirpath, pept):
    '''
    Input: directory with mapped peptides (has a .tsv file for each tissue), 
    list of peptide sequences, which tissue info should be found.
    Output: df with tissue info, has columns: peptide sequence, tissues.
    '''
    # Opening mapping files
    df = pd.DataFrame()
    for filename in os.listdir(dirpath):
        filepath = os.path.join(dirpath, filename)
        df_i = pd.read_csv(filepath, sep='\t', index_col=False, names=[
                           'sequence', 'protein', 'unique'])
        df_i.loc[:, 'tissue'] = filename[:-4]
        df = pd.concat([df, df_i], ignore_index=True)
    df = df[['sequence', 'tissue']]
    df = df.drop_duplicates(subset=['sequence', 'tissue'], keep='first')
    # Selecting only relevant peptides
    df = df[df['sequence'].isin(pept)]
    # Finding list of tissues for each peptide
    tsv = df.values.tolist()
    tsv.sort(key=lambda x: x[0])
    i = 0
    while i < len(tsv)-1:
        if tsv[i][0] == tsv[i+1][0]:
            tsv[i][1] = ",".join([tsv[i][1], tsv[i+1][1]])
            tsv.pop(i+1)
        else:
            i += 1
    df = pd.DataFrame(tsv, columns=['sequence', 'tissue'])
    return df


def open_pred(filename):
    '''
    Input: .gtf file
    Output: df with additional columns for protein ids and gene ids.
    '''
    df = pd.read_csv(filename, sep='\t', index_col=False, names=[
                     'Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl'])
    df = df[df['Type'] == 'CDS']
    suppl = df[['Suppl']].values.tolist()
    prot_ids = []
    gene_ids = []
    for i in range(len(suppl)):
        text = suppl[i][0]
        prot_id = re.findall("transcript_id \"\S+", text)[0]
        prot_id = re.findall("\"\S+", prot_id)[0][1:-2]
        prot_ids.append(prot_id)
        gene_id = re.findall("gene_id \"\S+", text)[0]
        gene_id = re.findall("\"\S+", gene_id)[0][1:-2]
        gene_ids.append(gene_id)
    df['protein'] = prot_ids
    df['gene'] = gene_ids
    return df


def make_tsv(df, filename):
    '''
    Input: df
    Output: .tsv file w/o header
    '''
    tsv = df.values.tolist()
    with open(filename, 'w') as f:
        for line in tsv:
            for word in line:
                f.write(str(word)+'\t')
            f.write('\n')


def save_new_gtf(proteins, file_pred, file_output):
    '''
    Input: 
    - proteins: protein ids,
    - file_pred: .gtf file with prediction,
    - file_output: output .gtf file
    Output: .gtf file with selected proteins
    '''
    df = open_pred(file_pred)
    df = df[df['protein'].isin(proteins)]
    df = df[['Seqid', 'Source', 'Type', 'Start',
             'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    make_tsv(df, file_output)


def main(argv):
    file_pred = argv[0]
    dir_mapped = argv[1]
    file_output_gtf = argv[2]

    df = find_highly_supp_prot(dir_mapped)
    make_tsv(df, file_output_gtf[:-4]+'.tsv')
    hsp = list(set(df['protein'].values.tolist()))
    save_new_gtf(hsp, file_pred, file_output_gtf)


if __name__ == "__main__":
    main(sys.argv[1:])
