#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes a .tsv file with gene-specific peptides from genes with >=2 such peptides
"""

import pandas as pd
import os
import re
import sys


def open_pred(filename):
    '''
    Input: .gtf file
    Output: df with additional columns for protein ids and gene ids
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


def find_gene_spec_pept(dirpath, df_prot_genes):
    '''
    Input: directory with mapped peptides (has a .tsv file for each tissue),
    df with protein/gene ids.
    Output: df with gene-specific peptides from genes with >= 2 gene-specific 
    peptides. It has following columns: peptide sequence, gene id, unique ('+' if 
    protein-specific), protein id.
    '''
    # Opening mapping results
    print('Opening mapping results...')
    df_mapping = pd.DataFrame()
    for filename in os.listdir(dirpath):
        filepath = os.path.join(dirpath, filename)
        df_i = pd.read_csv(filepath, sep='\t', index_col=False, names=[
                           'sequence', 'protein', 'unique'])
        print(filename+':', len(df_i))
        df_mapping = pd.concat([df_mapping, df_i], ignore_index=True)
    proteins = df_mapping['protein'].values.tolist()
    # Adding gene information
    df_mapping = pd.merge(df_mapping, df_prot_genes, how='inner', on='protein')
    # Finding all genes with at least one gene-specific peptide
    print('Finding all genes with at least one gene-specific peptide...')
    df = df_mapping.copy()
    df = df[['sequence', 'gene', 'unique']]
    df = df.drop_duplicates(subset=['sequence', 'gene'], keep='first')
    pept = list(set(df['sequence'].values.tolist()))
    print('Total number of all peptides (also unspecific):', len(pept))
    gene_unique = []
    for p in pept:
        df_i = df[df['sequence'] == p]
        if len(df_i) == 1:
            gene_unique.append(p)
    df = df[df['sequence'].isin(gene_unique)]
    print('All gene-specific peptides:', len(gene_unique))
    print('Genes with $\geq$1 gene-specific peptide:',
          len(set(df['gene'].values.tolist())))
    # Finding genes with >= 2 such peptides
    print('Finding genes with >= 2 such peptides...')
    genes = df['gene'].values.tolist()
    cand = []
    genes_set = list(set(genes))
    for gene in genes_set:
        if genes.count(gene) > 1:
            cand.append(gene)
    df = df[df['gene'].isin(cand)]
    print('Genes with $\geq$2 gene-specific peptides:', len(cand))
    print('Gene-specific peptides from genes with >= 2 gene-specific peptides',
          len(set(df['sequence'].values.tolist())))
    # Adding protein information
    df_mapping_info = df_mapping[['sequence', 'protein']]
    df_mapping_info = df_mapping_info.drop_duplicates(
        subset=['sequence', 'protein'], keep='first')
    df = pd.merge(df, df_mapping_info, how='inner', on=['sequence'])
    # Adding tissue info to the peptides
    print('Adding tissue info to the peptides...')
    pept = list(set(df['sequence'].values.tolist()))
    df_t = tissue_info(dirpath, pept)
    df = pd.merge(df, df_t, how='inner', on=['sequence'])
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


def main(argv):
    file_pred = argv[0]
    dir_mapped = argv[1]
    file_output_tsv = argv[2]

    df_pred = open_pred(file_pred)
    df_prot_genes = df_pred[['protein', 'gene']]
    df_prot_genes = df_prot_genes.drop_duplicates(
        subset=['protein'], keep='first')
    df = find_gene_spec_pept(dir_mapped, df_prot_genes)
    make_tsv(df, file_output_tsv)


if __name__ == "__main__":
    main(sys.argv[1:])
