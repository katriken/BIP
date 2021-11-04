#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find supported proteins based on BRAKER scores, gene-specific and 
protein-specific peptides. 
"""

import pandas as pd
import re
import sys


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


def find_highly_supported_proteins(file_gene_spec_pept):
    '''
    Input: 
    - file_gene_spec_pept: file with gene-specific peptides from genes with 
        >= 2 gene-specific peptides
    - t_info: dictionary with peptide sequences and tissues where the peptides
         were identified. 
    Output: df with unique peptides from proteins with >= 2 unique peptides 
    from 1 tissue.
    '''
    df1 = pd.read_csv(file_gene_spec_pept, sep='\t', index_col=False, names=[
                      'sequence', 'gene', 'unique', 'protein', 'tissue'])
    df2 = df1[df1['unique'] == '+']
    df2 = df2[['sequence', 'protein', 'tissue']]
    prot2 = df2['protein'].values.tolist()
    prot2_set = list(set(prot2))
    prot3 = []
    for p in prot2_set:
        if prot2.count(p) > 1:
            df2_i = df2[df2['protein'] == p]
            tissues = df2_i['tissue'].values.tolist()
            tissues = [x.split(',') for x in tissues]
            tissues = [y for x in tissues for y in x]
            if len(set(tissues)) < len(tissues):
                prot3.append(p)
    print('Number of proteins supported by >= 2 protein-specific peptides from 1 tissue:', len(prot3))
    df3 = df2[df2['protein'].isin(prot3)]
    return df3


def find_all_supp_proteins(file_gene_spec_pept, highly_supp_prot, t_info, tx_scores):
    '''
    Input: 
    - file_gene_spec_pept: file with gene-specific peptides from genes with 
        >= 2 gene-specific peptides,
    - highly_supp_prot: proteins with >= 2 unique peptides from 1 tissue,
    - t_info: dictionary with peptide sequences and tissues where the peptides
         were identified,
    - tx_scores: dictionary with BRAKER transcript scores, where keys are
          protein ids.
    Output: df with all supported proteins, has columns: peptide sequence, 
    gene id, unique (='+' for protein-specific peptides), protein id.
    '''
    df = pd.read_csv(file_gene_spec_pept, sep='\t', index_col=False, names=[
                     'sequence', 'gene', 'unique', 'protein', 'tissue'])
    genes = list(set(df['gene'].values.tolist()))
    supp_proteins = []
    supp_prot_pept = []
    for g in genes:
        df_i = df[df['gene'] == g]
        prot = list(set(df_i['protein'].values.tolist()))
        seq = [df_i[df_i['protein'] == x]['sequence'].values.tolist()
               for x in prot]
        prot, seq = leave_prot_2pept1tissue(prot, seq, t_info)
        prot, seq = find_supp_proteins(
            prot, seq, tx_scores, t_info, highly_supp_prot)
        if len(prot) != 0:
            supp_proteins.append(prot)
            supp_prot_pept.append(seq)
    supp_proteins = [y for x in supp_proteins for y in x]
    supp_prot_pept = [y for x in supp_prot_pept for y in x]
    df2 = pd.DataFrame()
    for i in range(len(supp_proteins)):
        df_i = df[df['protein'] == supp_proteins[i]]
        df_i = df_i[df_i['sequence'].isin(supp_prot_pept[i])]
        if len(df_i) == 0:
            print('ERROR')
        df2 = pd.concat([df2, df_i], ignore_index=True)
    return df2


def find_supp_proteins(proteins, peptides, tx_scores, t_info, confirmed_proteins):
    '''
    Input: 
    - proteins: list of protein ids of a gene that have >= 2 gene-specific
          peptides from 1 tissue,
    - peptides: 2D list with  gene-specific peptides found in corresponding 
        proteins,
    - tx_scores: dictionary with BRAKER transcript scores, where keys are
          protein ids,
    - t_info: dictionary with peptide sequences and tissues where the peptides
          were identified. All tissues for a peptide are joined in a
          string with ',' as a separator.
    - confirmed_proteins: list of all proteins that have >= 2 protein-specific
          peptides from 1 tissue       
    Output: list of supported proteins. Proteins are selected based on BRAKER 
    scores and gene-specific peptides.
    '''
    # first highly supported proteins selected (>2 protein-specific peptides
    # from 1 tissue)
    res = [x for x in proteins if x in confirmed_proteins]
    res_pept = []
    if len(res) != 0:
        pept_from_res = [peptides[proteins.index(x)] for x in res]
        res_pept += pept_from_res
        pept_from_res = [y for x in pept_from_res for y in x]
        # remove peptides present in 'res' proteins from the pool
        for p in res:
            del peptides[proteins.index(p)]
            del proteins[proteins.index(p)]
        peptides = [list(set(x)-set(pept_from_res)) for x in peptides]
    # select supported proteins, based on BRAKER score -> remove peptides ->
    # select again
    while len(proteins) > 0:
        best = find_best_protein(proteins, peptides, tx_scores, t_info)
        if best == False:
            break
        else:
            res.append(best)
            pept_from_best = peptides[proteins.index(best)]
            res_pept.append(pept_from_best)
            del peptides[proteins.index(best)]
            del proteins[proteins.index(best)]
            peptides = [list(set(x)-set(pept_from_best)) for x in peptides]
    return res, res_pept


def find_best_protein(proteins, peptides, tx_scores, t_info):
    '''
    Input:
    - proteins: list of protein ids of a gene,
    - peptides: 2D list with gene-specific peptides found in corresponding 
        proteins,
    - tx_scores: dictionary with BRAKER transcript scores, where keys are
          protein ids,
    - t_info: dictionary with peptide sequences and tissues where the peptides
          were identified. 
    Output: protein with the highest BRAKER score and >= 2 gene-specific 
    peptides from 1 tissue. If none of the proteins has >= 2 gene-specific 
    peptides from 1 tissue, returns False.
    '''
    proteins, peptides = leave_prot_2pept1tissue(proteins, peptides, t_info)
    if len(proteins) == 0:
        return False
    num_hits = [len(x) for x in peptides]
    proteins = [proteins[i] for i in range(len(proteins)) if num_hits[i] >= 2]
    peptides = [peptides[i] for i in range(len(peptides)) if num_hits[i] >= 2]
    if len(proteins) == 0:
        return False
    else:
        s = [tx_scores[x] for x in proteins]
        return proteins[s.index(max(s))]


def leave_prot_2pept1tissue(proteins, peptides, t_info):
    '''
    Input: 
    - proteins: list of protein ids of a gene,
    - peptides: 2D list with gene-specific peptides found in corresponding 
        proteins,
    - t_info: dictionary with peptide sequences and tissues where the peptides
          were identified. 
    Output: new lists of proteins and peptides for proteins with >= 2 
    gene-specific peptides from 1 tissue.
    '''
    new_prot = []
    new_pept = []
    for i in range(len(proteins)):
        pept = peptides[i]
        t = [t_info[x] for x in pept]
        t = [x.split(',') for x in t]
        t = [y for x in t for y in x]
        if len(set(t)) < len(t):
            new_prot.append(proteins[i])
            new_pept.append(peptides[i])
    return new_prot, new_pept


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
    file_prediction = argv[0]
    file_tx_scores = argv[1]
    file_gene_spec_pept = argv[2]
    file_output_gtf = argv[3]

    # prediction
    df = open_pred(file_prediction)
    prot_all = list(set(df['protein'].values.tolist()))
    print('Number of all proteins:', len(prot_all))
    # dict with sequence & tissue info for gene-specific peptides
    df1 = pd.read_csv(file_gene_spec_pept, sep='\t', index_col=False, names=[
                      'sequence', 'gene', 'unique', 'protein', 'tissue'])
    df1 = df1[['sequence', 'tissue']]
    df1 = df1.drop_duplicates(subset=['sequence'], keep='first')
    seq = df1['sequence'].values.tolist()
    tissue = df1['tissue'].values.tolist()
    t_info = dict(zip(seq, tissue))
    # tx scores
    df_scores = pd.read_csv(file_tx_scores, sep='\t',
                            index_col=False, names=['protein', 'score'])
    df_scores.loc[(df_scores.score == '.'), 'score'] = '0.001'
    df_scores['score'] = df_scores['score'].astype(str).astype(float)
    prot_scores = df_scores['score'].values.tolist()
    prot_names = df_scores['protein'].values.tolist()
    tx_scores = dict(zip(prot_names, prot_scores))
    # Finding highly supported proteins (>= 2 protein-specific peptides from 1 tissue)
    print('Finding highly supported proteins (>= 2 protein-specific peptides from 1 tissue)...')
    df2 = find_highly_supported_proteins(file_gene_spec_pept)
    prot2 = list(set(df2['protein'].values.tolist()))
    # Selecting supported proteins
    print('Selecting supported proteins...')
    df3 = find_all_supp_proteins(file_gene_spec_pept, prot2, t_info, tx_scores)
    make_tsv(df3, file_output_gtf[:-4]+'.tsv')
    prot3 = list(set(df3['protein'].values.tolist()))
    print('Number of supported proteins:', len(prot3))
    save_new_gtf(prot3, file_prediction, file_output_gtf)


if __name__ == "__main__":
    main(sys.argv[1:])
