#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unites two gene predictions
- It should be specified which one is the main prediction.
  All entries from the main prediction will have a corresponding label, but only
  new proteins from the additional prediction that are not present in the main
  prediction will have a label of the additional prediction.
"""

import pandas as pd
import re
import os
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
    df = df.assign(protein=prot_ids)
    df = df.assign(gene=gene_ids)
    df['Frame'] = df['Frame'].astype(str).astype(int)
    return df


def len_overlap(cds1, cds2):
    '''
    Input: 2 lists with cds positions. Each list has CDS from all proteins 
    of a gene.
    Output: CDS overlap between two genes in nucleotides.
    '''
    # overlaping fragments
    overlapping_pos = []
    for a in cds1:
        for b in cds2:
            if (a[0] >= b[0]) and (a[1] <= b[1]):
                overlapping_pos.append(a)
            elif (a[0] <= b[0]) and (a[1] >= b[1]):
                overlapping_pos.append(b)
            elif (a[1] >= b[0]) and (a[1] <= b[1]):
                overlapping_pos.append([b[0], a[1]])
            elif (a[0] >= b[0]) and (a[0] <= b[1]):
                overlapping_pos.append([a[0], b[1]])

    # combine overlapping fragments, if they overlap with each other
    overlapping_pos.sort()
    i = 0
    while i < len(overlapping_pos)-1:
        if overlapping_pos[i][1] > overlapping_pos[i+1][0]:
            overlapping_pos[i] = [overlapping_pos[i]
                                  [0], overlapping_pos[i+1][1]]
            overlapping_pos.pop(i+1)
        else:
            i += 1
    overlapping_pos_nt = sum([x[1]-x[0]+1 for x in overlapping_pos])
    return overlapping_pos_nt


def unite_pred(df1, df2, min_overlap, pred_name1, pred_name2):
    '''
    Input: main (df1) and additional (df2) predictions, minimal CDS
    overlap (nt) between proteins of a gene, labels for both predictions.
    Output: df with united prediction
    '''
    # change protein/gene names for one of the predictions, so that there are no common names
    prot2 = df2['protein'].values.tolist()
    prot2 = ['x'+x for x in prot2]
    genes2 = df2['gene'].values.tolist()
    genes2 = ['x'+x for x in genes2]
    df2['protein'] = prot2
    df2['gene'] = genes2

    # find overlapping genes
    new_genes2 = []
    overlapping_genes2 = []
    overlapping_genes1 = []
    genes2_set = list(set(genes2))
    for g in genes2_set:
        df2_i = df2[df2['gene'] == g]
        strand = df2_i.iat[0, 6]
        seqid = df2_i.iat[0, 0]
        start = min(df2_i['Start'].values.tolist())
        end = max(df2_i['End'].values.tolist())
        df2_i = df2_i[['Start', 'End', 'Frame']]
        df2_i = df2_i.drop_duplicates(subset=['Start', 'End', 'Frame'])
        df1_i = df1[df1['Seqid'] == seqid]
        df1_i = df1_i[df1_i['Strand'] == strand]
        df1_i = df1_i[df1_i['Start'] < end]
        df1_i = df1_i[df1_i['End'] > start]
        df1_i = df1_i[['Start', 'End', 'Frame', 'gene']]
        df1_i = df1_i.drop_duplicates(subset=['Start', 'End', 'Frame', 'gene'])
        cand_genes = list(set(df1_i['gene'].values.tolist()))
        df_common = pd.merge(df1_i, df2_i, how='inner',
                             on=['Start', 'End', 'Frame'])
        # no overlaps with other genes
        if len(df1_i) == 0:
            new_genes2.append(g)
        # one possible overlapping gene
        elif len(cand_genes) == 1:
            # all CDS from a new gene (df2_i) are in the main prediction(df1_i). Very common
            if len(df_common) == len(df2_i):
                overlapping_genes1.append(cand_genes[0])
                overlapping_genes2.append(g)
            # there is an overlap -> check overlap length
            else:
                cds1 = df1_i[['Start', 'End']].values.tolist()
                cds2 = df2_i[['Start', 'End']].values.tolist()
                if len_overlap(cds1, cds2) >= min_overlap:
                    overlapping_genes1.append(cand_genes[0])
                    overlapping_genes2.append(g)
                else:
                    new_genes2.append(g)
        # >one possible overlapping genes
        else:
            best_gene = 'x'
            best_overlap = 0
            for cand in cand_genes:
                cds1 = df1_i[df1_i['gene'] == cand][[
                    'Start', 'End']].values.tolist()
                cds2 = df2_i[['Start', 'End']].values.tolist()
                overlap = len_overlap(cds1, cds2)
                if (overlap >= min_overlap) and (overlap > best_overlap):
                    best_gene = cand
                    best_overlap = overlap
            if best_gene == 'x':
                new_genes2.append(g)
            else:
                overlapping_genes1.append(best_gene)
                overlapping_genes2.append(g)
    print('New genes:', len(new_genes2),
          'Overlapping genes:', len(overlapping_genes2))

    # change gene names of overlapping genes in df2 to corresponding
    # names from df1
    df2_overlap = df2[df2['gene'].isin(overlapping_genes2)]
    names_converter = dict(zip(overlapping_genes2, overlapping_genes1))
    genes2_overlap = df2_overlap['gene'].values.tolist()
    df2_overlap = df2_overlap.assign(
        gene=[names_converter[x] for x in genes2_overlap])
    df2_rest = df2[df2['gene'].isin(new_genes2)]
    df2_new = pd.concat([df2_overlap, df2_rest], ignore_index=True)
    print('Gene names changed.')

    # rm df2 proteins present in df1
    chrs = list(set(df2['Seqid'].values.tolist()))
    prot1_info = []
    for chrm in chrs:
        df1_i = df1[df1['Seqid'] == chrm]
        prot1_set = list(set(df1_i['protein'].values.tolist()))
        for p in prot1_set:
            df1_ii = df1_i[df1_i['protein'] == p]
            cds = df1_ii[['Start', 'End']].values.tolist()
            cds.sort()
            cds = [y for x in cds for y in x]
            prot1_info.append(chrm+'.'+'.'.join(map(str, cds)))
    df1_info = pd.DataFrame()
    df1_info = df1_info.assign(info=prot1_info)
    prot2_info = []
    prot2_names = []
    for chrm in chrs:
        df2_i = df2[df2['Seqid'] == chrm]
        prot2_set = list(set(df2_i['protein'].values.tolist()))
        for p in prot2_set:
            df2_ii = df2_i[df2_i['protein'] == p]
            cds = df2_ii[['Start', 'End']].values.tolist()
            cds.sort()
            cds = [y for x in cds for y in x]
            prot2_info.append(chrm+'.'+'.'.join(map(str, cds)))
            prot2_names.append(p)
    df2_info = pd.DataFrame()
    df2_info = df2_info.assign(protein=prot2_names)
    df2_info = df2_info.assign(info=prot2_info)
    df2_info = pd.merge(df2_info, df1_info, how='inner', on=['info'])
    prot2_common = list(set(df2_info['protein'].values.tolist()))
    prot2_all = list(set(df2['protein'].values.tolist()))
    prot2_new = list(set(prot2_all) - set(prot2_common))
    df2_new = df2_new[df2_new['protein'].isin(prot2_new)]

    # combine df1 and df2_new
    df1['prediction'] = pred_name1
    df2_new['prediction'] = pred_name2
    df = pd.concat([df1, df2_new], ignore_index=True)

    # Correct Suppl column
    df['Suppl'] = 'transcript_id "'+df['protein']+'"; gene_id "' + \
        df['gene']+'"; prediction "'+df['prediction']+'";'
    df = df[['Seqid', 'Source', 'Type', 'Start',
             'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    df = df.sort_values(['Seqid', 'Start', 'End'],
                        ascending=[True, True, True])
    return df


def make_tsv(df, filename):
    '''
    Input: df
    Output: .tsv file w/o header
    '''
    gff = df.values.tolist()
    with open(filename, 'w') as f:
        for line in gff:
            for word in line:
                f.write(str(word)+'\t')
            f.write('\n')


def main(argv):
    file_main_pred = argv[0]
    file_second_pred = argv[1]
    label_main_pred = argv[2]
    label_second_pred = argv[3]
    min_overlap = int(argv[4])
    file_output = argv[5]

    df1 = open_pred(file_main_pred)
    df2 = open_pred(file_second_pred)
    df = unite_pred(df1, df2, min_overlap, label_main_pred, label_second_pred)
    make_tsv(df, file_output)


if __name__ == "__main__":
    main(sys.argv[1:])
