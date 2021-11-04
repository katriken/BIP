#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Finds a BRAKER transcript score cutoff based on the scores of supported proteins, 
applies the score filter, and makes a .gtf file with high-scoring predictions 
and all supported proteins. 
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


def find_cutoff(df_scores, prot_sp):
    '''
    Input:
    - df_scores: df with protein ids and transcript BRAKER scores
    - prot_sp: list of supported proteins
    Output: float, BRAKER score cutoff
    '''
    # Scores of supported proteins
    df_scores1 = df_scores[df_scores['protein'].isin(prot_sp)]
    scores1 = df_scores1['score'].values.tolist()
    # Scores of all proteins
    scores0 = df_scores['score'].values.tolist()
    # Possible cutoff scores: from 0 to 1 with step 0.1
    cut_offs = list(range(0, 101))
    cut_offs = [x/100 for x in cut_offs]
    cut_offs = cut_offs[::-1]
    # ROC curve with true positive rate = sensitivity for supported proteins,
    # and false positive rate = 1 - specificity.
    sn1 = [sn(scores1, x) for x in cut_offs]
    sp1 = [sp(scores0, scores1, x) for x in cut_offs]
    tpr = sn1
    fpr = [1-x for x in sp1]
    for i in range(2, len(tpr)):
        d = (tpr[i]-tpr[i-1])/(fpr[i]-fpr[i-1])
        # break when slope of the ROC curve < 1.
        if d < 1:
            cutoff = cut_offs[i]
            fpr_cut = fpr[i]
            tpr_cut = tpr[i]
            break
    print('BRAKER score cutoff:', cutoff)
    print('TPR at cutoff:', tpr_cut)
    print('FPR at cutoff:', fpr_cut)
    return cutoff


def sn(s1, x):
    '''
    Input: 
    - s1: list of scores of supported proteins
    - x: float, possible cutoff score
    Output: float, sensitivity as fraction of supported proteins with a score
    >= x from all supported proteins.
    '''
    s1_passed = [i for i in s1 if i >= x]
    res = len(s1_passed)/len(s1)
    return res


def sp(s, s1, x):
    '''
    Input: 
    - s: list of scores of all proteins
    - s1: list of scores of supported proteins
    - x: float, possible cutoff score
    Output: float, specificity as fraction of supported proteins with a score
    >= x from all proteins with a score >= x.
    '''
    s1_passed = [i for i in s1 if i >= x]
    s_passed = [i for i in s if i >= x]
    res = len(s1_passed)/len(s_passed)
    return res


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
    file_supp_prot_gtf = argv[1]
    file_tx_scores = argv[2]
    file_output = argv[3]

    # tx scores
    df_scores = pd.read_csv(file_tx_scores, sep='\t',
                            index_col=False, names=['protein', 'score'])
    df_scores.loc[(df_scores.score == '.'), 'score'] = '0.001'
    df_scores['score'] = df_scores['score'].astype(str).astype(float)
    # supported proteins
    df_sp = open_pred(file_supp_prot_gtf)
    prot_sp = list(set(df_sp['protein'].values.tolist()))
    # Find BRAKER score cutoff
    print('Calculating the BRAKER transcript score cutoff...')
    cutoff = find_cutoff(df_scores, prot_sp)
    # Make .gtf file with new prediciton: high-scoring + supported proteins
    df_high = df_scores[df_scores['score'] >= cutoff]
    prot_high = df_high['protein'].values.tolist()
    df_low = df_scores[df_scores['score'] < cutoff]
    prot_low = df_low[df_low['protein'].isin(
        prot_sp)]['protein'].values.tolist()
    prot = prot_high+prot_low
    prot = list(set(prot))
    print('High scoring proteins:', len(prot_high))
    print('Low scoring sp proteins:', len(prot_low))
    print('Final:', len(prot))
    save_new_gtf(prot, file_prediction, file_output)


if __name__ == "__main__":
    main(sys.argv[1:])
