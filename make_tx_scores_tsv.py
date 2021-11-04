#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes .tsv file with BRAKER tx scores.
Important: BRAKER .gtf file should contain transcript entries.
"""

import pandas as pd
import re
import sys


def find_tx_scores(filename):
    '''
    Input: .gtf file.
    Output: df with protein ids and corresponding scores.
    '''
    df = pd.read_csv(filename, sep='\t', index_col=False, names=[
                     'Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl'])
    df = df[df['Type'] == 'transcript']
    df = df[['Suppl', 'Score']]
    return df


def make_tsv(df, filename):
    '''
    Input: df.
    Output: .tsv file w/o header.
    '''
    tsv = df.values.tolist()
    with open(filename, 'w') as f:
        for line in tsv:
            for word in line:
                f.write(str(word)+'\t')
            f.write('\n')


def main(argv):
    file_braker_gtf = argv[0]
    file_output = argv[1]

    df = find_tx_scores(file_braker_gtf)
    make_tsv(df, file_output)


if __name__ == "__main__":
    main(sys.argv[1:])
