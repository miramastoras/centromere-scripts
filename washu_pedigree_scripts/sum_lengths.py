#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 21:06:23 2022

Intended for use with the WashU pedigree
Add up the total size of identified centromeres

@author: Jordan
"""

import sys
import os
import pandas as pd
import numpy as np

def parse_fasta(fa):
    seq = ""
    with open(fa) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            if line.startswith(">"):
                continue
            seq += line.strip().upper()
    return seq


if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./sum_lengths.py hor_assignment_table.tsv sequences_dir", file = sys.stderr)
        exit(1)
    
    assignment_table = sys.argv[1]
    seq_dir = sys.argv[2]
    
    tab = pd.read_table(assignment_table, header = 0)
    
    skip = set([("chr13", 2), ("chr15", 2), ("chr17", 1)])
    
    columns = ["chrom", "hap", "parent_len"]
    print("\t".join(columns))
    total_len = 0
    for i in range(tab.shape[0]):
        if np.isnan(tab.parent_haplotype.values[i]):
            continue
        chrom = str(tab.chrom.values[i])
        hap = int(tab.haplotype.values[i])
        parent = str(tab.parent.values[i])
        parent_hap = int(tab.parent_haplotype.values[i])
        
        if (chrom, hap) in skip:
            continue
        
        fasta = os.path.join(seq_dir, "{}_{}_active_arrays".format(parent, parent_hap), chrom + ".fasta")
        seq = parse_fasta(fasta)
        
        print("\t".join(str(v) for v in [tab.chrom.values[i], tab.haplotype.values[i], len(seq)]))
        total_len += len(seq)
        
    print("\t".join(str(v) for v in ["TOTAL", "", total_len]))
        
    
