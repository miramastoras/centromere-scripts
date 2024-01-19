#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 17:07:44 2022

@author: Jordan
"""

import sys
import re
import pandas as pd
import numpy as np

def get_matched_chroms(tab):
    
    matched_chrom = {}
    
    child = None
    for i in range(tab.shape[0]):
        chrom = tab.chr.values[i]
        for hap in (1, 2):
            max_j = None
            max_identity = None
            for j in range(4 * (hap - 1) + 1, 4 * hap + 1):
                identity = tab.iloc[i, j]
                if not np.isnan(identity) and (max_identity is None or max_identity < identity):
                    max_identity = identity
                    max_j = j
            
            if max_identity is not None and max_identity > 0.9:
                
                m = re.match("^(\w+)_(\w+)_(\w+)_(\w+)$", tab.columns[max_j])
                parent = m.group(3)
                parent_hap = int(m.group(4))
                if child is None:
                    child = m.group(1)
                assert child == m.group(1)
                
                matched_chrom[(chrom, hap)] = (parent, parent_hap)
                
    return child, matched_chrom

    
if __name__ == "__main__":
    
    aln_identity_table = sys.argv[1]
    
    tab = pd.read_table(aln_identity_table, header = 0)
    
    child, matched_chrom = get_matched_chroms(tab)
    
    parents = sorted(set(par for par, hap in matched_chrom.values()))
    assert(len(parents) == 2)
    
    columns = ["chrom", "haplotype", "parent", "parent_haplotype"]
    
    print("\t".join(columns))
    for chrom in tab.chr.values:
        for hap in (1, 2):
            if (chrom, hap) in matched_chrom:
                parent, parent_hap = matched_chrom[(chrom, hap)]
            else:
                parent_hap = "NA"
                other_hap = 3 - hap
                if (chrom, other_hap) in matched_chrom:
                    other_parent, other_parent_hap = matched_chrom[(chrom, other_hap)]
                    if other_parent == parents[0]:
                        parent = parents[1]
                    else:
                        parent = parents[0]
                else:
                    parent = "NA"
                    
            
            print("\t".join(str(v) for v in [chrom, hap, parent, parent_hap]))
