#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:25:45 2024

@author: Jordan
"""

import sys
import pandas as pd
import re
import os

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage:\nidentify_pedigree_alignment_problems.py identity_table.txt output_prefix", file = sys.stderr)
        exit(1)
        
    
    tab = pd.read_table(sys.argv[1])
    out_prefix = sys.argv[2]
    if os.path.isdir(out_prefix) and not out_prefix.endswith("/"):
        out_prefix += "/"
    
    # print("dirname " + os.path.dirname(out_prefix))
    # print("basename " + os.path.basename(out_prefix))\\
        
    out_dir = os.path.dirname(out_prefix)
    if out_dir != "" and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if out_dir == "":
        out_dir = "."
    for fp in os.listdir(out_dir):
        if fp.startswith(os.path.basename(out_prefix)):
            print("error: file with prefix already exists: " + os.path.join(os.path.dirname(out_prefix), fp), file = sys.stderr)
            sys.exit(1)
    
    min_ident = 0.8
    
    num_rels = (tab.shape[1] - 1) // 4
    
    assignments = {}
    
    col_regex = "(\w+)_(\w+)_(\w+)_(\w+)"
        
    for i in range(tab.shape[0]):
        chrom = tab.chr.values[i]
        # list of matched chroms for hap1 and hap2
        assignments[chrom] = {}
        for relative in range(num_rels):
            
            hap1rel1 = 1 + 2 * relative
            hap1rel2 = hap1rel1 + 1
            hap2rel1 = hap1rel1 + 2 * num_rels
            hap2rel2 = hap2rel1 + 1
            
            rel_name = re.match(col_regex, tab.columns[hap1rel1]).group(3)
            
            
            if max(tab.iloc[i,hap1rel1], tab.iloc[i,hap2rel2]) > max(tab.iloc[i,hap1rel2], tab.iloc[i,hap2rel1]):
                if tab.iloc[i,hap1rel1] > min_ident:
                    m = re.match(col_regex, tab.columns[hap1rel1])
                    k = (m.group(1), m.group(2))
                    if k not in assignments[chrom]:
                        assignments[chrom][k] = []
                    assignments[chrom][k].append((m.group(3), m.group(4)))
                if tab.iloc[i,hap2rel2] > min_ident:
                    m = re.match(col_regex, tab.columns[hap2rel2])
                    k = (m.group(1), m.group(2))
                    if k not in assignments[chrom]:
                        assignments[chrom][k] = []
                    assignments[chrom][k].append((m.group(3), m.group(4)))
            else:
                if tab.iloc[i,hap1rel2] > min_ident:
                    m = re.match(col_regex, tab.columns[hap1rel2])
                    k = (m.group(1), m.group(2))
                    if k not in assignments[chrom]:
                        assignments[chrom][k] = []
                    assignments[chrom][k].append((m.group(3), m.group(4)))
                if tab.iloc[i,hap2rel1] > min_ident:
                    m = re.match(col_regex, tab.columns[hap2rel1])
                    k = (m.group(1), m.group(2))
                    if k not in assignments[chrom]:
                        assignments[chrom][k] = []
                    assignments[chrom][k].append((m.group(3), m.group(4)))
    
    
    for chrom in assignments:
        for hap in assignments[chrom]:
            row = [chrom, hap[0], hap[1]]
            for assigned in assignments[chrom][hap]:
                row.extend(assigned)
            
            
            with open(out_prefix + "_array_matches_size_{}.txt".format(len(assignments[chrom][hap]) + 1), "a") as f:
                print("\t".join(str(v) for v in row), file = f)
            
            
                
    