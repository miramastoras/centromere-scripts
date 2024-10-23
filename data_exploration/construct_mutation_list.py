#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:46:35 2022

Intended for use on the WashU pedigree
Identifies which centromeres are inherited copies in the pedigree and
constructs a table of mutations inferred from their alignments

@author: Jordan
"""

import sys
import re
import os
import numpy as np
import pandas as pd

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed


def get_matched_chroms(aln_identity_table):
    tab = pd.read_table(aln_identity_table, header = 0)
    
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
                
                m = re.match("^(\w+)_(\d)_(\w+)_(\d)$", tab.columns[max_j])
                parent = m.group(3)
                parent_hap = int(m.group(4))
                if child is None:
                    child = m.group(1)
                assert child == m.group(1)
                
                matched_chrom[(chrom, hap)] = (parent, parent_hap)
                
    return child, matched_chrom

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./construct_mutation_list.py aln_identity.tsv aln_dir > muts.tsv", file = sys.stderr)
        sys.exit(1)
    
    aln_identity_table = sys.argv[1]
    aln_dir = sys.argv[2]
        
    filter_edge_alignments = True
    min_op_len_to_pass_filter = sys.maxsize
    
    tab = pd.read_table(aln_identity_table, header = 0)
    
    child, matched_chrom = get_matched_chroms(aln_identity_table)
                
#    for c in matched_chrom:
#        print(c, matched_chrom[c])
    
    columns = ["chrom", "child_hap", "child_pos", "child_len", "parent", "parent_hap", "parent_pos", "parent_len", "type", "length"]
    print("\t".join(columns))

    for chrom, child_hap in sorted(matched_chrom):
        parent, parent_hap = matched_chrom[(chrom, child_hap)]
        
        filename = os.path.join(aln_dir, "_".join((child, str(child_hap), parent, str(parent_hap), chrom)), "cigar.txt")
        cigar = parse_cigar(open(filename).read().strip())
        
        child_len = 0
        parent_len = 0
        for op, l in cigar:
            if op in ("M", "X", "=", "D"):
                child_len += l
            if op in ("M", "X", "=", "I", "S", "H"):
                parent_len += l

        
        child_pos = 0
        parent_pos = 0
        for op, l in cigar:
            
            # the child is treated like the reference in this pipeline, so you have to "reverse"
            # the orientation of operations
            mut_type = None
            if op == "X":
                mut_type = "SUB"
            elif op == "I":
                mut_type = "DEL"
            elif op == "D":
                mut_type = "INS"
                
            
            if op in ("M", "X", "=", "D"):
                child_op_len = l
            else:
                child_op_len = 0
            if op in ("M", "X", "=", "I", "S", "H"):
                parent_op_len = l
            else:
                parent_op_len = 0
            
            if mut_type is not None:
                if (not filter_edge_alignments or 
                    not (l < min_op_len_to_pass_filter and ((parent_pos == 0 and child_pos == 0) or
                                                            (parent_pos + parent_op_len == parent_len and child_pos + child_op_len == child_len)))):
                    print("\t".join(str(v) for v in (chrom, child_hap, child_pos, child_len, parent, parent_hap, parent_pos, parent_len, mut_type, l)))
            
            parent_pos += parent_op_len
            child_pos += child_op_len
                
    
