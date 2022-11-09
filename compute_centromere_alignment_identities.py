#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:18:12 2022

@author: Jordan
"""

import sys
import re
import os

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def compute_identity(cigar):
    
    seq1_len = 0
    seq2_len = 0
    match_len = 0
    
    for op, l in cigar:
        # FIXME: tandem aligner makes the classic M/= misuse of CIGAR
        if op == "M" or op == "=":
            seq1_len += l
            seq2_len += l
            match_len += l
        elif op == "I":
            seq1_len += l
        elif op == "D":
            seq2_len += l
        elif op == "X":
            seq1_len += l
            seq2_len += l
        else:
            assert(False)
        
    return float(match_len) / max(seq1_len, seq2_len)

def parse_aln_dir_name(dir_name):
    
    match = re.match("^(\w+)_([12])_(\w+)_([12])_([\w\.]+)$", dir_name)
    
    
    return (match.group(i) for i in range(1, 6))
    

if __name__ == "__main__":
    
    aln_dir = sys.argv[1]
    
    aln_output_dirs = os.listdir(aln_dir)
    
    chroms = set()
    parents = set()
    for output_dir in aln_output_dirs:
        child, child_hap, parent, parent_hap, chrom = parse_aln_dir_name(output_dir)
        chroms.add(chrom)
        parents.add(parent)
    
    chroms = sorted(chroms)
    parents = sorted(parents)
    
    chrom_idx = {chroms[i] : i for i in range(len(chroms))}
    
    table = [["NA" for j in range(4 * len(parents))] for i in range(len(chroms))]
    
    for output_dir in aln_output_dirs:
        child, child_hap, parent, parent_hap, chrom = parse_aln_dir_name(output_dir)
        
        cigar_path = os.path.join(aln_dir, output_dir, "cigar.txt")
        
        cigar = parse_cigar(open(cigar_path).read().strip())
        
        identity = compute_identity(cigar)
        
        i = chrom_idx[chrom]
        j = 4 * "12".index(child_hap) + 2 * parents.index(parent) + "12".index(parent_hap)
        
        table[i][j] = str(identity)
        
    
    for i in range(len(chroms)):
        table[i].insert(0, chroms[i])
    
    columns = ["chr"]
    for c_hap in range(1, 3):
        for parent in parents:
            for p_hap in range(1, 3):
                columns.append("{}_{}_{}_{}".format(child, c_hap, parent, p_hap))
                
    print("\t".join(columns))
    for row in table:
        print("\t".join(row))
        
        
        
        
        
        
    
    
    
        