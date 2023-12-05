#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:53:20 2023

@author: Jordan
"""

import sys
import os
import re

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSNMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def to_simple_cigar(cigar):
    simple = []
    for m in re.finditer("(\d+)([HSNMIDX=])", cigar):
        op = m.group(2)
        l = int(m.group(1))
        if op == "=" or op == "X":
            op = "M"
        if len(simple) != 0 and op == simple[-1][1]:
            simple[-1][0] += l
        else:
            simple.append([l, op])
    return "".join(str(l) + op for l, op in simple)
        

def count_aligned(cigar):
    aligned = 0
    for op, l in cigar:
        if op == "M" or op == "X" or op == "=":
            aligned += l
            
    return aligned

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: estimate_nuc_diversity.py snp_mat.tsv pairwise_aln_prefix")
        exit(1)
        

    prefix = os.path.abspath(sys.argv[2])
    cigar_dir = os.path.dirname(prefix)
    cigar_file_prefix = os.path.basename(prefix)

    cigar_files = [os.path.join(cigar_dir, f) for f in os.listdir(cigar_dir) if f.startswith(cigar_file_prefix)]


    total_aligned = 0
    for cigar_file in cigar_files:
        with open(cigar_file) as f:
            cigar = parse_cigar(f.read().strip())
            total_aligned += count_aligned(cigar)

    
    column_count = None
    alleles = None
    
    with open(sys.argv[1]) as f:
        
        # discard the header
        next(f)
        
        # parse the data
        for line in f:
            tokens = line.strip().split()
            # init the 
            if column_count is None:
                # records of (allele1 count, allele2 count)
                column_count = [[0, 0] for i in range(len(tokens) - 1)]
                alleles = [[None, None] for i in range(len(tokens) - 1)]
            # skip the sample identifier
            for i in range(1, len(tokens)):
                if tokens[i] != "?":
                    # this SNP is observed
                    
                    # make sure we know the alleles
                    if alleles[i - 1][0] is None:
                        alleles[i - 1][0] = tokens[i]
                    elif alleles[i - 1][1] is None and tokens[i] != alleles[i - 1][0]:
                        alleles[i - 1][1] = tokens[i]
                        
                    # add count to the allele
                    if tokens[i] == alleles[i - 1][0]:
                        column_count[i - 1][0] += 1
                    else:
                        column_count[i - 1][1] += 1
                        
    
    transitions = {("C", "T"), ("T", "C"), ("A", "G"), ("G", "A")}
    
    num_transitions = 0
    num_transversions = 0
    for i in range(len(column_count)):
        count1, count2 = column_count[i]
        num_diffs = count1 * count2
        if tuple(alleles[i]) in transitions:
            num_transitions += num_diffs
        else:
            num_transversions += num_diffs
                    
    print("total aligned pairs: {}".format(total_aligned))
    print("total nucleotide differences: {}".format(num_transitions + num_transversions))
    print("nucleotide diversity: {}".format(float(num_transitions + num_transversions) / total_aligned))
    print("Ti/Tv ratio: {}".format(float(num_transitions) / num_transversions))
    
        