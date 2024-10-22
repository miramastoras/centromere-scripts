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

def parse_aln_name(filename):
    
    match = re.search("(PAN\d+)_haplotype([12])_(PAN\d+)_haplotype([12]).(chr\w+)", filename)
    
    return (match.group(i) for i in range(1, 6))
    

if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print("usage:\ncompute_new_centromere_alignment_identities.py aln_dir")
        exit(1)
    
    aln_dir = sys.argv[1]
        
    chroms = set()
    children = set()
    parents = set()
    child_haps = set()
    parent_haps = {}
    for aln in os.listdir(aln_dir):
        
        child, child_hap, parent, parent_hap, chrom = parse_aln_name(aln)
        
        chroms.add(chrom)
        parents.add(parent)
        child_haps.add(child_hap)
        if parent not in parent_haps:
            parent_haps[parent] = set()
        parent_haps[parent].add(parent_hap)
        children.add(child)
    
    chroms = sorted(chroms)
    child_haps = sorted(child_haps)
    parents = sorted(parents)
    for parent in parent_haps:
        parent_haps[parent] = sorted(parent_haps[parent])
    children = sorted(children)
    
    assert(len(child_haps) == 2)
    assert(len(children) == 1)
    # assert(len(parents) == 2)
    # assert(len(parent_haps) == 2)
    
    chrom_idx = {chroms[i] : i for i in range(len(chroms))}
    
    table = [["NA" for j in range(4 * len(parents))] for i in range(len(chroms))]
    
    for aln in os.listdir(aln_dir):
        
        child, child_hap, parent, parent_hap, chrom = parse_aln_name(aln)
        
        cigar_path = os.path.join(aln_dir, aln)
        
        cigar = parse_cigar(open(cigar_path).read().strip())
        
        identity = compute_identity(cigar)
        
        i = chrom_idx[chrom]
        j = 2 * len(parents) * child_haps.index(child_hap) + 2 * parents.index(parent) + parent_haps[parent].index(parent_hap)
        
        table[i][j] = str(identity)
        
    
    for i in range(len(chroms)):
        table[i].insert(0, chroms[i])
    
    columns = ["chr"]
    for c_hap in child_haps:
        for parent in parents:
            for p_hap in parent_haps[parent]:
                columns.append("{}_{}_{}_{}".format(children[0], c_hap, parent, p_hap))
                
    print("\t".join(columns))
    for row in table:
        print("\t".join(row))
        
        
        
        
        
        
    
    
    
        
