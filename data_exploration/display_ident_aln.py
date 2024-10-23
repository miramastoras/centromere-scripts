#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 13:35:19 2023

Print out the truth identity of two simulated sequences and an alignment of them
for manual/visual comparison

@author: Jordan
"""

import sys
import re


def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print("usage: display_ident_aln.py identity1.txt identity2.txt aln.txt", file = sys.stderr)
        sys.exit(1);
    
    ident1_fp = sys.argv[1]
    ident2_fp = sys.argv[2]
    cigar_fp = sys.argv[3]
    
    ident1 = [int(v) for v in open(ident1_fp).read().strip().split()]
    ident2 = [int(v) for v in open(ident2_fp).read().strip().split()]
    
    cigar = parse_cigar(open(cigar_fp).read().strip())
    
    
    print("\t".join(["i", "j", "id1", "id2"]))
    i = 0
    j = 0
    for op, l in cigar:
        if op == "M" or op == "X" or op == "=":
            for k in range(l):
                
                print("\t".join(str(v) for v in [i, j, ident1[i], ident2[j]]))
                i += 1
                j += 1
        elif op == "D":
            for k in range(l):
                
                print("\t".join(str(v) for v in [i, j, ident1[j], "-"]))
                i += 1
        elif op == "I":
            for k in range(l):
                
                print("\t".join(str(v) for v in [i, j, "-", ident2[j]]))
                j += 1
    
