#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 19:55:13 2022

@author: Jordan
"""

import sys
import collections

if __name__ == "__main__":
    
    fasta = sys.argv[1]
    k = int(sys.argv[2])
    
    seq = ""
    with open(sys.argv[1]) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            if line.startswith(">"):
                continue
            seq += line.strip().upper()
            
    cntr = collections.Counter()
    
    for i in range(len(seq) - k + 1):
        cntr[seq[i:i+k]] += 1
    
    for kmer in sorted(cntr):
        print("{}\t{}".format(kmer, cntr[kmer]))