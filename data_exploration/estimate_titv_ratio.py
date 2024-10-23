#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:16:38 2024

Estimate transition-transversion ratio using a mismatch table

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print("usage: ./estimate_titv_ratio.py mismatch_table.tsv", file = sys.stderr)
        sys.exit(1)
    
    mismatch_table = sys.argv[1]
    
    transitions = {("C", "T"), ("T", "C"), ("A", "G"), ("G", "A")}
    
    num_transitions = 0
    num_transversions = 0
    for line in open(mismatch_table):
        
        tokens = line.strip().split()
        assert(len(tokens) == 6)
        
        if (tokens[2], tokens[5]) in transitions:
            num_transitions += 1
        else:
            num_transversions += 1
            
    
    print("transitions: {}".format(num_transitions))
    print("transversions: {}".format(num_transversions))
    print("Ti/Tv ratio: {}".format(num_transitions / num_transversions))
        
