#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:45:58 2022

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    assert(len(sys.argv) == 1 or len(sys.argv) == 2)
    if len(sys.argv) == 1:
        bed = sys.stdin
    else:
        bed = open(sys.argv[1])
        
    prev_chrom = None
    prev_pos = None
    for line in bed:
        if type(line) == bytes:
            line = line.decode("utf-8")
            
        if line.startswith("track"):
            continue
        
        tokens = line.strip().split("\t")
        b = int(tokens[1])
        e = int(tokens[2])
        l = str(e - b)
        
        if prev_chrom != tokens[0]:
            gap = ""
        else:
            gap = str(b - prev_pos)
            
        tokens.insert(3, gap)
        tokens.insert(3, l)
        
        print("\t".join(tokens))
        
        prev_chrom = tokens[0]
        prev_pos = e