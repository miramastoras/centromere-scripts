#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:48:48 2023

Count the number and size of small variant (mismatches and indels) relative to
the length of matching sequence

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
    
    if len(sys.argv) != 2:
        print("usage: ./variant_density.py cigar.txt", file = sys.stderr)
        sys.exit(1)
    
    cigar = parse_cigar(open(sys.argv[1]).read().strip())
    
    max_allele_size = 32
    
    i = 0
    total_len = 0
    num_runs = 0
    num_vars = 0
    while i < len(cigar):
        num_runs += 1
        while i < len(cigar) and (cigar[i][0] == '=' or cigar[i][0] == 'M' or cigar[i][1] <= max_allele_size):
            total_len += cigar[i][1]
            if cigar[i][0] != '=' and cigar[i][0] != 'M':
                num_vars += 1
            i += 1
        while i < len(cigar) and cigar[i][0] != '=' and cigar[i][0] != 'M' and cigar[i][1] > max_allele_size:
            i += 1
    
    print("{} variants in {} short variant runs of total length {} for {} variants/run, {} bases/run, and {} bases/variant".format(num_vars, num_runs, total_len, num_vars / num_runs, total_len / num_runs, total_len / num_vars))
    
