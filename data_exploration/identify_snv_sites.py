#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 17:26:09 2024

Identify SNVs as alignment mismatches using the intermediate subalignments that
centrolign can optionally produce as output
Does some heuristic filtering for artifacts

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print("usage: ./identify_snv_sites.py subalns.txt", file = sys.stderr)
        sys.exit(1)
    
    sub_aln_fp = sys.argv[1]
    
    indel_sv_min_size = 50
    radius_multiplier = 5.0
    max_radius = 50000
    do_filter = True
    boundary_buffer = 20
    
    
    table = []
    
    sub_aln_group = 0
    aln_pos_1 = 0
    aln_pos_2 = 0
    curr_indel_len = 0
    min_unblocked_anti_diag = boundary_buffer
    aln_buffer = []
    for line in open(sub_aln_fp):
        if line.startswith("#"):
            sub_aln_group = (sub_aln_group + 1) % 3
            
            if aln_pos_1 != 0 or aln_pos_2 != 0:
                while len(aln_buffer) != 0 and aln_buffer[-1][0] >= aln_pos_1 + aln_pos_2 - boundary_buffer:
                    aln_buffer.pop() 
                for r in aln_buffer:
                    table.append(r[1:])            
                
                aln_pos_1 = 0
                aln_pos_2 = 0
                curr_indel_len = 0
                min_unblocked_anti_diag = boundary_buffer
                aln_buffer = []
            
        elif sub_aln_group == 0:
            
            seq1, pos1, base1, seq2, pos2, base2 = line.strip().split()
            
            if seq1 == "-":
                aln_pos_2 += 1
                curr_indel_len += 1
            elif seq2 == "-":
                aln_pos_1 += 1
                curr_indel_len += 1
            else:
                if curr_indel_len > indel_sv_min_size:
                    
                    radius = min(round(radius_multiplier * curr_indel_len), max_radius)
                    
                    while len(aln_buffer) != 0 and aln_buffer[-1][0] >= aln_pos_1 + aln_pos_2 - curr_indel_len - radius:
                        aln_buffer.pop()
                    
                    min_unblocked_anti_diag = max(min_unblocked_anti_diag, aln_pos_1 + aln_pos_2 + radius)
                    
                curr_indel_len = 0
                
                if aln_pos_1 + aln_pos_2 >= min_unblocked_anti_diag and base1 != base2:
                    # we're at a 
                    aln_buffer.append((aln_pos_1 + aln_pos_2, seq1, pos1, base1, seq2, pos2, base2))
                
                aln_pos_1 += 1
                aln_pos_2 += 1
    
    while len(aln_buffer) != 0 and aln_buffer[-1][0] >= aln_pos_1 + aln_pos_2 - boundary_buffer:
        aln_buffer.pop()  
    
    for r in aln_buffer:
        table.append(r[1:])
    
    for r in table:
        print("\t".join(str(v) for v in r))
