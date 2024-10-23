#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 12:06:21 2024

Intended for the WashU pedigree
Identify putative mutations/misassemblies as differences between assemblies
of centromeres that we've traced through the pedigree

@author: Jordan
"""

import sys
import os
import re

def parse_cigar(cigar_str):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar_str):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def ref_len(cigar):
    ref = 0
    for op, l in cigar:
        if op != "I":
            ref += l
    return ref

def reverse_cigar(cigar):
    for i in range(len(cigar)):
        if cigar[i][0] == "I":
            cigar[i] = ("D", cigar[i][1])
        elif cigar[i][0] == "D":
            cigar[i] = ("I", cigar[i][1])

if __name__ == "__main__":
    
    induced_prefix = sys.argv[1]
    ref_seq_prefix = sys.argv[2]
    
    # clear out point variants with a neighborhood of an SV equal to X times its length 
    rad_multiplier = 5.0
    max_window = 30000
    do_filter = False
    boundary_buffer = 10
    
    queries = set()
    
    muts = {}
    
    refl = -1
    
    for fp in os.listdir(os.path.dirname(induced_prefix)):
        if fp.startswith(os.path.basename(induced_prefix)):
            if ref_seq_prefix not in fp:
                continue
            
            cigar = parse_cigar(open(os.path.join(os.path.dirname(induced_prefix), fp)).read().strip())
            m = re.match(os.path.basename(induced_prefix) + "_([^_]+)(_rc)?_([^_]+)(_rc)?\.txt", fp)
            assert(m is not None)
            query = m.group(3)
            
            if query.startswith(ref_seq_prefix):
                reverse_cigar(cigar)
                query = m.group(1)
            refl = ref_len(cigar)
            
            queries.add(query)
            
            p = 0
            for op, l in cigar:
                if op == "X" or op == "I" or op == "D":
                    mut = (p, op, l)
                    if mut not in muts:
                        muts[mut] = []
                    muts[mut].append(query)
                        
                
                if op != "I":
                    p += l
     
    query_col = {q : i for i, q in enumerate(sorted(queries))}
    
    table  = []
    
    
    print("\t".join(str(v) for v in ["pos", "edit", "type"] + sorted(queries)))
    for mut in sorted(muts):
        pos, op, l = mut
        row = [0 for q in range(len(queries) + 3)]
        row[0] = pos
        row[1] = l
        row[2] = op
        for q in muts[mut]:
            row[query_col[q] + 3] = 1
        table.append(row)
    
    
    keep = [True for i in range(len(table))]
    for i in range(len(table)):
        row = table[i]
        l = row[1]
        if row[0] < boundary_buffer or ((row[2] == "I" and row[0] > refl - boundary_buffer) or 
                                        (row[2] != "I" and row[0] + l > refl - boundary_buffer)):
            keep[i] = False
        if l >= 50 and do_filter:
            min_pos = row[0] - min(max_window, rad_multiplier * l)
            max_pos = row[0] + min(max_window, rad_multiplier * l)
            if row[2] == 'D':
                max_pos += l
            
            j = i - 1
            while j >= 0 and table[j - 1][0] >= min_pos:
                j -= 1
            
            while j < len(table) and table[j][0] < max_pos:
                row_len = table[j][1]
                if row_len <= 15:
                    keep[j] = False
                
                j += 1
            
    
    table = [table[i] for i in range(len(table)) if keep[i]]
    
    for row in table:
        print("\t".join(str(v) for v in row))
            
            
