#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:12:11 2022

@author: Jordan
"""

import sys
import re

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def aligned_interval(cigar):
    
    i = 0
    
    read_start = 0
    if len(cigar) > 0 and cigar[0][0] in ["H", "S"]:
        read_start += cigar[0][1]
        i += 1
        
    read_end = read_start
    aligned_length = 0
    
    while i < len(cigar):
        op, l = cigar[i]
        
        if op in ["I", "X", "=", "M"]:
            read_end += l
        if op in ["D", "X", "=", "M"]:
            aligned_length += l
        
        i += 1
        
    return read_start, read_end, aligned_length
       
def num_matches(cigar):
    m = 0
    for op, l in cigar:
        if op in ["M", "="]:
            m += l
    return m
    

if __name__ == "__main__":
    
    sam = sys.argv[1]
    
    records = []
    
    contig_lens = {}
    
    for line in open(sam):
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith("@"):
            if line.startswith("@SQ"):
                ms = re.search("SN:(\S+)", line)
                ml = re.search("LN:(\d+)", line)
                contig_lens[ms.group(1)] = int(ml.group(1))
            continue
        
        tokens = line.strip().split("\t")
        
        name = tokens[0]
        flag = int(tokens[1])
        contig = tokens[2]
        ref_start = int(tokens[3]) - 1
        mapq = tokens[4]
        cigar = parse_cigar(tokens[5])
        seq_len = len(tokens[9])
        if contig == '*':
            contig_len = 0
        else:
            contig_len = contig_lens[contig]
        
        if flag & 16:
            strand = '-'
        else:
            strand = '+'
        
        
        read_start, read_end, aligned_length = aligned_interval(cigar)
        ref_end = ref_start + aligned_length
        matches = num_matches(cigar)
        block_length = sum(c[1] for c in cigar)
        
        #               0     1        2           3         4       5       6           7          8        9        10            11 
        records.append((name, seq_len, read_start, read_end, strand, contig, contig_len, ref_start, ref_end, matches, block_length, mapq))
        
    records.sort(key = lambda r : (r[5], r[7], r[8], r[0], r[2], r[3], r))
        
    for record in records:
        print("\t".join(str(v) for v in record))
        