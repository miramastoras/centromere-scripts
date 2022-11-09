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
        

if __name__ == "__main__":
    
    sam = sys.argv[1]
    
    records = []
    
    for line in open(sam):
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith("@"):
            continue
        
        tokens = line.strip().split("\t")
        
        name = tokens[0]
        flag = int(tokens[1])
        contig = tokens[2]
        ref_start = int(tokens[3]) - 1
        
        cigar = parse_cigar(tokens[5])
        
        read_start, read_end, aligned_length = aligned_interval(cigar)
        ref_end = ref_start + aligned_length
        
        records.append((contig, ref_start, ref_end, name, read_start, read_end))
        
    records.sort(key = lambda r: (r[4], r[3], r[1], r[2], r[0]))
        
    for contig, ref_start, ref_end, name, read_start, read_end in records:
        print("{}\t{}\t{}\t{}\t{}\t{}".format(contig, ref_start, ref_end, name, read_start, read_end))
        