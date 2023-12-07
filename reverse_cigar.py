#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 18:46:12 2023

@author: Jordan
"""

import re
import sys

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def reverse_cigar(parsed):
    for i in range(len(parsed)):
        if parsed[i][0] == "I":
            parsed[i] = ("D", parsed[i][1])
        elif parsed[i][0] == "D":
            parsed[i] = ("I", parsed[i][1])

def to_cigar_string(cigar):
    return "".join(str(l) + op for op, l in cigar)

if __name__ == "__main__":
    
    cig_fp = sys.argv[1]
    cigar = parse_cigar(open(cig_fp).read().strip())
    reverse_cigar(cigar)
    print(to_cigar_string(cigar))