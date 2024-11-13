#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:25:21 2022

@author: Jordan
"""

import sys


def print_rc(seq):
    comp = {"A":"T", "C":"G", "G":"C", "T":"A","N":"N"}
    rc = "".join(comp[c] for c in reversed(seq))
    for i in range(0, len(seq), 80):
        print(rc[i:i+80])

if __name__ == "__main__":
    
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("usage: ./fasta_to_rev_comp.py seq.fasta > rev.fasta\nor./fasta_to_rev_comp.py < seq.fasta > rev.fasta\n", file = sys.stderr)
        sys.exit(1)
    
    f = None
    if len(sys.argv) == 1:
        f = sys.stdin
    elif len(sys.argv) == 2:
        f = open(sys.argv[1])
    else:
        print("error: invalid input", file = sys.stderr)
        sys.exit()
    
    seq = ""
    for line in f:
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith(">"):
            if seq != "":
                print_rc(seq)
                seq = ""
            print(line.strip() + " rc")
        else:
            seq += line.strip().upper()
            
    print_rc(seq)
