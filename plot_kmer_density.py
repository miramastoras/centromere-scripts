#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 11:03:54 2022

@author: Jordan
"""


import sys
import matplotlib.pyplot as plt

def parse_fasta(fa):
    seq = ""
    with open(fa) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            if line.startswith(">"):
                continue
            seq += line.strip().upper()
    return seq

def shorten(kmer, l = 40):
    if len(kmer) < l:
        return kmer
    else:
        return kmer[:l] + "_trunc"

if __name__ == "__main__":
    
    fasta = sys.argv[1]
    kmer_file = sys.argv[2]
    outplot = sys.argv[3]
    
    print("reading kmers", file = sys.stderr)
    
    kmers = set()
    for line in open(kmer_file):
        if type(line) == bytes:
            line = line.decode("utf-8")
        kmer = line.strip()
        kmers.add(kmer)
        
    print("reading fasta", file = sys.stderr)
    seq = parse_fasta(fasta)
    
    k = len(next(iter(kmers)))
    
    positions = []
    
    print("finding positions", file = sys.stderr)
    for i in range(len(seq) - k + 1):
        if seq[i:i+k] in kmers:
            positions.append(i)
            
    nbins = 250
    plt.hist(positions, nbins)
    plt.savefig(outplot, dpi = 250, format = "png")
    
    