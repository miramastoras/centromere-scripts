#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 20:29:26 2022

Script for visualizing k-mer matches between two sequences

@author: Jordan
"""

import sys
import re
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
    
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    kmers = sys.argv[3]
    outdir = sys.argv[4]    
            
    seq_height = .25
    y1 = 1.0 - seq_height
    y2 = 0.0
    kmer_width = .001
    resolution = 250
    
    
    seq1 = parse_fasta(fasta1)
    seq2 = parse_fasta(fasta2)
    
    for line in open(kmers):
        if type(line) == bytes:
            line = line.decode("utf-8")
        kmer = line.strip()
        
        len1 = len(seq1)
        len2 = len(seq2)
        
        positions1 = [m.span()[0] for m in re.finditer(kmer, seq1)]
        positions2 = [m.span()[0] for m in re.finditer(kmer, seq2)]
        
        f, ax = plt.subplots(1, 1, figsize = (8, 2), dpi = resolution)
        ax.set_xlim(-.05, 1.05)
        ax.set_ylim(-.05, 1.05)
        
        relative_len1 = float(len1) / max(len1, len2)
        relative_len2 = float(len2) / max(len1, len2)
        x1 = 0.5 - 0.5 * relative_len1
        x2 = 0.5 - 0.5 * relative_len2
        
        seq_rect1 = patches.Rectangle((x1, y1), relative_len1, seq_height, linewidth=0, edgecolor=None, facecolor='k')
        seq_rect2 = patches.Rectangle((x2, y2), relative_len2, seq_height, linewidth=0, edgecolor=None, facecolor='k')
        
        ax.add_patch(seq_rect1)
        ax.add_patch(seq_rect2)
        
        for pos1 in positions1:
            xp = x1 + relative_len1 * (float(pos1) / len1) - .5 * kmer_width
            pos_rect = patches.Rectangle((xp, y1), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
            ax.add_patch(pos_rect)
        for pos2 in positions2:
            xp = x2 + relative_len2 * (float(pos2) / len2) - .5 * kmer_width
            pos_rect = patches.Rectangle((xp, y2), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
            ax.add_patch(pos_rect)
            
        for pos1 in positions1:
            xp1 = x1 + relative_len1 * (float(pos1) / len1) - .5 * kmer_width
            for pos2 in positions2:
                xp2 = x2 + relative_len2 * (float(pos2) / len2) - .5 * kmer_width
                ax.plot([xp1, xp2], [y1, y2 + seq_height], "r-,", linewidth = 1)
        
        
        f.savefig(os.path.join(outdir, shorten(kmer) + ".png"), format = "png", bbox_inches = "tight")
        
        plt.close(f)
        
