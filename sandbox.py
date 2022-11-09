#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:52:48 2022

@author: Jordan
"""

def split_to_approx_monomers(hor):
    n = int(round(len(s) / 171))
    
    monomers = []
    for i in range(n):
        monomers.append(hor[(i * len(hor)) // n : ((i+1) * len(hor)) // n])
    return monomers

if __name__ == "__main__":
    
    seqs = []
    s = ""
    for line in open("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/stv_seqs/chm_stv.fasta"):
        if line.startswith(">"):
            if s != "":
                seqs.append(s)
            s = ""
        else:
            s += line.strip()
            