#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 18:01:16 2023

@author: Jordan
"""

import sys
import os
import re
import skbio

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def cigar_to_dist(cigar):
    
    query_len = 0
    ref_len = 0
    matches = 0
    for op, op_len in cigar:
        if op in "MX=":
            query_len += op_len
            ref_len += op_len
            if op != "X":
                matches += op_len
        elif op == "D":
            ref_len += op_len
        elif op in "IHS":
            query_len += op_len
        else:
            assert(False)
    
    return 1.0 - (2.0 * matches) / (ref_len + query_len)

if __name__ == "__main__":
    
    aln_dir = sys.argv[1]
    
    mat = {}
    
    for fp in os.listdir(aln_dir):
        with open(os.path.join(aln_dir, fp)) as f:
            cigar = parse_cigar(f.read().strip())
        m = re.match("([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_cigar.txt", fp)
        assert(m is not None)
        samp1 = m.group(1)
        samp2 = m.group(2)
        dist = cigar_to_dist(cigar)
        mat[(samp1, samp2)] = dist
        mat[(samp2, samp1)] = dist
        
    samps = sorted(set(s[0] for s in mat))
    
    # reorganize as an array
    D = []
    for samp1 in samps:
        D.append([])
        for samp2 in samps:
            if samp1 == samp2:
                D[-1].append(0.0)
            else:
                D[-1].append(mat[(samp1, samp2)])
    
    
    # make skbio type
    dist_mat = skbio.DistanceMatrix(D, samps)
    print(dist_mat.to_data_frame(), file = sys.stderr)
    
    tree = skbio.tree.nj(dist_mat)
    tree = tree.root_at_midpoint()
    
    print(tree.ascii_art(), file = sys.stderr)
    print(tree)
        