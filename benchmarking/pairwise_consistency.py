#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:37:40 2023

Quantify the consistency of direct pairwise alignments and pairwise alignments
induced from a multiple sequence alignment

@author: Jordan
"""

import sys
import os
import re

def map_to_samples(fps):
    suff_regex = "([a-zA-Z0-9]+.[0-9])_([a-zA-Z0-9]+.[0-9]).txt$"

    sample_map = {}
    for fp in fps:
        m = re.search(suff_regex, fp)
        assert(m is not None)
        key = tuple([m.group(1), m.group(2)])
        sample_map[key] = fp
    
    return sample_map

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

def cigar_to_positions(cigar):
    
    positions = []
    
    r = 0
    q = 0
    for op, op_len in cigar:
        if op in "MX=":
            for k in range(op_len):
                positions.append((r,q))
                r += 1
                q += 1
        elif op == "D":
            for k in range(op_len):
                positions.append((r,-1))
                r += 1
        elif op in "IHS":
            for k in range(op_len):
                positions.append((-1,q))
                q += 1
        else:
            assert(False)
    
    return positions, r, q



if __name__ == "__main__":
    
    
    if len(sys.argv) != 3:
        print("usage:\npairwise_consistency.py induced_prefix direct_prefix > consistency.txt")
        exit(1)
    
    induced_prefix = os.path.abspath(sys.argv[1])
    direct_prefix = os.path.abspath(sys.argv[2])
    
    induced_dir = os.path.dirname(induced_prefix)
    direct_dir = os.path.dirname(direct_prefix)
    
    induced_file_prefix = os.path.basename(induced_prefix)
    direct_file_prefix = os.path.basename(direct_prefix)
    
    induced_files = [os.path.join(induced_dir, f) for f in os.listdir(induced_dir) if f.startswith(induced_file_prefix)]
    direct_files = [os.path.join(direct_dir, f) for f in os.listdir(direct_dir) if f.startswith(direct_file_prefix)]
    
    samples_to_induced = map_to_samples(induced_files)
    samples_to_direct = map_to_samples(direct_files)
    
    header = ["sample1", "sample2", "intersection", "union", "aligned_intersection", "aligned_union", "jaccard", "aligned_jaccard", "num_pos_ind", "num_pos_dir"]
    print("\t".join(header))
    for pair in sorted(samples_to_induced, key = lambda k : sorted(k)):
        
        induced_cigar = parse_cigar(open(samples_to_induced[pair]).read())

        if pair not in samples_to_direct:
            rev_pair = (pair[1], pair[0])
            assert(rev_pair in samples_to_direct)
            direct_cigar = parse_cigar(open(samples_to_direct[rev_pair]).read())
            reverse_cigar(direct_cigar)
            
        else:
            direct_cigar = parse_cigar(open(samples_to_direct[pair]).read())
        
        induced_positions, rl1, ql1 = cigar_to_positions(induced_cigar)
        direct_positions, rl2, ql2 = cigar_to_positions(direct_cigar)
        induced_positions = set(induced_positions)
        direct_positions = set(direct_positions)
        
        # ensure that we have consistent assignment of ref and query
        assert(rl1 == rl2)
        assert(ql1 == ql2)
        
        induced_aligned_positions = set(p for p in induced_positions if p[0] != -1 and p[1] != -1)
        direct_aligned_positions = set(p for p in direct_positions if p[0] != -1 and p[1] != -1)
        
        inter = len(induced_positions.intersection(direct_positions))
        union = len(induced_positions.union(direct_positions))
        inter_aligned = len(induced_aligned_positions.intersection(direct_aligned_positions))
        union_aligned = len(induced_aligned_positions.union(direct_aligned_positions))
        
        jacc = inter / union
        if union_aligned == 0:
            jacc_aln = "NA"
        else:
            jacc_aln = inter_aligned / union_aligned
        
        
        row = [pair[0], pair[1], inter, union, inter_aligned, union_aligned, jacc, jacc_aln, len(induced_positions), len(direct_positions)]
        print("\t".join(str(v) for v in row))
    
