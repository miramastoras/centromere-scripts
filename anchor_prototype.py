#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 11:27:12 2022

@author: Jordan
"""

import sys
import collections
import bisect

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

def longest_increasing_subsequence(arr):
    
    # we could actually use only dp_idx, but python's bisect doesn't allow
    # a key function
    dp = [-1]
    dp_idx = [len(arr)]
    pred = [None for i in range(len(arr))]
    
    for i in range(len(arr)):
        p = bisect.bisect_right(dp, arr[i])
        if p == len(dp):
            dp.append(arr[i])
            dp_idx.append(i)
        elif dp[p - 1] < arr[i] and arr[i] < dp[p]:
            dp[p] = arr[i]
            dp_idx[p] = i
        pred[i] = dp_idx[p - 1]
    
    traceback = [dp_idx[-1]]
    while pred[traceback[-1]] != len(arr):
        traceback.append(pred[traceback[-1]])
        
    return traceback[::-1]

def count_filter(c1, c2, filt):
    if filt < 0:
        return True
    elif filt == 0:
        return c1 > 1 or c2 > 1
    else:
        return max(c1, c2) > 1 and min(c1, c2) <= filt

if __name__ == "__main__":
    
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    k = int(sys.argv[3])
    # <0 : all kmers, 0: all multi, n>0, multi with min occurences <= n
    multi_filter = int(sys.argv[4])
    #only_multi = bool(int(sys.argv[4]))
    rejected_pairs = None
    if len(sys.argv) >= 6:
        rejected_pairs = sys.argv[5]
    dump_all = False
    if len(sys.argv) >= 7:
        dump_all = bool(int(sys.argv[6]))
    unique_pairs = None
    if len(sys.argv) >= 8:
        unique_pairs = open(sys.argv[7], "w")
    
    print("loading FASTAs", file = sys.stderr)
    
    seq1 = parse_fasta(fasta1)
    seq2 = parse_fasta(fasta2)
    
    
    print("counting k-mers", file = sys.stderr)
    
    seq1_pos = collections.defaultdict(list)
    for i in range(len(seq1) - k + 1):
        seq1_pos[seq1[i:i+k]].append(i)
        
    seq2_pos = collections.defaultdict(list)
    for j in range(len(seq2) - k + 1):
        seq2_pos[seq2[j:j+k]].append(j)
    

    print("constructing bipartite graph", file = sys.stderr)
    
    graph = [[] for i in range(len(seq1) - k + 1)]
    
    it = 0
    
    
    num_edges = 0
    for kmer in seq1_pos:
        if kmer in seq2_pos:
            c1 = len(seq1_pos[kmer])
            c2 = len(seq2_pos[kmer])
            
#            if it <= 50:
#                print("c1 {} c2 {} om {}".format(c1, c2, only_multi), file = sys.stderr)
#                it += 1
            
            if count_filter(c1, c2, multi_filter):
#            if (only_multi and ((c1 == 1) + (c2 == 1)) == 1) or ((not only_multi) and (c1 == 1 or c2 == 1)):
                for i in seq1_pos[kmer]:
                    for j in seq2_pos[kmer]:
                        graph[i].append(j)
                        num_edges += 1
            if unique_pairs is not None and c1 == 1 and c2 == 1:
                print("{}\t{}".format(seq1_pos[kmer][0], seq2_pos[kmer][0]), file = unique_pairs)
            
                        
#    print("added {} edges".format(num_edges), file = sys.stderr)
    
    
    if dump_all:
        print("outputting graph", file = sys.stderr)
        for i in range(len(graph)):
            for j in graph[i]:
                print("{}\t{}".format(i, j))
        
        sys.exit(0)
    
    print("constructing permutation graph", file = sys.stderr)

    
    in_degree = [0 for j in range(len(seq2) - k + 1)]
    for l in graph:
        for j in l:
            in_degree[j] += 1
            
    
#    print("average in-degree: {}".format(sum(in_degree) / len(in_degree)), file = sys.stderr)
            
    ranges_top = [0]
    for l in graph:
        ranges_top.append(ranges_top[-1] + len(l))
    
    ranges_bottom = [0]
    for d in in_degree:
        ranges_bottom.append(ranges_bottom[-1] + d)
        
    
    next_bottom = ranges_bottom[:]
    
    permutation = [None for i in range(ranges_top[-1])]
    
    
    print("permutation graph has {} nodes".format(len(permutation)), file = sys.stderr)
    
    p = 0
    for i in range(len(graph)):
        for j in graph[i]:
            permutation[p] = next_bottom[j]
            next_bottom[j] += 1
            p += 1
            
    print("minimizing edge crossings", file = sys.stderr)

            
    
    edges = []
    for e in longest_increasing_subsequence(permutation):
        i = bisect.bisect(ranges_top, permutation[e]) - 1
        j = bisect.bisect(ranges_bottom, permutation[e]) - 1
        print("{}\t{}".format(i, j))
        
    
    if rejected_pairs is not None:
    
        edge_set = set(edges)
        
        with open(rejected_pairs, "w") as f:
        
            for i in range(len(graph)):
                for j in graph[i]:
                    if (i, j) not in edge_set:
                        print("{}\t{}".format(i, j), file = f)
        
        
        