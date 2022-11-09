#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:01:05 2022

@author: Jordan
"""

import sys

mems = "/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/chm_hg002_mems_342.txt"
mums = "/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/chm_hg002_mums_342.txt"

inf = sys.maxsize
    
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

class MatchNode:
    
    def __init__(self, b1, b2, length, node_id):
        self.begin1 = b1
        self.begin2 = b2
        self.length = length
        self.node_id = node_id
        self.edges = []
        
        
def construct_graph_simple(mems):
    
    nodes = []
    
    for i in range(len(mems)):
        v, h, length = mems[i]
        nodes.append(MatchNode(v, h, length, i))
        
    # TODO: i don't think you need an n log n sort here
    nodes.sort(key = lambda n: min(n.begin1, n.begin2))
    
    for i in range(len(nodes)):
        
        node_from = nodes[i]
        
        for j in range(i + 1, len(nodes)):
            
            node_to = nodes[j]
            
            if (node_from.begin1 < node_to.begin1 and node_from.begin2 < node_to.begin2 
                and node_from.begin1 + node_from.length < node_to.begin1 + node_to.length 
                and node_from.begin2 + node_from.length < node_to.begin2 + node_to.length):
                
                node_from.edges.append(j)
                
    
    return nodes


def graph_dp(nodes, match, gap_open, gap_extend):
    
    fwd = [-inf for i in range(len(nodes))]
    bwd = [-inf for i in range(len(nodes))]
    
    is_source = [True for i in range(len(nodes))]
    
    for i in range(len(nodes)):
        node = nodes[i]
        if len(node.edges) == 0:
            bwd[i] = 0
        for j in node.edges:
            is_source[j] = False
            
    
    for i in range(len(nodes)):
        node = nodes[i]
        if is_source[i]:
            fwd[i] = match * node.length
            
        for j in node.edges:
            next_node = nodes[j]
            
            diag_offset = abs((node.begin1 - node.begin2) - (next_node.begin1 - next_node.begin2))
            
            edge_score = 0
            if diag_offset != 0:
                edge_score = diag_offset * gap_extend + gap_open
                
            overlap = max(node.begin1 + node.length - next_node.begin1,
                          node.begin2 + node.length - next_node.begin2)
            
            if overlap > 0:
                edge_score -= overlap * match
                
            fwd[j] = max(fwd[j], fwd[i] + edge_score + match * next_node.length)

    for i in range(len(nodes) - 1, -1, -1):
        node = nodes[i]
        for j in next_node.edges:
            next_node = nodes[j]
            
            diag_offset = abs((node.begin1 - node.begin2) - (next_node.begin1 - next_node.begin2))
            
            edge_score = 0
            if diag_offset != 0:
                edge_score = diag_offset * gap_extend + gap_open
                
            overlap = max(node.begin1 + node.length - next_node.begin1,
                          node.begin2 + node.length - next_node.begin2)
            
            if overlap > 0:
                edge_score -= overlap * match
                
            bwd[i] = max(bwd[i], bwd[j] + edge_score + match * next_node.length)
            
    
    return fwd, bwd
                
            
            
def plot_lines(drawing, filename, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
               line_width, color):
    coord_multiplier = float(long_side) / max(seq1_len, seq2_len)
    num_lines = 0
    for line in open(filename):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        if line.startswith(">"):
            continue
        
        i, j, length = map(int, line.strip().split())
        
        if seq1_begin is not None and seq1_end is not None:
            if i < seq1_begin or i >= seq1_end:
                continue
        
        if seq2_begin is not None and seq2_end is not None:
            if j < seq2_begin or j >= seq2_end:
                continue
        
        
        x = (i - seq1_begin) * coord_multiplier
        y = (j - seq2_begin) * coord_multiplier
        l = length * coord_multiplier
        drawing.add(drawing.line((x, y), (x + l, y + l), 
                                 stroke = color, stroke_width = line_width))
        
        num_lines += 1
        if num_lines % 10000 == 0:
            print(f"added {num_lines} lines")
            
    print(f"added {num_lines} lines")
    return num_lines

if __name__ == "__main__":
    
    
    
    mems = []
    for line in open(filename):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        if line.startswith(">"):
            continue
        
        mems.append(tuple(map(int, line.strip().split())))
        
    endpoints1 = set(item[0] for item in mems)
    endpoints2 = set(item[1] for item in mems)
    
    print("unique endpoints:")
    print("\tseq1: {}".format(len(endpoints1)))
    print("\tseq2: {}".format(len(endpoints2)))