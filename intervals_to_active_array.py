#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:45:44 2022

@author: Jordan
"""

import sys
import re

def construct_chain_graph(left_intervals, right_intervals, max_gap):
    
    #TODO: should i also have a max overlap?
    
    # even entries are forward, odd entries are reverse
    graph = [[] for i in range(2 * (len(left_intervals) + len(right_intervals)))]
    
    # add edges within the left side
    for i in range(len(left_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1 = left_intervals[i]
        for j in range(i + 1, len(left_intervals)):
            contig2, contig_b2, contig_e2, flank_b2, flank_e2 = left_intervals[j]
#            print("from left", left_intervals[i])
#            print("to left", left_intervals[j])
            if contig1 == contig2:
                # they are mapped to the same contig
                if (contig_e1 <= contig_b2 and flank_e1 <= flank_b2 and 
                    contig_b2 - contig_e1 <= max_gap and flank_b2 - flank_e1 <= max_gap):
                    # there is a forward edge
                    graph[2 * i].append(2 * j)
                    graph[2 * j + 1].append(2 * i + 1)
#                    print("add forward edges")
                if (contig_e2 <= contig_b1 and flank_e2 <= flank_b1 and
                    contig_b1 - contig_e2 <= max_gap and flank_b1 - flank_e2 <= max_gap):
                    # there is a reverse edge
                    graph[2 * j].append(2 * i)
                    graph[2 * i + 1].append(2 * j + 1)
#                    print("add reverse edge")
              
    # add edges between the left and right sides      
    for i in range(len(left_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1 = left_intervals[i]
        for j in range(len(right_intervals)):
#            print("from left", left_intervals[i])
#            print("to right", right_intervals[j])
            contig2, contig_b2, contig_e2, flank_b2, flank_e2 = right_intervals[j]
            if contig1 == contig2:
                if contig_e1 <= contig_b2:
                    graph[2 * i].append(2 * len(left_intervals) + 2 * j)
#                    print("add forward edge")
                if contig_e2 <= contig_b1:
                    graph[2 * i + 1].append(2 * len(left_intervals) + 2 * j + 1)
#                    print("add reverse edge")
            
    for i in range(len(right_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1 = right_intervals[i]
        # add edges within the right side
        for j in range(i + 1, len(right_intervals)):
            contig2, contig_b2, contig_e2, flank_b2, flank_e2 = right_intervals[j]
#            print("from right", right_intervals[i])
#            print("to right", right_intervals[j])
            if contig1 == contig2:
                # they are mapped to the same contig
                if (contig_e1 <= contig_b2 and flank_e1 <= flank_b2 and 
                    contig_b2 - contig_e1 <= max_gap and flank_b2 - flank_e1 <= max_gap):
                    # there is a forward edge
                    graph[2 * len(left_intervals) + 2 * i].append(2 * len(left_intervals) + 2 * j)
                    graph[2 * len(left_intervals) + 2 * j + 1].append(2 * len(left_intervals) + 2 * i + 1)
#                    print("add forward edge")
                if (contig_e2 <= contig_b1 and flank_e2 <= flank_b1 and
                    contig_b1 - contig_e2 <= max_gap and flank_b1 - flank_e2 <= max_gap):
                    # there is a reverse edge
                    graph[2 * len(left_intervals) + 2 * j].append(2 * len(left_intervals) + 2 * i)
                    graph[2 * len(left_intervals) + 2 * i + 1].append(2 * len(left_intervals) + 2 * j + 1)
#                    print("add reverse edge")
            
    return graph

def topological_order(graph):
    
    stack = []
    
    in_degree = [0 for i in range(len(graph))]
    
    for adj in graph:
        for j in adj:
            in_degree[j] += 1
            
    for i in range(len(graph)):
        if in_degree[i] == 0:
            stack.append(i)
        
    order = []
    while len(stack) != 0:
        i = stack.pop()
        order.append(i)
        for j in graph[i]:
            in_degree[j] -= 1
            if in_degree[j] == 0:
                stack.append(j)
                
    return order
    


def graph_dp(graph, left_intervals, right_intervals):
    
    dp = [-sys.maxsize for i in range(len(graph))]
    backpointer = [None for i in range(len(graph))]
    
    # must start on a left forward or right reverse interval
    for i in range(len(left_intervals)):
        dp[2 * i] = 0
        dp[2 * i + 1] = 0
        
    # do dynamic programming
    for i in topological_order(graph):
        if i < 2 * len(left_intervals):
            contig, contig_b, contig_e, flank_b, flank_e = left_intervals[i // 2]
        else:
            contig, contig_b, contig_e, flank_b, flank_e = right_intervals[(i - 2 * len(left_intervals)) // 2]
        
        total_len = contig_e - contig_b + flank_e - flank_b
        
#        print("DP at", i, "with score", dp[i], "extending to", dp[i] + total_len)
        for j in graph[i]:
            if dp[i] + total_len > dp[j]:
#                print("\textend to", j)
                dp[j] = dp[i] + total_len
                backpointer[j] = i
                
    opt_value = 0
    opt_index = None
    # check for opt on forward right intervals
    for i in range(len(right_intervals)):
        contig, contig_b, contig_e, flank_b, flank_e = right_intervals[i]
        interval_value = contig_e - contig_b + flank_e - flank_b
#        print("forward", i, 2 * len(left_intervals) + 2 * i)
        for j in range(2 * len(left_intervals) + 2 * i, 2 * len(left_intervals) + 2 * (i + 1)):
            value = dp[j] + interval_value
#            print("right", i, "->", 2 * len(left_intervals) + 2 * i, "value", value)
            if value > opt_value:
                opt_value = value
                opt_index = j

            
#    print(left_intervals)
#    print(right_intervals)
#    print(graph)
#    print(dp)
#    print(backpointer)
#    print(opt_index)
    
    if opt_index is None:
        return None, None, None
    
    left_traceback = []
    right_traceback = []
    is_reverse = (opt_index % 2 == 1)
    
    i = opt_index
    while i is not None:
        if i < 2 * len(left_intervals):
            left_traceback.append(left_intervals[i // 2])
        else:
            right_traceback.append(right_intervals[(i - 2 * len(left_intervals)) // 2])
        i = backpointer[i]
        
    left_traceback.reverse()
    right_traceback.reverse()
    
    return left_traceback, right_traceback, is_reverse
    

def align_stats(selection, flank_len):
    
    contig_aligned_begin = sys.maxsize
    contig_aligned_end = -sys.maxsize
    
    contig_aligned_len = 0
    flank_aligned_len = 0
    
    for contig, contig_b, contig_e, flank_b, flank_e in selection:
        contig_aligned_begin = min(contig_aligned_begin, contig_b)
        contig_aligned_end = max(contig_aligned_end, contig_e)
        
        contig_aligned_len += contig_e - contig_b
        flank_aligned_len += flank_e - flank_b
        
    contig_len = contig_aligned_end - contig_aligned_begin
    return contig_aligned_len / contig_len, flank_aligned_len / flank_len


if __name__ == "__main__":
    
    intervals = sys.argv[1]
    
    flank_regex = "\S+:(\d+)\-(\d+)\((\S+):(\d+)\-(\d+)([LR])\)"
    
    centromere_flank_alns = {}
    flank_lens = {}
    
    for line in open(intervals):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        # these are 0-based, half-open
        tokens = line.strip().split()
        
        target_contig = tokens[0]
        target_begin = int(tokens[1])
        target_end = int(tokens[2])
        flank = tokens[3]
        flank_aligned_begin = int(tokens[4])
        flank_aligned_end = int(tokens[5])
        
        m = re.match(flank_regex, flank)
        
        assert(m is not None)
        
        # these are 1-based, closed
        flank_begin = int(m.group(1)) - 1
        flank_end = int(m.group(2))
                
        centro_chr = m.group(3)
        centro_begin = int(m.group(4)) - 1
        centro_end = int(m.group(5))
        side = 0 if m.group(6) == "L" else 1
        
        
        centro = (centro_chr, centro_begin, centro_end)
        
        flank_lens[(centro, side)] = flank_end - flank_begin
        
        if centro not in centromere_flank_alns:
            centromere_flank_alns[centro] = [[], []]
            
        centromere_flank_alns[centro][side].append((target_contig, target_begin, target_end, flank_aligned_begin, flank_aligned_end))
    
    # header
    columns = ["centro_chrom", "centro_begin", "centro_end", "flank_contig", 
               "left_flank_begin", "left_flank_end", "right_flank_begin", "right_flank_end",
               "left_contig_begin", "left_contig_end", "right_contig_begin", "right_contig_end",
               "left_contig_coverage", "right_contig_coverage",
               "left_flank_coverage", "right_flank_coverage", "is_reverse"]
    print("\t".join(columns))
    
    for centro in sorted(centromere_flank_alns):
        for side in centromere_flank_alns[centro]:
            side.sort()
            
        centro_chrom, centro_begin, centro_end = centro
    
        flank_len_left = flank_lens[(centro, 0)]
        flank_len_right = flank_lens[(centro, 1)]
        
        # TODO: separate max gaps for each side?
        max_gap = max(flank_len_left, flank_len_right) // 2
    
        left_intervals, right_intervals = centromere_flank_alns[centro]
        chain_graph = construct_chain_graph(left_intervals, right_intervals, max_gap)
        
        left_selection, right_selection, is_reverse = graph_dp(chain_graph, left_intervals, right_intervals)
        
        if left_selection is None:
            print("\t".join([str(v) for v in [centro_chrom, centro_begin, centro_end]] + (len(columns) - 3) * ["NA"]))
            continue
        
#        print(centro, "rev", is_reverse)
#        print("\tleft")
#        for v in left_selection:
#            print("\t{}".format(v))
#        print("\tright")
#        for v in right_selection:
#            print("\t{}".format(v))
        
        left_contig_coverage, left_flank_coverage = align_stats(left_selection, flank_len_left)
        right_contig_coverage, right_flank_coverage = align_stats(right_selection, flank_len_right)
        
        flank_contig = left_selection[0][0]
        
        left_flank_begin = min(rec[1] for rec in left_selection)
        left_flank_end = max(rec[2] for rec in left_selection)
        right_flank_begin = min(rec[1] for rec in right_selection)
        right_flank_end = max(rec[2] for rec in right_selection)
        
        left_contig_begin = min(rec[3] for rec in left_selection)
        left_contig_end = max(rec[4] for rec in left_selection)
        right_contig_begin = min(rec[3] for rec in right_selection)
        right_contig_end = max(rec[4] for rec in right_selection)
        
        
        print("\t".join(str(v) for v in [centro_chrom, centro_begin, centro_end, flank_contig,
                                         left_flank_begin, left_flank_end, right_flank_begin, right_flank_end,
                                         left_contig_begin, left_contig_end, right_contig_begin, right_contig_end,
                                         left_contig_coverage, right_contig_coverage,
                                         left_flank_coverage, right_flank_coverage, int(is_reverse)]))
        
        
        
        
        
        
        
        