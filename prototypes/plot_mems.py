#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:18:22 2022

@author: Jordan
"""

import sys
import svgwrite as svg
import matplotlib.pyplot as plt
import heapq
import rbtree
import rangetree
import math
import collections

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

def length_score_heuristic(seq_len1, seq_len2, matches):
    
    lengths = [l for x, y, l in matches]
    lengths.sort(reverse = True)
    i = 0
    s = 0
    while i < len(lengths) and s < min(seq_len1, seq_len2):
        s += lengths[i]
        i += 1
        
    approx_avg_len = s / i
    
    alpha = math.log(2 * seq_len1 * seq_len2 / 3) / math.log(approx_avg_len) - 2.0
    return alpha
    


def split_matches(matches):
    
    xtree = rbtree.RedBlackTree()
    ytree = rbtree.RedBlackTree()
    
    
    for x, y, l in matches:
        
        for v in (x, x + l):
            if xtree.find(v) is None:
                xtree.insert(v, None)
        for v in (y, y + l):
            if ytree.find(v) is None:
                ytree.insert(v, None)
        
    
    split = []
    
    i = 0
    
    for x, y, l in matches:
        
        i += 1
        if i % 50000 == 0:
            print("split {} matches".format(i))
        
        
        xnode = xtree.successor(xtree.find(x))
        ynode = ytree.successor(ytree.find(y))
        
        cursor = 0
        
        while xnode is not None or ynode is not None:
            
            if xnode is None:
                if ynode.key < y + l:
                    
                    split_len =  ynode.key - (y + cursor)
                    split.append((x + cursor, y + cursor, split_len))
                    cursor += split_len
                    
                    ynode = ytree.successor(ynode)
                else:
                    break
            elif ynode is None:
                if xnode.key < x + l:
                    
                    split_len =  xnode.key - (x + cursor)
                    split.append((x + cursor, y + cursor, split_len))
                    cursor += split_len
                    
                    xnode = xtree.successor(xnode)
                else:
                    break
            else:
                if xnode.key == ynode.key and xnode.key < x + l:
                    
                    split_len =  ynode.key - (y + cursor)
                    split.append((x + cursor, y + cursor, split_len))
                    cursor += split_len
                    
                    xnode = xtree.successor(xnode)
                    ynode = ytree.successor(ynode)
                elif xnode.key < ynode.key and xnode.key < x + l:
                    
                    split_len =  xnode.key - (x + cursor)
                    split.append((x + cursor, y + cursor, split_len))
                    cursor += split_len
                    
                    xnode = xtree.successor(xnode)
                elif ynode.key < xnode.key and ynode.key < y + l:
                    
                    split_len =  ynode.key - (y + cursor)
                    split.append((x + cursor, y + cursor, split_len))
                    cursor += split_len
                    
                    ynode = ytree.successor(ynode)
                else:
                    break
        
        split.append((x + cursor, y + cursor, l - cursor))
        
    return split
        
 
    


class MatchNode:
    
    def __init__(self, b1, b2, length, node_id):
        self.begin1 = b1
        self.begin2 = b2
        self.length = length
        self.node_id = node_id
        self.edges = []
        
        
def construct_graph_simple(matches):
    
    nodes = []
    
    for i in range(len(matches)):
        v, h, length = matches[i]
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

def compute_edge_score(node, next_node, match, gap_opens, gap_extends):
    
    diag_offset = abs((node.begin1 - node.begin2) - (next_node.begin1 - next_node.begin2))
    
    edge_score = 0
    if diag_offset != 0:
        edge_score = -min(diag_offset * gap_extend + gap_open for gap_extend, gap_open in zip(gap_extends, gap_opens))
        
    overlap = max(node.begin1 + node.length - next_node.begin1,
                  node.begin2 + node.length - next_node.begin2)
    
    if overlap > 0:
        edge_score -= overlap * match
    
    return edge_score


def graph_dp(graph, match, gap_opens, gap_extends):
    
    fwd = [-inf for i in range(len(graph))]
    bwd = [-inf for i in range(len(graph))]
    
    is_source = [True for i in range(len(graph))]
    
    for i in range(len(graph)):
        node = graph[i]
        if len(node.edges) == 0:
            bwd[i] = 0
        for j in node.edges:
            is_source[j] = False
            
    
    for i in range(len(graph)):
        node = graph[i]
        if is_source[i]:
            fwd[i] = match * node.length
            
        for j in node.edges:
            next_node = graph[j]
            
            edge_score = compute_edge_score(node, next_node, match, gap_opens, gap_extends)
                
#            if j == 51 and fwd[j] < fwd[i] + edge_score + match * next_node.length:
#                print("F: update 51 from {} to {}, fwd {}, edge {}, match {}".format(i, fwd[i] + edge_score + match * next_node.length, fwd[i], edge_score, match * next_node.length))
            fwd[j] = max(fwd[j], fwd[i] + edge_score + match * next_node.length)

    for i in range(len(graph) - 1, -1, -1):
        node = graph[i]
        for j in node.edges:
            next_node = graph[j]
            
            edge_score = compute_edge_score(node, next_node, match, gap_opens, gap_extends)
                
#            if i == 51 and bwd[i] < bwd[j] + edge_score + match * next_node.length:
#                print("B: update 51 from {} to {}, bwd {}, edge {}, match {}".format(j, bwd[j] + edge_score + match * next_node.length, bwd[j], edge_score, match * next_node.length))
                
            #print("i {}, j {}, bwd[i] {}, score {}".format(i, j, bwd[i], bwd[j] + edge_score + match * next_node.length))
            bwd[i] = max(bwd[i], bwd[j] + edge_score + match * next_node.length)
            
    
    return fwd, bwd
    

def simple_sweep_line_chain(matches, range_score = False):
    
    order = sorted(range(len(matches)), key = lambda i: matches[i])
    
    sweep_line = rbtree.RedBlackTree()
    
    end_heap = []
    
    backpointer = [None for i in range(len(matches))]
    match_score = [None for i in range(len(matches))]
    
    print_num = 0
    
    range_tree_begins = None
    range_tree_ends = None
    
    if range_score:
        print("making range tree for beginning points")
        range_tree_begins = rangetree.RangeTree2D(matches, True)
        print("making range tree for ending points")
        range_tree_ends = rangetree.RangeTree2D([(x + l, y + l, l) for x, y, l in matches], True)
    
    def dump_tree():
        sweep_line.prettyPrint()
        if sweep_line.root is not None:
            n = sweep_line.minimum(sweep_line.root)
            s = "scores:"
            while n is not None:
                s += " {} {},".format(n.value, match_score[n.value])
                n = sweep_line.successor(n)
            print(s)
    
    for i in order:
        
        if print_num % 50000 == 0:
            print("iter {} of {}, heap size {}, sweep line size {}, height {}".format(print_num, len(order), len(end_heap), sweep_line.size(), sweep_line.height()))
        
        y, x, length = matches[i]
        
        if range_score:
            score = (length * length * length
                    + range_tree_ends.range_weight_query(0, y, 0, x)
                    + range_tree_begins.range_weight_query(y + length, inf, x + length, inf))
        else:
            score = length
        
        while len(end_heap) != 0 and end_heap[0][0] <= y:
            
            end_y, end_x, idx = heapq.heappop(end_heap)
            
            node_prev = sweep_line.lower_bound(end_x)
            
            
            if node_prev is not None:
                if match_score[node_prev.value] >= match_score[idx]:
                    # automatically dominated, don't need it 
                    continue
                elif node_prev.key == end_x and match_score[idx] > match_score[node_prev.value]:
                    # we'll be replacing the previous node
                    sweep_line.delete(node_prev.key)
            
            node = sweep_line.insert(end_x, idx)
            # TODO: this is only O(1) amortized under a particular access pattern that 
            # i don't guarantee, could use level links for O(1) guaranteed
            succ = sweep_line.successor(node)
            while succ is not None and match_score[succ.value] <= match_score[idx]:
                tmp = sweep_line.successor(succ)
                sweep_line.delete(succ.key)
                succ = tmp
                
        
        node_prev = sweep_line.lower_bound(x)
        if node_prev is not None:
            backpointer[i] = node_prev.value
            dp_score = match_score[node_prev.value] + score
        else:
            dp_score = score
        
        match_score[i] = dp_score
        
        heapq.heappush(end_heap, (y + length, x + length, i))
        
        print_num += 1
        
    i = max(range(len(matches)), key = lambda j: match_score[j])
    
    print("opt score occurs at {} with score {}".format(i, match_score[i]))
    
    opt_anchors = []
    while i is not None:
#        print("backtrace to {}, y {}, x {}, l {}, score {}".format(i, matches[i][0], matches[i][1], matches[i][2], match_score[i]))
        opt_anchors.append(i)
        i = backpointer[i]
        
    opt_anchors.reverse()
        
    anchor_graph = []
    for i in range(len(opt_anchors)):
        idx = opt_anchors[i]
        v, h, l = matches[idx]
        anchor_graph.append(MatchNode(v, h, l, idx))
        if i + 1 != len(opt_anchors):
            anchor_graph[i].edges.append(i + 1)

    return anchor_graph

         
                
                
    
    

def select_near_optimal_matches(matches, match, gap_opens, gap_extends, seq1_len, seq2_len, opt_frac):
    
    
    graph = construct_graph_simple(matches)
    
    print("constructed graph of size {} with {} edges".format(len(graph), sum(len(n.edges) for n in graph)))
    
    fwd, bwd = graph_dp(graph, match, gap_opens, gap_extends)
    
#    print("fwd")
#    for f in fwd:
#        print("\t{}".format(f))
#    print("bwd")
#    for f in bwd:
#        print("\t{}".format(f))
    
    opt = -inf
    
    for i in range(len(matches)):
        if len(graph[i].edges) == 0:
            opt = max(opt, fwd[i])
    
    diff = int(opt_frac * match * (seq1_len + seq2_len) / 2)
    min_score = opt - diff
            
    print("opt score is {}, accepting down to score {}".format(opt, min_score))
    
#    printed = 0
    
    rejected_before = [0]
    selected = []
    for i in range(len(graph)):
        
        node = graph[i]
        
        if fwd[i] + bwd[i] >= min_score:
            selected.append(graph[i])
            rejected_before.append(rejected_before[-1])
            
            edges_rejected = 0
            
            for j in range(len(selected[-1].edges)):
                next_idx = selected[-1].edges[j]
                next_node = graph[next_idx]
                
                edge_score = compute_edge_score(node, next_node, match, gap_opens, gap_extends)
                
#                if printed < 1000:
#                    print("edge {} to {}, fwd {}, es {}, m {}, bwd {}, tot {}, keep? {}".format(i, j, fwd[i], edge_score, match * next_node.length, bwd[next_idx], fwd[i] + edge_score + match * next_node.length + bwd[next_idx], fwd[i] + edge_score + match * next_node.length + bwd[next_idx] >= min_score))
#                    printed += 1
                if fwd[i] + edge_score + match * next_node.length + bwd[next_idx] >= min_score:
                    if edges_rejected != 0:
                        selected[-1].edges[j - edges_rejected] = next_idx
                else:
                    edges_rejected += 1
                    
            for j in range(edges_rejected):
                selected[-1].edges.pop()
                
        else:
            rejected_before.append(rejected_before[-1] + 1)
            
            
    for node in selected:
        for i in range(len(node.edges)):
            node.edges[i] -= rejected_before[node.edges[i]]
            
            
    print("filtered to graph of size {} with {} edges".format(len(selected), sum(len(n.edges) for n in selected)))
            
    return selected

def graph_to_matches(graph):
    return [(n.begin1, n.begin2, n.length) for n in graph]


def flush_rectangle(expanded, point_tree, midpoint_to_match, min_node, max_node, 
                    prev_ending, ending, radius):
    
    d_begin = min_node.key - radius
    d_end = max_node.key + radius + 1 # past-the-last
    
    # query the midpoints
    midpoints = point_tree.range_query(d_begin, d_end, prev_ending, ending)
    
    for midpoint in midpoints:
        
        # record the match
        expanded.append(midpoint_to_match[midpoint])
    
def process_ending(expanded, point_tree, midpoint_to_match, diag_num_active, active,
                   prev_ending, endings, radius):
    
    ending, end_diag = heapq.heappop(endings)
            
    diag_num_active[end_diag] -= 1
    
    if diag_num_active[end_diag] == 0:
        # this diagonal is leaving the active set
        
        min_node = active.get_minimum()
        max_node = active.get_maximum()
        
        if min_node is not None and (end_diag == min_node.key or end_diag == max_node.key):
            # the outer boundary will shrink after this, so we need a rectangle
            
            flush_rectangle(expanded, point_tree, midpoint_to_match, min_node, max_node, 
                            prev_ending, ending, radius)
            
            
            # the next rectangle will start here
            prev_ending = ending
        
        # record that it is no longer active
        active.delete(end_diag)
    
    return prev_ending

def expand_graph(graph, matches, radius):
    
    print("expanding match graph of size {} with radius {}".format(len(graph), radius))
    
    # represent each match by its midpoint in (diagonal, antidiagonal) coordinates
    midpoint_to_match = {}
    ad_coord_matches = []
    for k in range(len(matches)):
        
        i, j, l = matches[k]
        
        mid_d = i - j
        mid_a = i + j + l
        
        midpoint_to_match[(mid_d, mid_a)] = k
        
        # TODO: silly to have the weight here, but it's what the interface expects
        ad_coord_matches.append((mid_d, mid_a, l))
        
    print("constructing midpoint range tree")
    # index the coordinates for range queries
    point_tree = rangetree.RangeTree2D(ad_coord_matches, True)

    print("extending matches using edges")
    
    # extend the matches forward and backward using their edges
    rev_edges = [[] for i in range(len(graph))]
    for i in range(len(graph)):
        for j in graph[i].edges:
            rev_edges[j].append(i)
            
    extended_matches = []
    for i in range(len(graph)):
        
        node = graph[i]
        
        d = node.begin1 - node.begin2
        a_min = node.begin1 + node.begin2
        a_max = a_min + 2 * node.length
        
        for j in rev_edges[i]:
            node_prev = graph[j]
            a_min = min(a_min, node_prev.begin1 + node_prev.begin2 + 2 * node_prev.length)
        for j in node.edges:
            node_next = graph[j]
            a_max = max(a_max, node_next.begin1 + node_next.begin2)
            
        extended_matches.append((a_min, a_max, d))
    
    extended_matches.sort()
    
    # the expanded set of matches
    expanded = []
    
    # state tracking structures for the sweep line
    active = rbtree.RedBlackTree()
    diag_num_active = collections.defaultdict(int)
    endings = []
    prev_ending = 0
    
    it = 0
    for a_min, a_max, d in extended_matches:
        
        it += 1
        if it % 1000 == 0:
            print("match iteration {} of {}".format(it, len(extended_matches)))
        
        while len(endings) != 0 and endings[0][0] <= a_min:
            
            # the next event is the end of a match
            
            prev_ending = process_ending(expanded, point_tree, midpoint_to_match, diag_num_active, active,
                                         prev_ending, endings, radius)
        
        if diag_num_active[d] == 0:
            # this diag is not currently in the state trackers
            
            min_node = active.get_minimum()
            max_node = active.get_maximum()
            
            if min_node is not None and (d < min_node.key or d > max_node.key):
                # we're expanding the rectangle's boundary, so finish the existing rectangle
                
                flush_rectangle(expanded, point_tree, midpoint_to_match, min_node, max_node, 
                                prev_ending, a_min, radius)
                
                prev_ending = a_min
            
            active.insert(d)
            
        # make updates for the ending
        heapq.heappush(endings, (a_max, d))
        diag_num_active[d] += 1
        
    # clear out the remaining endings
    while len(endings) != 0:
        prev_ending = process_ending(expanded, point_tree, midpoint_to_match, diag_num_active, active,
                                     prev_ending, endings, radius)
    
    return [matches[i] for i in expanded]


def contruct_colinearity_graph(matches, seq_len1, seq_len2):
    
    print("construct colinearity graph")
    
    beginning_to_match = {}
    ending_to_match = {}
    endings = []
    for i, j, l in matches:
        beginning_to_match[(i, j)] = l
        ending_to_match[(i + l, j + l)] = l
        
        endings.append((i + l, j + l, l))
        
#    print("creating range trees")
#    
#    beginning_tree = rangetree.RangeTree2D(matches, True)
#    ending_tree = rangetree.RangeTree2D(endings, True)
        
    
    
    # get a set of long matches that 
    size_order = sorted(range(len(matches)), key = lambda i: matches[i][2], reverse = True)
    total_len = 0
    for i in range(len(size_order)):
        total_len += matches[size_order[i]][2]
        if total_len > seq_len1 + seq_len2:
            size_order = size_order[:i+1]
            break



def load_matches(filename):
    matches = []
    for line in open(filename):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        if line.startswith(">"):
            continue
        
        matches.append(tuple(map(int, line.strip().split())))
    return matches
        

def filter_matches(matches, seq1_begin, seq1_end, seq2_begin, seq2_end, min_length = None):
    
    filtered = []
    for i, j, length in matches:
        if seq1_begin is not None and seq1_end is not None:
            if i < seq1_begin or i >= seq1_end:
                continue
        
        if seq2_begin is not None and seq2_end is not None:
            if j < seq2_begin or j >= seq2_end:
                continue
            
        if min_length is not None:
            if length < min_length:
                continue
        
        filtered.append((i, j, length))
    
    return filtered
    


def plot_lines(drawing, matches, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
               line_width, color, alpha = 1.0):
    coord_multiplier = float(long_side) / max(seq1_end - seq1_begin, seq2_end - seq2_begin)
    num_lines = 0
    for i, j, length in matches:
        
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


def antidiagonal_coverage(matches, seq1_len, seq2_len):
    
    events = []
    for i, j, l in matches:
        events.append((i + j, True))
        events.append((i + j + 2 * l, False))
    events.sort()
        
    curr_cov = 0
    cov_profile = []
    
    for a, is_begin in events:
        while len(cov_profile) < a:
            cov_profile.append(curr_cov)
            
        if is_begin:
            curr_cov += 1
        else:
            curr_cov -= 1
    
    while len(cov_profile) < seq1_len + seq2_len:
        cov_profile.append(0)
        
    
    return cov_profile


if __name__ == "__main__":
    
#    tree = SBBST()
#    tree.insert(1)
#    tree.insert(3)
#    tree.insert(5)
#    for i in range(7):
#        #print(tree.search(tree.head, i))
#        node = (bst_upper_bound(tree, i))
#        if node:
#            print(node.val)
#        else:
#            print("None")
#        
#    quit()
    
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    mems_file = sys.argv[3]
    mums_file = sys.argv[4]
    out = sys.argv[5]
    
    long_side = 3000
    line_width = 5
    cluster_line_width = 5
    cluster_line_alpha = 0.75
    
    seq1_begin = 1000000
    seq1_end = 2250000
    seq2_begin = 1000000
    seq2_end = 2250000
    
    min_length = 513
    
    print("loading sequences")
    
    seq1 = parse_fasta(fasta1)
    seq2 = parse_fasta(fasta2)
    
    print("loading matches")
    mems = load_matches(mems_file)
    mums = load_matches(mums_file)
        
    seq1_end = min(seq1_end, len(seq1))
    seq2_end = min(seq2_end, len(seq2))
    
    seq1_len = seq1_end - seq2_begin
    seq2_len = seq2_end - seq2_begin
    
    print(f"plotting in intervals {seq1_begin}:{seq1_end} and {seq2_begin}:{seq2_end}")
    
    filtered_mems = filter_matches(mems, seq1_begin, seq1_end, seq2_begin, seq2_end, min_length)
    filtered_mums = filter_matches(mums, seq1_begin, seq1_end, seq2_begin, seq2_end, min_length)
    
#    print("applying splits to matches")
#    split_mems = split_matches(mems)
#    print("splitting increases number of matches from {} to {}".format(len(mems), len(split_mems)))
    
    print("finding anchor set")
    
    match = 10
    gap_opens = [5, 100000]
    gap_extends = [10, 2]
    opt_frac = .025
    
    # attempt 1: clustering MUMs
#    cluster = select_near_optimal_matches(mums, match, gap_opens, gap_extends, seq1_len, seq2_len, opt_frac)
    
    # attempt 2: cluster MEMs that are long enough to nearly guarantee uniqueness
#    length_filtered_matches = filter_matches(mems, None, None, None, None, 1500)
#    cluster_graph = select_near_optimal_matches(length_filtered_matches, match, gap_opens, gap_extends, seq1_len, seq2_len, opt_frac)
#    
#    cluster = graph_to_matches(cluster_graph)
    
    # attempt 3: simple maximum colinear chain
#    alpha = length_score_heuristic(seq1_len, seq2_len, mems)
#    print("estimate range scoring exponent as {}".format(alpha))
    range_scoring = False
    cluster_graph = simple_sweep_line_chain(mems, range_scoring)
    #cluster = graph_to_matches(cluster_graph)
    
    expansion_radius = 10000
    cluster = expand_graph(cluster_graph, mems, expansion_radius)
    
    
    # filter for plotting
    filtered_cluster = filter_matches(cluster, seq1_begin, seq1_end, seq2_begin, seq2_end, min_length)
    
    
#    cov_profile = antidiagonal_coverage(cluster, seq1_len, seq2_len)
#    filtered_cov_profile = cov_profile[seq1_begin + seq2_begin:seq1_end + seq2_end]
#    height = 5
#    width = 50
#    f, ax = plt.subplots(1, 1, figsize = (width, height))
#    ax.tick_params(axis = 'x', labelsize = int(height * 2))
#    ax.plot(list(range(len(filtered_cov_profile))), filtered_cov_profile, '-', marker = None, linewidth = 1)
#    f.savefig("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/antidiagonal_coverage.svg", format = "svg", bbox_inches = "tight")
    
    multiplier1 = float(seq1_len) / max(seq1_len, seq2_len)
    multiplier2 = float(seq2_len) / max(seq1_len, seq2_len)
    
    height = long_side * multiplier1
    width = long_side * multiplier2
    
    print("making drawing of size {} x {}".format(height, width))
    
    drawing = svg.Drawing(out, size = (height, width))
    
    
    print("drawing MEMs")
    
    num_lines = 0
    num_lines += plot_lines(drawing, filtered_mems, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                            line_width, "black")
    
    print("drawing MUMs")
    
    num_lines += plot_lines(drawing, filtered_mums, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                            line_width, "red")
    
        
    print("drawing clustered matches")
    
    num_lines += plot_lines(drawing, filtered_cluster, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                            cluster_line_width, "lime", cluster_line_alpha)
    
    print(f"added {num_lines} lines in total")
    
    drawing.save()
    