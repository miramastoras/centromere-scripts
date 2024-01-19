#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 14:06:11 2023

@author: Jordan
"""

import subprocess
import sys
import re
import os
import pandas as pd
import numpy as np
import tempfile
import collections

import argparse

tmp_files = []

def tmp_file_name(suffix = None):
    tmp = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    if suffix is not None:
        tmp += "." + suffix
    tmp_files.append(tmp)
    return tmp

def get_contig_lens(fasta):
    contig_lens = {}
    contig = None
    for line in open(fasta):
        if line.startswith(">"):
            contig = line[1:].split()[0]
            contig_lens[contig] = 0
        else:
            contig_lens[contig] += len(line.strip())
    return contig_lens

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def aligned_interval(cigar):
    
    i = 0
    
    read_start = 0
    if len(cigar) > 0 and cigar[0][0] in ["H", "S"]:
        read_start += cigar[0][1]
        i += 1
        
    read_end = read_start
    aligned_length = 0
    
    while i < len(cigar):
        op, l = cigar[i]
        
        if op in ["I", "X", "=", "M"]:
            read_end += l
        if op in ["D", "X", "=", "M"]:
            aligned_length += l
        
        i += 1
        
    return read_start, read_end, aligned_length
       
def num_matches(cigar):
    m = 0
    for op, l in cigar:
        if op in ["M", "="]:
            m += l
    return m


def get_paf(sam):
    
    records = []
    
    contig_lens = {}
    
    for line in open(sam):
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith("@"):
            if line.startswith("@SQ"):
                ms = re.search("SN:(\S+)", line)
                ml = re.search("LN:(\d+)", line)
                contig_lens[ms.group(1)] = int(ml.group(1))
            continue
        
        tokens = line.strip().split("\t")
        
        name = tokens[0]
        flag = int(tokens[1])
        contig = tokens[2]
        ref_start = int(tokens[3]) - 1
        mapq = int(tokens[4])
        cigar = parse_cigar(tokens[5])
        seq_len = len(tokens[9])
        if contig == '*':
            contig_len = 0
        else:
            contig_len = contig_lens[contig]
        
        if flag & 16:
            strand = '-'
        else:
            strand = '+'
        
        
        read_start, read_end, aligned_length = aligned_interval(cigar)
        ref_end = ref_start + aligned_length
        matches = num_matches(cigar)
        block_length = sum(c[1] for c in cigar)
        
        #               0     1        2           3         4       5       6           7          8        9        10            11 
        records.append((name, seq_len, read_start, read_end, strand, contig, contig_len, ref_start, ref_end, matches, block_length, mapq))
        
    records.sort(key = lambda r : (r[5], r[7], r[8], r[0], r[2], r[3], r))
        
    return records


def construct_chain_graph(left_intervals, right_intervals, max_gap):
    
    #TODO: should i also have a max overlap?
    
    # even entries are forward, odd entries are reverse
    graph = [[] for i in range(2 * (len(left_intervals) + len(right_intervals)))]
    
    # add edges within the left side
    for i in range(len(left_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1, rev1 = left_intervals[i]
        for j in range(i + 1, len(left_intervals)):
            contig2, contig_b2, contig_e2, flank_b2, flank_e2, rev2 = left_intervals[j]
            if rev1 != rev2 or contig1 != contig2:
                continue
#            print("from left", left_intervals[i])
#            print("to left", left_intervals[j])
            if (contig_e1 <= contig_b2 and flank_e1 <= flank_b2 and 
                contig_b2 - contig_e1 <= max_gap and flank_b2 - flank_e1 <= max_gap):
                # there is a forward edge
                graph[2 * i].append(2 * j)
                graph[2 * j + 1].append(2 * i + 1)
#                print("add forward edges")
            if (contig_e2 <= contig_b1 and flank_e2 <= flank_b1 and
                contig_b1 - contig_e2 <= max_gap and flank_b1 - flank_e2 <= max_gap):
                # there is a reverse edge
                graph[2 * j].append(2 * i)
                graph[2 * i + 1].append(2 * j + 1)
#                print("add reverse edge")
              
    # add edges between the left and right sides      
    for i in range(len(left_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1, rev1 = left_intervals[i]
        for j in range(len(right_intervals)):
#            print("from left", left_intervals[i])
#            print("to right", right_intervals[j])
            contig2, contig_b2, contig_e2, flank_b2, flank_e2, rev2 = right_intervals[j]
            if rev1 != rev2 or contig1 != contig2:
                continue
            if contig_e1 <= contig_b2:
                graph[2 * i].append(2 * len(left_intervals) + 2 * j)
#                print("add forward edge")
            if contig_e2 <= contig_b1:
                graph[2 * i + 1].append(2 * len(left_intervals) + 2 * j + 1)
#                print("add reverse edge")
            
    for i in range(len(right_intervals)):
        contig1, contig_b1, contig_e1, flank_b1, flank_e1, rev1 = right_intervals[i]
        # add edges within the right side
        for j in range(i + 1, len(right_intervals)):
            contig2, contig_b2, contig_e2, flank_b2, flank_e2, rev1 = right_intervals[j]
            if rev1 != rev2 or contig1 != contig2:
                continue
#            print("from right", right_intervals[i])
#            print("to right", right_intervals[j])
            # they are mapped to the same contig
            if (contig_e1 <= contig_b2 and flank_e1 <= flank_b2 and 
                contig_b2 - contig_e1 <= max_gap and flank_b2 - flank_e1 <= max_gap):
                # there is a forward edge
                graph[2 * len(left_intervals) + 2 * i].append(2 * len(left_intervals) + 2 * j)
                graph[2 * len(left_intervals) + 2 * j + 1].append(2 * len(left_intervals) + 2 * i + 1)
#                print("add forward edge")
            if (contig_e2 <= contig_b1 and flank_e2 <= flank_b1 and
                contig_b1 - contig_e2 <= max_gap and flank_b1 - flank_e2 <= max_gap):
                # there is a reverse edge
                graph[2 * len(left_intervals) + 2 * j].append(2 * len(left_intervals) + 2 * i)
                graph[2 * len(left_intervals) + 2 * i + 1].append(2 * len(left_intervals) + 2 * j + 1)
#                print("add reverse edge")
            
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
            contig, contig_b, contig_e, flank_b, flank_e, rev = left_intervals[i // 2]
        else:
            contig, contig_b, contig_e, flank_b, flank_e, rev = right_intervals[(i - 2 * len(left_intervals)) // 2]
        
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
        contig, contig_b, contig_e, flank_b, flank_e, rev = right_intervals[i]
        interval_value = contig_e - contig_b + flank_e - flank_b
#        print("forward", i, 2 * len(left_intervals) + 2 * i)
        for j in range(2 * len(left_intervals) + 2 * i, 2 * len(left_intervals) + 2 * (i + 1)):
            value = dp[j] + interval_value
#            print("right", i, "->", 2 * len(left_intervals) + 2 * i, "value", value)
            if value > opt_value:
                opt_value = value
                opt_index = j

    
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
    
    for contig, contig_b, contig_e, flank_b, flank_e, rev in selection:
        contig_aligned_begin = min(contig_aligned_begin, contig_b)
        contig_aligned_end = max(contig_aligned_end, contig_e)
        
        contig_aligned_len += contig_e - contig_b
        flank_aligned_len += flank_e - flank_b
        
    contig_len = contig_aligned_end - contig_aligned_begin
    return contig_aligned_len / contig_len, flank_aligned_len / flank_len


def paf_to_active_array(paf, min_cov):
    
    intervals = paf
    
    flank_regex = "\S+:(\d+)\-(\d+)\((\S+)([\+\-]):(\d+)\-(\d+)([LR])\)"
    
    centromere_flank_alns = {}
    flank_lens = {}
    
    for tokens in paf:
        
        # if (tokens[5] == "HG002#2#JAHKSD010000046.1" or tokens[5] == "HG002#2#JAHKSD010000013.1"):
        #     print("\t".join(str(v) for v in tokens), file = sys.stderr)
        
        target_contig = tokens[5]
        
        target_begin = int(tokens[7])
        target_end = int(tokens[8])
        if tokens[4] == '+':
            rev = False
        else:
            rev = True
        flank = tokens[0]
        flank_aligned_begin = int(tokens[2])
        flank_aligned_end = int(tokens[3])
        
        m = re.match(flank_regex, flank)
        
        assert(m is not None)
        
        # these are 1-based, closed
        flank_begin = int(m.group(1)) - 1
        flank_end = int(m.group(2))
                
        centro_chr = m.group(3)
        centro_rev = m.group(4) == '-'
        centro_begin = int(m.group(5)) - 1
        centro_end = int(m.group(6))
        side = 0 if m.group(7) == "L" else 1
        
        
        centro = (centro_chr, centro_begin, centro_end, centro_rev)
        
        flank_lens[(centro, side)] = flank_end - flank_begin
        
        if centro not in centromere_flank_alns:
            centromere_flank_alns[centro] = [[], []]
            
        centromere_flank_alns[centro][side].append((target_contig, target_begin, target_end, flank_aligned_begin, flank_aligned_end, rev))
    
    # header
    columns = ["centro_chrom", "centro_begin", "centro_end", "centro_rev", "flank_contig", 
                "left_flank_begin", "left_flank_end", "right_flank_begin", "right_flank_end",
                "left_contig_begin", "left_contig_end", "right_contig_begin", "right_contig_end",
                "left_contig_coverage", "right_contig_coverage",
                "left_flank_coverage", "right_flank_coverage", "is_reverse"]
    
    
    active_arrays = []
    
    for centro in sorted(centromere_flank_alns):
        for side in centromere_flank_alns[centro]:
            side.sort()
            
        centro_chrom, centro_begin, centro_end, centro_rev = centro
    
        flank_len_left = flank_lens[(centro, 0)]
        flank_len_right = flank_lens[(centro, 1)]
        
        # TODO: separate max gaps for each side?
        max_gap = max(flank_len_left, flank_len_right) // 2
    
        left_intervals, right_intervals = centromere_flank_alns[centro]
        chain_graph = construct_chain_graph(left_intervals, right_intervals, max_gap)
        
        # if centro_chrom == "chrX":
        #     for i in range(len(chain_graph)):
        #         print("{}: {}".format(i, ", ".join(str(v) for v in chain_graph[i])), file = sys.stderr)
        
        
        left_selection, right_selection, is_reverse = graph_dp(chain_graph, left_intervals, right_intervals)
        
        if left_selection is None:
            # print("\t".join([str(v) for v in [centro_chrom, centro_begin, centro_end]] + (len(columns) - 3) * ["NA"]))
            continue
        
        
        left_contig_coverage, left_flank_coverage = align_stats(left_selection, flank_len_left)
        right_contig_coverage, right_flank_coverage = align_stats(right_selection, flank_len_right)
        
        if (left_contig_coverage < min_cov or left_flank_coverage < min_cov
            or right_contig_coverage < min_cov or right_flank_coverage < min_cov):
            continue
        
        flank_contig = left_selection[0][0]
        
        left_flank_begin = min(rec[3] for rec in left_selection)
        left_flank_end = max(rec[4] for rec in left_selection)
        right_flank_begin = min(rec[3] for rec in right_selection)
        right_flank_end = max(rec[4] for rec in right_selection)
        
        left_contig_begin = min(rec[1] for rec in left_selection)
        left_contig_end = max(rec[2] for rec in left_selection)
        right_contig_begin = min(rec[1] for rec in right_selection)
        right_contig_end = max(rec[2] for rec in right_selection)
        
        
        active_arrays.append([centro_chrom, centro_begin, centro_end, centro_rev, flank_contig,
                             left_flank_begin, left_flank_end, right_flank_begin, right_flank_end,
                             left_contig_begin, left_contig_end, right_contig_begin, right_contig_end,
                             left_contig_coverage, right_contig_coverage,
                             left_flank_coverage, right_flank_coverage, is_reverse])
    
    
    return pd.DataFrame(active_arrays, columns = columns)

def annotation_to_active_array(hum_as_bed, fasta, slop_distance, min_length, min_frac, min_flank, strand_scope,
                               mergeable_chroms):
    
    curr_contig = None
    
    curr_intervals = {}
    
    live_array_regex = "(S\d+C[\dXY/]+H\d+L)\.[\d/]+"
    
    chrom_regex = "C([\dXY/]+)H"
    
    columns = ["contig", "begin", "end", "chrom_set", "strand"]
    
    get_chrom_set = lambda hor: ",".join("chr" + v for v in re.search(chrom_regex, hor).group(1).split("/"))
    
    def determine_strand(left_revs, right_revs):
        num_rev_left = 0
        num_rev_right = 0
        for rev in left_revs:
            if rev:
                num_rev_left += 1
        for rev in right_revs:
            if rev:
                num_rev_right += 1
        
        frac_left = num_rev_left / len(left_revs)
        frac_right = num_rev_right / len(right_revs)
        if (frac_left > .5) != (frac_right > .5):
            return None
        else:
            return frac_left > .5
    
    table = []
    for line in open(hum_as_bed):
        tokens = line.strip().split()
        contig = tokens[0]
        begin = int(tokens[1])
        end = int(tokens[2])
        rev = (tokens[5] == "-")
        
        if contig != curr_contig:
            # flush the current intervals
            for hor in curr_intervals:
                curr_begin, curr_end, cov, left_revs, right_revs = curr_intervals[hor]
                if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
                    strand = determine_strand(left_revs, right_revs)
                    if strand is not None:
                        table.append([curr_contig, curr_begin, curr_end, get_chrom_set(hor), '-' if strand else '+'])
            curr_intervals.clear()
            curr_contig = contig
        
        
        m = re.match(live_array_regex, tokens[3])

        if m is None:
            continue
        
        hor = m.group(1)
        
        if hor in curr_intervals:
            curr_begin, curr_end, cov, left_revs, right_revs = curr_intervals[hor]
            if curr_end + slop_distance >= begin:
                curr_intervals[hor][1] = end
                curr_intervals[hor][2] += end - begin
                if len(left_revs) < strand_scope:
                    left_revs.append(rev)
                    right_revs.append(rev)
                elif len(right_revs) == strand_scope:
                    right_revs.popleft()
                    right_revs.append(rev)
            else:
                if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
                    strand = determine_strand(left_revs, right_revs)
                    if strand is not None:
                        table.append([curr_contig, curr_begin, curr_end, get_chrom_set(hor), '-' if strand else '+'])
                left_revs = collections.deque()
                right_revs = collections.deque()
                left_revs.append(rev)
                right_revs.append(rev)
                curr_intervals[hor] = [begin, end, end - begin, left_revs, right_revs]
        else:
            left_revs = collections.deque()
            right_revs = collections.deque()
            left_revs.append(rev)
            right_revs.append(rev)
            curr_intervals[hor] = [begin, end, end - begin, left_revs, right_revs]
        
        
    
    for hor in curr_intervals:
        curr_begin, curr_end, cov, left_revs, right_revs = curr_intervals[hor]
        if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
            strand = determine_strand(left_revs, right_revs)
            if strand is not None:
                table.append([curr_contig, curr_begin, curr_end, get_chrom_set(hor), '-' if strand else '+'])
    
    contig_lens = get_contig_lens(fasta)
    table = [row for row in table if row[1] >= min_flank and contig_lens[row[0]] - row[2] >= min_flank]        


    # merge the chromosome that we know have interrupted HOR arrays
    merged_table = []
    b = 0
    while b < len(table):
        e = b + 1
        while e < len(table) and table[e][0] == table[b][0]:
            e += 1
        
        to_remove = set()
        
        for i in range(b, e):
            if i in to_remove:
                continue
            for j in range(i + 1, e):
                if table[i][3] == table[j][3] and table[i][4] == table[j][4] and table[i][3] in mergeable_chroms:
                    # enlarge to include all of them
                    table[i][1] = min(table[i][1], table[j][1])
                    table[i][2] = max(table[i][2], table[j][2])
                    to_remove.add(j)
                    
        
        for i in range(b, e):
            if i not in to_remove:
                merged_table.append(table[i][:])
        
        b = e
        
    
    return pd.DataFrame(merged_table, columns = columns)

def get_inferred_array(flank_active_array_table, row):
    if flank_active_array_table.is_reverse.values[row]:
        return (flank_active_array_table.right_contig_end.values[row], 
                flank_active_array_table.left_contig_begin.values[row])
    else:
        return (flank_active_array_table.left_contig_end.values[row],
                flank_active_array_table.right_contig_begin.values[row])

def get_centro_strands(bed):
    
    strands = {}
    
    for line in open(bed):
        tokens = line.strip().split()
        assert(len(tokens) >= 6)
        assert(tokens[5] == '-' or tokens[5] == '+')
        rev = tokens[5] == '-'
        strands[tokens[0]] = rev
    
    return strands
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly", type=str, required = True,
                        help="FASTA to look for HOR arrays in")
    parser.add_argument("-f", "--flank_files", type=str, required = True,
                        help="comma separated list of reference flank FASTA files (from extract_flanks.py)")
    parser.add_argument("-b", "--hum_as_bed", type=str, required = True,
                        help="alpha satellite annotations from Hum-AS_HMMER")
    parser.add_argument("-H", "--hor_strand_bed", type=str, required = True,
                        help="BED file that gives centromere strand by chromosome")
    parser.add_argument("-m", "--min_flank_cov", type=float, default = 0.25,
                        help="minimum proportion of mapped flanks to identify an array from flank input")
    parser.add_argument("-d", "--merge_dist", type=int, default = 20000,
                        help="max distance to merge arrays identified from annotation input")
    parser.add_argument("-l", "--min_length", type=int, default = 25000,
                        help="minimum length to identify an array from annotation input")
    parser.add_argument("-c", "--min_anno_cov", type=float, default = 0.33,
                        help="minimum proportion live satellite to identify an array from annotation input")
    parser.add_argument("-s", "--min_anno_shoulder", type=float, default = 20000,
                        help="minimum distance from array to ends of contig to identify array from annotation input")
    parser.add_argument("-w", "--strand_window", type=int, default = 50,
                        help="determine strand of annotation-derived arrays with this many monomers on each end")
    parser.add_argument("-M", "--mergeable_chroms", type=str, default = "chr3,chr4",
                        help="combine annotation-derived arrays from the same contig for these chromosomes")
    
    args = parser.parse_args()
        
    flank_active_array_tables = []
    
    debug = False
    
    centro_strands = get_centro_strands(args.hor_strand_bed)
    
    if debug:
        print(centro_strands, file = sys.stderr)
    
    for flank_file in args.flank_files.split(","):
        
        print(f"mapping HOR flanks from {flank_file}", file = sys.stderr)
        
        # TODO: this is wasteful cruft, i could now run without the -a and directly use the paf
        # sam_file = "hg002.sam"
        sam_file = tmp_file_name("sam")
        subprocess.check_call(f"winnowmap -x asm20 -a {args.assembly} {flank_file} > {sam_file}", shell = True)

        print(f"identifying HOR arrays from mapped flanks", file = sys.stderr)
        paf = get_paf(sam_file)
        
        flank_active_array_tables.append(paf_to_active_array(paf, args.min_flank_cov))
    
    mergeable_chroms = set(args.mergeable_chroms.split(","))
    
    print(f"identifying HOR arrays from Hum-AS_HMMER annotations", file = sys.stderr)
    anno_active_arrays = annotation_to_active_array(args.hum_as_bed, args.assembly, args.merge_dist, args.min_length, 
                                                    args.min_anno_cov, args.min_anno_shoulder, args.strand_window,
                                                    mergeable_chroms)
    
    if debug:
        pd.set_option('display.max_columns', None)
        for i in range(len(flank_active_array_tables)):
            print("flank mapping table {}".format(i), file = sys.stderr)
            print(flank_active_array_tables[i], file = sys.stderr)
        print("annotation table {}".format(i), file = sys.stderr)
        print(anno_active_arrays, file = sys.stderr)
    
    print(f"integrating HOR evidence", file = sys.stderr)
    
    contig_to_table_row = collections.defaultdict(list)
    for i in range(len(flank_active_array_tables)):
        table = flank_active_array_tables[i]
        for j in range(len(table)):
            contig_to_table_row[table.flank_contig.values[j]].append((i, j))
    
    for j in range(len(anno_active_arrays)):
        contig_to_table_row[anno_active_arrays.contig.values[j]].append((len(flank_active_array_tables), j))


    # FIXME: this will fail on instances of broken-up HOR arrays like chr3 and chr4, where we could have more
    # than one annotation array on a contig
    num_inconsistent = 0
    num_incomplete = 0
    num_nonoverlapping = 0
    num_filtered = 0
    num_ambigous_strand = 0
    
    ambiguous_buffer = []
    
    chroms_found = set()
    
    chroms_filtered = set()
    
    for contig in sorted(contig_to_table_row):
        anno_row = -1
        flank_rows = []
        for idx, row in contig_to_table_row[contig]:
            if idx < len(flank_active_array_tables):
                flank_rows.append((idx, row))
            else:
                anno_row = row
        
        if anno_row != -1 and len(flank_rows) != 0:
            # we have both annotation and flank derived HOR arrays on this contig
            
            chrom_set = set(anno_active_arrays.chrom_set.values[anno_row].split(","))
            max_begin = anno_active_arrays.begin.values[anno_row]
            min_end = anno_active_arrays.end.values[anno_row]
            strand = anno_active_arrays.strand.values[anno_row]
            
            for flank_idx, flank_row in flank_rows:
                flank_set = set([ flank_active_array_tables[flank_idx].centro_chrom.values[flank_row] ])
                chrom_set.intersection_update(flank_set)
                begin, end = get_inferred_array(flank_active_array_tables[flank_idx], flank_row)
                max_begin = max(max_begin, begin)
                min_end = min(min_end, end)
                rev = (flank_active_array_tables[flank_idx].is_reverse.values[flank_row] != 
                       flank_active_array_tables[flank_idx].centro_rev.values[flank_row])
                flank_strand = '-' if rev else '+'
                if flank_strand != strand:
                    strand = None
            
            
            if len(chrom_set) == 1 and max_begin < min_end and strand is not None:
                # there is one unambiguous chromosome assignment and strand across all data types
                chrom = next(iter(chrom_set))
                chroms_found.add(chrom)
                #print("chrom {} array found on {} from integrated input".format(chrom, contig), file = sys.stderr)
                print("\t".join(str(v) for v in [contig, max_begin, min_end, chrom, 3, strand]))
            else:
                num_filtered += 1
                for chrom in chrom_set:
                    chroms_filtered.add(chrom)
                if len(chrom_set) > 1:
                    num_incomplete += 1
                if len(chrom_set) == 0:
                    num_inconsistent += 1
                if max_begin >= min_end:
                    num_nonoverlapping += 1
                if strand is None:
                    num_ambigous_strand += 1
                if debug:
                    print("filter flank-annotation contig {}, chroms {}, missing rev? {}, begin {}, end {} ".format(contig, chrom_set,
                                                                                                  strand is None, 
                                                                                                  max_begin, min_end), file = sys.stderr)
            
        elif len(flank_rows) != 0:
            # we have only flank derived HOR arrays on this contig
            flank_idx, flank_row = flank_rows[0]
            chrom = flank_active_array_tables[flank_idx].centro_chrom.values[flank_row]
            max_begin, min_end = get_inferred_array(flank_active_array_tables[flank_idx], flank_row)
            is_reverse = (flank_active_array_tables[flank_idx].is_reverse.values[flank_row] != 
                          flank_active_array_tables[flank_idx].centro_rev.values[flank_row])
            
            for i in range(1, len(flank_rows)):
                flank_idx, flank_row = flank_rows[i]
                if flank_active_array_tables[flank_idx].centro_chrom.values[flank_row] != chrom:
                    chrom = None
                begin, end = get_inferred_array(flank_active_array_tables[flank_idx], flank_row)
                max_begin = max(max_begin, begin)
                min_end = min(min_end, end)
                rev = (flank_active_array_tables[flank_idx].is_reverse.values[flank_row] != 
                       flank_active_array_tables[flank_idx].centro_rev.values[flank_row])
                if rev != is_reverse:
                    is_reverse = None
            
            if chrom is not None and is_reverse is not None and max_begin < min_end:
                # the flank mappings are consistent
                #print("chrom {} array found on {} from flank input".format(chrom, contig), file = sys.stderr)
                chroms_found.add(chrom)
                print("\t".join(str(v) for v in [contig, max_begin, min_end, chrom, 1, '-' if is_reverse else '+']))
            else:
                num_filtered += 1
                if chrom is None:
                    num_inconsistent += 1
                else:
                    chroms_filtered.add(chrom)
                if max_begin >= min_end:
                    num_nonoverlapping += 1
                if debug:
                    print("filter flank contig {}, inconsistent chrom? {}, missing rev? {}, begin {}, end {} ".format(contig,
                                                                                                               chrom is None,
                                                                                                               is_reverse is None,
                                                                                                               max_begin, min_end), file = sys.stderr)
            
        else:
            # we have only annotation derived HOR arrays on this contig
            chrom_set = anno_active_arrays.chrom_set.values[anno_row].split(",")
            if len(chrom_set) != 1:
                # this is a cross-chromosome family and we can't disambiguate it with the flanks
                ambiguous_buffer.append(contig)
                if debug:
                    print("buffering annotation contig {} with annotation chrom set {}".format(contig, chrom_set), file = sys.stderr)
                continue
            chrom = chrom_set[0]
            chroms_found.add(chrom)
            begin = anno_active_arrays.begin.values[anno_row]
            end = anno_active_arrays.end.values[anno_row]
            strand_rev = anno_active_arrays.strand.values[anno_row]
            strand = '+' if (strand_rev == centro_strands[chrom]) else '-'
            #print("chrom {} array found on {} from annotation input".format(chrom, contig), file = sys.stderr)
            print("\t".join(str(v) for v in [contig, begin, end, chrom, 2, strand]))
    
    # check if we ruled out all of the other chromosomes in an ambiguous assignment
    for contig in ambiguous_buffer:
        anno_row = -1
        for idx, row in contig_to_table_row[contig]:
            anno_row = row
            assert(idx == len(flank_active_array_tables))
            
        chrom_set = [c for c in anno_active_arrays.chrom_set.values[anno_row].split(",") if c not in chroms_found]
        if len(chrom_set) != 1:
            # this is a cross-chromosome family and we can't disambiguate it with the flanks
            for chrom in chrom_set:
                chroms_filtered.add(chrom)
            num_filtered += 1
            num_incomplete += 1
            if debug:
                print("filtering annotation contig {} with annotation chrom set {}".format(contig, chrom_set), file = sys.stderr)
            continue
        chrom = chrom_set[0]
        chroms_found.add(chrom)
        begin = anno_active_arrays.begin.values[anno_row]
        end = anno_active_arrays.end.values[anno_row]
        strand_rev = anno_active_arrays.strand.values[anno_row]
        strand = '+' if (strand_rev == centro_strands[chrom]) else '-'
        #print("chrom {} array found on {} from annotation input".format(chrom, contig), file = sys.stderr)
        print("\t".join(str(v) for v in [contig, begin, end, chrom, 2, strand]))
    
    print(f"filtered out arrays on {num_filtered} contigs: {num_inconsistent} inconsistent chromosome, {num_nonoverlapping} non-overlapping, {num_ambigous_strand} ambiguous strand, {num_incomplete} incomplete chromosome assignment", file = sys.stderr)
    if num_filtered != 0:
        print("filtering may have affected: {}".format(", ".join(sorted(chroms_filtered))), file = sys.stderr)

    for tmp_file in tmp_files:
        os.remove(tmp_file)