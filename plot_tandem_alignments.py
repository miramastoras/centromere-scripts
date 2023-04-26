#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 16:54:32 2022

@author: Jordan
"""

import sys
import svgwrite as svg
import tempfile
import subprocess
import re
import os
import pandas as pd
import numpy as np
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


def parse_matches(filename):
    matches = []
    for line in open(filename):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        if line.startswith(">"):
            continue
        
        matches.append(tuple(map(int, line.strip().split())))
    return matches
     

def get_matches(fasta1, fasta2, min_length = 300):
    
    sa_file = get_temp_file_name()
    mem_file = get_temp_file_name()
    mum_file = get_temp_file_name()
    
    print("getting MEMs", file = sys.stderr)
    subprocess.check_call(f"mummer -maxmatch -l {min_length} -L -q 8 -save {sa_file} {fasta1} {fasta2} > {mem_file} 2> /dev/null", shell = True)
    print("getting MUMs", file = sys.stderr)
    subprocess.check_call(f"mummer -mum -l {min_length} -L -q 8 -load {sa_file} {fasta1} {fasta2} > {mum_file} 2> /dev/null", shell = True)
    
    mums = parse_matches(mum_file)
    mems = parse_matches(mem_file)
    return mems, mums


def get_temp_file_name():
    return os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def plot_lines(drawing, matches, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
               line_width, color):
    
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

def plot_cigar(drawing, cigar, ref_begin, ref_end, query_begin, query_end, long_side,
               line_width, match_color, mismatch_color):
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    num_segments = 0

    curr_ref_pos = 0
    curr_query_pos = 0
    
    for op, op_len in cigar:
        
        if op in "MX=":
            next_ref_pos = curr_ref_pos + op_len
            next_query_pos = curr_query_pos + op_len
        elif op == "D":
            next_ref_pos = curr_ref_pos
            next_query_pos = curr_query_pos + op_len
        elif op in "IHS":
            next_ref_pos = curr_ref_pos + op_len
            next_query_pos = curr_query_pos
        else:
            assert(False)
            
        x_begin = (curr_query_pos - query_begin) * coord_multiplier
        x_end = (next_query_pos - query_begin) * coord_multiplier
        y_begin = (curr_ref_pos - ref_begin) * coord_multiplier
        y_end = (next_ref_pos - ref_begin) * coord_multiplier
        if op in "M=":
            color = match_color
        else:
            color = mismatch_color
        
        # only draw inside the 
        if ((curr_ref_pos >= ref_begin and curr_ref_pos < ref_end and curr_query_pos >= query_begin and curr_query_pos < query_end) or
            (next_ref_pos >= ref_begin and next_ref_pos < ref_end and next_query_pos >= query_begin and next_query_pos < query_end)):
        
            drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                     stroke = color, stroke_width = line_width))
            
            num_segments += 1
            
        curr_ref_pos = next_ref_pos
        curr_query_pos = next_query_pos
        
        
    print(f"added {num_segments} alignment segments", file = sys.stderr)
    
    return num_segments

def plot_axis_ticks(drawing, x_axis, begin, end, long_side, coord_multiplier):
    
    main_vals = [5, 10, 25]
    schedule = [main_vals[i % len(main_vals)] * (10**(i // len(main_vals))) for i in range(len(main_vals) * 8)]
    
    pixels = long_side / 75.0
    
    seq_len = end - begin
    
    optimum_num_ticks = 10
    tick_interval = 1
    for interval in schedule:
        if abs(seq_len / interval - optimum_num_ticks) < abs(seq_len / tick_interval - optimum_num_ticks):
            tick_interval = interval
            
    tick_begin = begin - (begin % tick_interval)
    
    for tick in range(tick_begin, end, tick_interval):
        if tick < begin:
            continue
     
        if x_axis:
            y_begin = 0
            y_end = long_side / 500.0
            x_begin = tick * coord_multiplier
            x_end = x_begin
        else:
            x_begin = 0
            x_end = long_side / 500.0
            y_begin = tick * coord_multiplier
            y_end = y_begin
            
        drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                 stroke = "black", stroke_width = long_side / 1000.0))
        
        if x_axis:
            x = x_end
            y = 2 * y_end + pixels
        else:
            x = 2 * x_end
            y = y_end
        
        drawing.add(drawing.text(str(int(tick)), insert = (x, y), style = "font-size:{}px".format(pixels)))
        

def plot_ticks(drawing, ref_begin, ref_end, query_begin, query_end, long_side):
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    
    plot_axis_ticks(drawing, False, query_begin, query_end, long_side, coord_multiplier)
    plot_axis_ticks(drawing, True, ref_begin, ref_end, long_side, coord_multiplier)
    
    
    

def get_matched_chroms(aln_identity_table):
    
    tab = pd.read_table(aln_identity_table, header = 0)
    
    matched_chrom = {}
    
    child = None
    for i in range(tab.shape[0]):
        chrom = tab.chr.values[i]
        for hap in (1, 2):
            max_j = None
            max_identity = None
            for j in range(4 * (hap - 1) + 1, 4 * hap + 1):
                identity = tab.iloc[i, j]
                if not np.isnan(identity) and (max_identity is None or max_identity < identity):
                    max_identity = identity
                    max_j = j
            
            if max_identity is not None and max_identity > 0.9:
                
                m = re.match("^(\w+)_(\d)_(\w+)_(\d)$", tab.columns[max_j])
                parent = m.group(3)
                parent_hap = int(m.group(4))
                if child is None:
                    child = m.group(1)
                assert child == m.group(1)
                
                matched_chrom[(chrom, hap)] = (parent, parent_hap)
                
    return child, matched_chrom

def parse_anchoring(table):
    
    anchors = []
    
    with open(table) as f:
        f.readline() # skip the header
        for line in f:
            tokens = line.strip().split()
            i = int(tokens[0])
            j = int(tokens[1])
            l = int(tokens[2])
            anchors.append((i, j, l))
    return anchors

def plot_anchoring(drawing, anchors, ref_begin, ref_end, query_begin, query_end, long_side,
                   line_width, connector_width, color):
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    num_segments = 0

    
    prev_x_end = None
    prev_y_end = None
    prev_inside = False
    for i, j, l in anchors:
        
            
        x_begin = (i - query_begin) * coord_multiplier
        x_end = (i + l - query_begin) * coord_multiplier
        y_begin = (j - ref_begin) * coord_multiplier
        y_end = (j + l - ref_begin) * coord_multiplier
    
        inside = ((j >= ref_begin and j < ref_end and i >= query_begin and i < query_end) or
                  (j + l >= ref_begin and j + l < ref_end and i + l >= query_begin and i + l < query_end))
        
        if inside:
            drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                     stroke = color, stroke_width = line_width))
            
            num_segments += 1
        
        if prev_x_end is not None and (inside or prev_inside):
            drawing.add(drawing.line((prev_x_end, prev_y_end), (x_begin, y_begin), 
                                 stroke = color, stroke_width = connector_width))
            num_segments += 1
            
        prev_x_end = x_end
        prev_y_end = y_end
        prev_inside = inside
            
        
    print(f"added {num_segments} anchoring segments", file = sys.stderr)  
    return num_segments      

if __name__ == "__main__":
    
    # ref or query according to the cigar operations
    ref_fasta = sys.argv[1]
    query_fasta = sys.argv[2]
    alignments = sys.argv[3]
    out = sys.argv[4]
    
    # i'm also repurposing this to plot centrolign anchoring
    is_cigar = True
    if len(sys.argv) >= 6:
        is_cigar = bool(int(sys.argv[5]))
    
    ref = parse_fasta(ref_fasta)
    query = parse_fasta(query_fasta)
    
    min_length = 768
    
    seq1_begin = 0
    seq1_end = len(ref)
    seq2_begin = 0
    seq2_end = len(query)
    
    long_side = 4000
    mem_line_width = 1
    mum_line_width = 2
    alignment_line_width = 1.5
    connector_width = 0.5
    
    mems, mums = get_matches(ref_fasta, query_fasta, min_length)
    
    print("loading sequences", file = sys.stderr)
    
    seq1_end = min(seq1_end, len(ref))
    seq2_end = min(seq2_end, len(query))
    
    seq1_len = seq1_end - seq2_begin
    seq2_len = seq2_end - seq2_begin
    
    
    
    multiplier1 = float(seq1_len) / max(seq1_len, seq2_len)
    multiplier2 = float(seq2_len) / max(seq1_len, seq2_len)
    
    height = long_side * multiplier1
    width = long_side * multiplier2
    
    print("initializing drawing", file = sys.stderr)
    drawing = svg.Drawing(out, size = (height, width))
    
    print("plotting MEMs", file = sys.stderr)
    
    num_lines = 0
    num_lines += plot_lines(drawing, mems, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                            mem_line_width, "black")
    
    print("plotting MUMs", file = sys.stderr)
    
    num_lines += plot_lines(drawing, mums, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                            mum_line_width, "lightcoral")
    
    alignment_colors = [("darkturquoise", "springgreen"), ("mediumpurple", "violet"),
                        ("peru", "sandybrown"), ("lightslategray", "lightsteelblue"), ("goldenrod", "gold")]
    
    if is_cigar:
        alignment_names = alignments.split(',')
        for i in range(len(alignment_names)):
            alignment = alignment_names[i]
            cigar = parse_cigar(open(alignment).read().strip())
        
            print("plotting alignment", file = sys.stderr)
            
            match_color, mismatch_color = alignment_colors[i % len(alignment_colors)]
            num_lines += plot_cigar(drawing, cigar, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                    alignment_line_width, match_color, mismatch_color)
    else:
        anchors = parse_anchoring(alignment)
        num_lines += plot_anchoring(drawing, anchors, seq1_begin, seq1_end, seq1_begin, seq1_end, long_side,
                                    alignment_line_width, connector_width, "forestlengreen")
    
    plot_ticks(drawing, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side)
    
    drawing.save()