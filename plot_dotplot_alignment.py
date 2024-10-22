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
import collections
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


def group_matches(matches):
    
    end_to_group = {}
    groups = []
    for i, j, l in sorted(matches):
        
        key_i = (i, True, l)
        key_j = (j, False, l)
        
        group = -1
        if key_i in end_to_group:
            group = end_to_group[key_i]
        elif key_j in end_to_group:
            group = end_to_group[key_j]
            
        if group == -1:
            group = len(groups)
            groups.append([])
        
        groups[group].append((i, j, l))
        end_to_group[key_i] = group
        end_to_group[key_j] = group
        
    
    return groups
            
        
    

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

def plot_lines_with_count_shading(drawing, matches, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                  line_width, count_shading_power, white_cutoff = 0):
    
    groups = group_matches(matches)
    
    print("{} matches come in {} groups".format(len(matches), len(groups)), file = sys.stderr)
    #print(sum(len(g) for g in groups))
    
    coord_multiplier = float(long_side) / max(seq1_end - seq1_begin, seq2_end - seq2_begin)
    num_lines = 0
    
    for group in groups:
        
        opacity = 1.0 / len(group)**count_shading_power
        
        val = round((1.0 - opacity) * 255)
        color = svg.rgb(val, val, val)
        
        if val < white_cutoff:
            for i, j, length in group: 
                x = (i - seq1_begin) * coord_multiplier
                y = (j - seq2_begin) * coord_multiplier
                l = length * coord_multiplier
                
                
                drawing.add(drawing.line((x, y), (x + l, y + l), 
                                         stroke = color, stroke_width = line_width, opacity = opacity))
                
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
    
    d_max = -1
    i_max = -1
    for op, op_len in cigar:
        
        if op in "MX=":
            next_ref_pos = curr_ref_pos + op_len
            next_query_pos = curr_query_pos + op_len
        elif op == "D":
            next_ref_pos = curr_ref_pos
            next_query_pos = curr_query_pos + op_len
            d_max = max(op_len, d_max)
        elif op in "IHS":
            next_ref_pos = curr_ref_pos + op_len
            next_query_pos = curr_query_pos
            i_max = max(op_len, i_max)
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
        
        # # only draw inside the 
        # if ((curr_ref_pos >= ref_begin and curr_ref_pos < ref_end and curr_query_pos >= query_begin and curr_query_pos < query_end) or
        #     (next_ref_pos >= ref_begin and next_ref_pos < ref_end and next_query_pos >= query_begin and next_query_pos < query_end)):
        
        drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                 stroke = color, stroke_width = line_width))
        
        num_segments += 1
            
        curr_ref_pos = next_ref_pos
        curr_query_pos = next_query_pos
        
        
    print(f"added {num_segments} alignment segments", file = sys.stderr)
    print(f"max deletion length {d_max}, max insertion length {i_max}", file = sys.stderr)
    
    return num_segments

def plot_windowed_identity_cigar(drawing, cigar, ref_begin, ref_end, query_begin, query_end, long_side,
                                 line_width, window_len, block_size = 10):
    
        
    segments = []
    curr_ref_pos = 0
    curr_query_pos = 0
    for op, op_len in cigar:
        match = (op in "M=")
        if op in "MX=":
            incr_ref = 1
            incr_query = 1
        elif op == "D":
            incr_ref = 0
            incr_query = 1
        elif op in "IHS":
            incr_ref = 1
            incr_query = 0
        else:
            assert(False)
            
        for i in range(op_len):
            segments.append([(curr_query_pos, curr_ref_pos), 
                             (curr_query_pos + incr_query, curr_ref_pos + incr_ref), match])
            curr_query_pos += incr_query
            curr_ref_pos += incr_ref
            
    
    print("parsed {} alignment positions".format(len(segments)), file = sys.stderr)
    
    # init the window
    window = collections.deque()
    window_end = window_len // 2 - 1
    window_count = 0
    for i in range(max(0, min(window_end, len(segments)))):
        window.append(segments[i][-1])
        if segments[i][-1]:
            window_count += 1
            
    
    for i in range(len(segments)):
        
        if window_end < len(segments):
            window.append(segments[window_end][-1])
            if segments[window_end][-1]:
                window_count += 1
        
        window_end += 1
        
        if window_end - window_len > 0:
            is_match = window.popleft()
            if is_match:
                window_count -= 1
                
        window_prop_match = window_count / len(window)
        segments[i][-1] = window_prop_match
        
    print("identified windowed match proportion", file = sys.stderr)
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    
    r_lo = 10
    g_lo = 15
    b_lo = 180
    r_hi = 50
    g_hi = 100
    b_hi = 10
    
    num_lines = 0
    for i in range(0, len(segments), block_size):
        j = i
        avg_prop = 0.0
        while j < min(len(segments), i + block_size):
            start, end, prop_match = segments[i]
            avg_prop += prop_match
            j += 1
        avg_prop /= (j - i)
        
        x_begin = (segments[i][0][0] - query_begin) * coord_multiplier
        x_end = (segments[j-1][1][0] - query_begin) * coord_multiplier
        y_begin = (segments[i][0][1] - ref_begin) * coord_multiplier
        y_end = (segments[j-1][1][1] - ref_begin) * coord_multiplier
        
        r = round(r_lo + (r_hi - r_lo) * avg_prop)
        g = round(g_lo + (g_hi - g_lo) * avg_prop)
        b = round(b_lo + (b_hi - b_lo) * avg_prop)
    
        drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                 stroke = svg.rgb(r, g, b), stroke_width = line_width, opacity = 1, fill_opacity = 1))
        
        # if num_lines % 10 == 0:
        #     print(svg.rgb(r, g, b), avg_prop , file = sys.stderr)
        
        num_lines += 1
        if num_lines % 10000 == 0:
            print("added {} lines".format(num_lines, svg.rgb(r, g, b)), file = sys.stderr)
        
    print("added {} lines".format(num_lines), file = sys.stderr)

def plot_axis_ticks(drawing, x_axis, begin, end, other_begin, other_end, long_side, coord_multiplier, do_grid = False):
    
    main_vals = [5, 10, 25]
    schedule = [main_vals[i % len(main_vals)] * (10**(i // len(main_vals))) for i in range(len(main_vals) * 8)]
    
    pixels = long_side / 75.0
    
    seq_len = max(end - begin, other_end - other_begin)
    
    major_tick_width = long_side / 1000.0
    minor_tick_width = major_tick_width / 10.0
    subminor_tick_width = minor_tick_width / 5.0
    major_grid_line_width = minor_tick_width * 2.0
    major_tick_length = long_side / 500.0
    if do_grid:
        minor_tick_length = (other_end - other_begin) * coord_multiplier
    else:
        minor_tick_length = major_tick_length / 1.5
    
    optimum_num_ticks = 10
    major_tick_interval = 1
    
    for interval in schedule:
        if abs(seq_len / interval - optimum_num_ticks) < abs(seq_len / major_tick_interval - optimum_num_ticks):
            major_tick_interval = interval
            
    minor_tick_interval = major_tick_interval // 5
    subminor_tick_interval = minor_tick_interval // 5
    assert(major_tick_interval % minor_tick_interval == 0)
    assert(minor_tick_interval % subminor_tick_interval == 0)
    
    
    major_tick_begin = begin - (begin % major_tick_interval)
    subminor_tick_begin = begin - (begin % subminor_tick_interval)
    
    for tick in range(subminor_tick_begin, end, subminor_tick_interval):
        
        if tick < begin:
            continue
        
        if tick == 0:
            continue
        
        if x_axis:
            y_begin = 0
            y_end = minor_tick_length
            x_begin = tick * coord_multiplier
            x_end = x_begin
        else:
            x_begin = 0
            x_end = minor_tick_length
            y_begin = tick * coord_multiplier
            y_end = y_begin
            
        width = subminor_tick_width
        if tick % minor_tick_interval == 0:
            width = minor_tick_width
        if do_grid and tick % major_tick_interval == 0:
            width = major_grid_line_width
        
        drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                 stroke = "black", stroke_width = width))
    
    for tick in range(major_tick_begin, end, major_tick_interval):
        if tick < begin:
            continue
     
        if x_axis:
            y_begin = 0
            y_end = major_tick_length
            x_begin = tick * coord_multiplier
            x_end = x_begin
        else:
            x_begin = 0
            x_end = major_tick_length
            y_begin = tick * coord_multiplier
            y_end = y_begin
            
            
        drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                 stroke = "black", stroke_width = major_tick_width))
        
        if x_axis:
            x = x_end
            y = 2 * y_end + pixels
        else:
            x = 2 * x_end
            y = y_end
        
        drawing.add(drawing.text(str(int(tick)), insert = (x, y), style = "font-size:{}px".format(pixels)))
    
    

def plot_ticks(drawing, ref_begin, ref_end, query_begin, query_end, long_side, add_grid = False):
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    
    plot_axis_ticks(drawing, False, query_begin, query_end, ref_begin, ref_end, long_side, coord_multiplier, add_grid)
    plot_axis_ticks(drawing, True, ref_begin, ref_end, query_begin, query_end, long_side, coord_multiplier, add_grid)
    
    
    

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
            level = 0
            if len(tokens) > 3:
                level = int(tokens[3])
            anchors.append((i, j, l, level))
    return anchors

def plot_anchoring(drawing, anchors, ref_begin, ref_end, query_begin, query_end, long_side,
                   line_width, connector_width, bond_shrinkage, colors, bond_colors):
    
    coord_multiplier = float(long_side) / max(ref_end - ref_begin, query_end - query_begin)
    num_segments = 0

    
    prev_x_end = None
    prev_y_end = None
    prev_inside = False
    prev_level = None
    for i, j, l, level in anchors:
        local_line_width = line_width
        local_connector_width = connector_width
        if level >= 0:
            color = colors[level % len(colors)]
        else:
            color = bond_colors[(-level - 1) % len(bond_colors)]
            local_line_width *= bond_shrinkage
            local_connector_width *= bond_shrinkage
        
            
        x_begin = (i - query_begin) * coord_multiplier
        x_end = (i + l - query_begin) * coord_multiplier
        y_begin = (j - ref_begin) * coord_multiplier
        y_end = (j + l - ref_begin) * coord_multiplier
    
        inside = ((j >= ref_begin and j < ref_end and i >= query_begin and i < query_end) or
                  (j + l >= ref_begin and j + l < ref_end and i + l >= query_begin and i + l < query_end))
        
        # TODO: i'm not actually sure what the bug here is, but it sometimes removes segments
        inside = True
        if inside:
            drawing.add(drawing.line((x_begin, y_begin), (x_end, y_end), 
                                     stroke = color, stroke_width = local_line_width))
            
            num_segments += 1
        
        if prev_x_end is not None and (inside or prev_inside) and level == prev_level:
            drawing.add(drawing.line((prev_x_end, prev_y_end), (x_begin, y_begin), 
                                 stroke = color, stroke_width = local_connector_width))
            num_segments += 1
            
        prev_x_end = x_end
        prev_y_end = y_end
        prev_inside = inside
        prev_level = level
            
        
    print(f"added {num_segments} anchoring segments", file = sys.stderr)  
    return num_segments      

if __name__ == "__main__":
    
    if len(sys.argv) not in [5, 6]:
        print("usage:", file = sys.stderr)
        print("./plot_dotplot_alignment.py fasta1 fasta2 cigar[,cigar2,cigar3,...] svg_out_name", file = sys.stderr)
        print("or", file = sys.stderr)
        print("./plot_dotplot_alignment.py fasta1 fasta2 \"\" svg_out_name", file = sys.stderr)
        exit(1)
    
    # ref or query according to the cigar operations
    ref_fasta = sys.argv[1]
    query_fasta = sys.argv[2]
    alignments = sys.argv[3] # leave empty ("") to only plot MEMs
    out = sys.argv[4]
    
    # i'm also repurposing this to plot centrolign anchoring
    is_cigar = True
    if len(sys.argv) >= 6:
        is_cigar = bool(int(sys.argv[5]))
    
    ref = parse_fasta(ref_fasta)
    query = parse_fasta(query_fasta)
    
    # include sequence position ticks
    add_ticks = True
    
    # the window to plot
    seq1_begin = 0
    seq1_end = len(ref)
    seq2_begin = 0
    seq2_end = len(query)
    
    # size of the image in pixels
    long_side = 3000
    
    # min length match (for MEMs and MUMs) to be plotted
    min_length = 700
    
    # plot MUMs in a separate color from MEMs
    plot_mums = True
    mum_color = "crimson"
    
    # sizes for lines corresponding to matches and alignments
    mem_line_width = 2.1
    mum_line_width = 2.1
    alignment_line_width = 1.5
    connector_width = 0.75 # for connecting anchors in anchoring plot
    bond_shrinkage = 0.5
    
    # shade matches according to their uniqueness
    use_count_shading = True
    # how steeply to shade matches (closer to 0 -> less steep)
    count_shading_power = 0.5
    
    add_grid = True
    
    # shade the alignment acording to its % identity within a window
    do_windowed = False
    window = 2
    
    mems, mums = get_matches(ref_fasta, query_fasta, min_length)
    
    print("loading sequences", file = sys.stderr)
    
    seq1_end = min(seq1_end, len(ref))
    seq2_end = min(seq2_end, len(query))
    
    
    seq1_len = seq1_end - seq2_begin
    seq2_len = seq2_end - seq2_begin
    
    print("plotting sequences of length {} (x-axis) and {} (y-axis)".format(seq1_len, seq2_len), file = sys.stderr)
    
    multiplier1 = float(seq1_len) / max(seq1_len, seq2_len)
    multiplier2 = float(seq2_len) / max(seq1_len, seq2_len)
    
    height = long_side * multiplier1
    width = long_side * multiplier2
    
    print("initializing drawing", file = sys.stderr)
    drawing = svg.Drawing(out, size = (height, width))
    
    print("plotting MEMs", file = sys.stderr)
    
    num_lines = 0
    if use_count_shading:
        num_lines += plot_lines_with_count_shading(drawing, mems, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                                   mem_line_width, count_shading_power, 240)
    else:
        num_lines += plot_lines(drawing, mems, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                mem_line_width, "black")
    
    
    if plot_mums:
        
        print("plotting MUMs", file = sys.stderr)
        
        num_lines += plot_lines(drawing, mums, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                mum_line_width, mum_color)
    
    alignment_colors = [("darkturquoise", "springgreen"), ("mediumpurple", "violet"),
                        ("peru", "sandybrown"), ("goldenrod", "gold"), ("lightslategray", "lightsteelblue")]
    
    if len(alignments.strip()) != 0:
        alignment_names = alignments.strip().split(',')
        if is_cigar:
            if len(alignment_names) > 1 or not do_windowed:
                for i in range(len(alignment_names)):
                    alignment = alignment_names[i]
                    cigar = parse_cigar(open(alignment).read().strip())
                
                    print("plotting alignment", file = sys.stderr)
                    
                    match_color, mismatch_color = alignment_colors[i % len(alignment_colors)]
                    num_lines += plot_cigar(drawing, cigar, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                            alignment_line_width, match_color, mismatch_color)
            else:
                print("plotting windowed alignment", file = sys.stderr)
                cigar = parse_cigar(open(alignment_names[0]).read().strip())
                plot_windowed_identity_cigar(drawing, cigar, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                             alignment_line_width, window)
                
        else:
            assert(len(alignment_names) == 1)
            alignment = alignment_names[0]
            anchors = parse_anchoring(alignment)
            anchor_colors = ["forestgreen", "darkturquoise", "mediumpurple", "peru", "coral"]
            bond_colors = ["deeppink", "lawngreen", "dodgerblue", "darkkhaki", "crimson", "aqua", "red", "greenyellow", 
                           "olivedrab", "aquamarine", "hotpink", "lightblue", "darksalmon", "midnightblue", "blueviolet"]
            num_lines += plot_anchoring(drawing, anchors, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side,
                                        alignment_line_width, connector_width, bond_shrinkage, anchor_colors, bond_colors)
    
    if add_ticks:
        plot_ticks(drawing, seq1_begin, seq1_end, seq2_begin, seq2_end, long_side, add_grid)
    
    drawing.save()