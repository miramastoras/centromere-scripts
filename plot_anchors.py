#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:17:25 2022

@author: Jordan
"""

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
import random

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

# adjusts for the fact that angled lines bunch up together more
def weight(xs, ys):
    return math.sin(math.atan(abs(ys[1] - ys[0]) / (abs(xs[1] - xs[0]) + .00000000001)))

def filter_lines(lines, xmin, xmax):
    removed = 0
    for i in range(len(lines)):
        x_bottom, x_top = lines[i]
        if (x_top < xmin and x_bottom < xmin) or (x_top > xmax and x_bottom > xmax):
            removed += 1
        elif removed > 0:
            lines[i - removed] = lines[i]
    
    for i in range(removed):
        lines.pop()


def parse_pairs(filename):
    pairs = []
    for line in open(filename):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        pairs.append(tuple(map(int, line.strip().split())))
    return pairs

def plot_lines(ax, y1, y2, line_xs, linewidth, col, alpha, image_width, image_height, xlims, ylims, short_edge_len = None, short_edge_alpha_multiplier = 1.0, zorder = None):

    plotted = 0
    ys = [y1, y2]
    scaled_ys = [image_height * (y - ylims[0]) / (ylims[1] - ylims[0]) for y in ys]
    for xs in line_xs:
        scaled_xs = [image_width * (x - xlims[0]) / (xlims[1] - xlims[0]) for x in xs]
        scaled_alpha = weight(scaled_xs, scaled_ys) * alpha
        if short_edge_len is not None:
            if abs(xs[0] - xs[1]) < short_edge_len:
                scaled_alpha *= short_edge_alpha_multiplier
        ax.plot(xs, ys, "-", linewidth = linewidth, marker = None, zorder = zorder, alpha = scaled_alpha, c = col)
        plotted += 1
        if plotted % 100 == 0:
            print("plotted {} pairs".format(plotted), file = sys.stderr)

def filter_pairs(pairs, min_x, max_x, shift_bottom, shift_top):
    lines = []
    for i, j in pairs:
        ish = i - shift_bottom
        jsh = j - shift_top
        if not ((ish < min_x and jsh < min_x) or (ish > max_x and jsh > max_x)):
            lines.append([ish, jsh])
    return lines
    
    
if __name__ == "__main__":
    
    rejected_pairs = []
    unique_pairs = []
    
    
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    anchors = sys.argv[3]
    out = sys.argv[4]
    if len(sys.argv) >= 6:
        rejected_pairs = parse_pairs(sys.argv[5])
        print("loaded {} rejected pairs".format(len(rejected_pairs)), file = sys.stderr)
    if len(sys.argv) >= 7:
        unique_pairs = parse_pairs(sys.argv[6])
        print("loaded {} unique pairs".format(len(unique_pairs)), file = sys.stderr)
    self1_pairs = []
    self2_pairs = []
    self1_unique_pairs = []
    self2_unique_pairs = []
    if len(sys.argv) >= 8:
        self1_pairs = parse_pairs(sys.argv[7])
        self2_pairs = parse_pairs(sys.argv[8])
        self1_unique_pairs = parse_pairs(sys.argv[9])
        self2_unique_pairs = parse_pairs(sys.argv[10])
            
    seq_height = .1
    y0 = 2.0 - 2 * seq_height
    y1 = 1.0 - seq_height
    y2 = 0.0
    y3 = -1.0 + seq_height
    kmer_width = .0001
    resolution = 400
    height = 12
    width = 80
    connecting_line_width = 0.15
    line_alpha_retained = 0.5
    line_alpha_removed = 0.5
    line_alpha_unique = 0.2
    line_alpha_self = 0.6
    line_alpha_self_unique = 0.5
    
    short_edge_alpha_multiplier = 0.1
    short_edge_len = None#20000
    
    max_lines = 12000
    max_unique_lines = 2000
    max_self_nonunique_lines = 8000
    max_self_unique_lines = 2000
    filter_nonunique_self_identities = True
    
    view_min = 900000
    view_len = 50000
        
    only_unique = False
    view_percentile = 20
    hard_xlim = (910000, 990000)
    override_shift = 147500#None
    
    retained_col = "r"
    removed_col = "r"
    unique_col = "0.25"
    
    seq1 = parse_fasta(fasta1)
    seq2 = parse_fasta(fasta2)
    len1 = len(seq1)
    len2 = len(seq2)
    
    pairs = []
    
    for line in open(anchors):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        pairs.append(tuple(map(int, line.strip().split())))
        
    print("loaded {} selected pairs".format(len(pairs)), file = sys.stderr)
    
    f, ax = plt.subplots(1, 1, figsize = (width, height), dpi = resolution)
    
    x1 = 0.0
    x2 = 0.0
    
    
    selected = 0
    
    num_pairs = 0
    
    removed_lines = []
    retained_lines = []
    
    print("seeding with lines from top interval {} to {}".format(view_min, view_min + view_len), file = sys.stderr)    

    min_seed = len(seq1) + len(seq2) + 1
    max_seed = -1
    for i, j in pairs:    
            
#            ax.plot([xp1, xp2], [y1, y2 + seq_height], "r-,", linewidth = connecting_line_width, marker = None, zorder = 0, alpha = line_alpha_removed)
        
        
        xp1 = i
#        pos_rect = patches.Rectangle((xp1, y1), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        xp2 = j
#        pos_rect = patches.Rectangle((xp2, y2), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        if xp1 >= view_min and xp1 <= view_min + view_len:
        
#        if not ((xp1 < view_min and xp2 < view_min) or (xp1 > view_min + view_len and xp2 > view_min + view_len)):
            retained_lines.append([xp1, xp2])
            min_seed = min(min_seed, i)
            max_seed = max(max_seed, i)
                
            selected += 1
            
    
    for ci, cj in rejected_pairs:
        
        xp1 = ci
#        pos_rect = patches.Rectangle((xp1, y1), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        xp2 = cj
#        pos_rect = patches.Rectangle((xp2, y2), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        if xp1 >= view_min and xp1 <= view_min + view_len:
#        if not ((xp1 < view_min and xp2 < view_min) or (xp1 > view_min + view_len and xp2 > view_min + view_len)):
        
            removed_lines.append([xp1, xp2])
            min_seed = min(min_seed, ci)
            max_seed = max(max_seed, ci)
            
            selected += 1
    print("seeded with {} retained pairs and {} rejected pairs from seq1 interval {}:{}".format(len(retained_lines), len(removed_lines), min_seed, max_seed), file = sys.stderr)
    
    # just to make the formula work
    assert(len(removed_lines) + len(retained_lines) >= 2)
    all_seed_lines = sorted(removed_lines + retained_lines)
    num_end_shift = 1000
    bnd_1 = min(num_end_shift, len(all_seed_lines)//2)
    bnd_2 = max(len(all_seed_lines) - num_end_shift, len(all_seed_lines)//2)
    shift_left_end = np.median([xs[1] - xs[0] for xs in all_seed_lines[:bnd_1]])
    shift_right_end = np.median([xs[1] - xs[0] for xs in all_seed_lines[bnd_2:]])
    
    print("left and right ends are shifted {} and {} relative to top".format(shift_left_end, shift_right_end), file = sys.stderr)
    print("keeping lines crossing into window {} and {}".format(view_min + shift_left_end, view_min + view_len + shift_right_end), file = sys.stderr)    

    for i, j in pairs:    
            
#            ax.plot([xp1, xp2], [y1, y2 + seq_height], "r-,", linewidth = connecting_line_width, marker = None, zorder = 0, alpha = line_alpha_removed)
        
        
        xp1 = i
#        pos_rect = patches.Rectangle((xp1, y1), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        xp2 = j
#        pos_rect = patches.Rectangle((xp2, y2), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        
        
        if (xp1 < view_min and xp2 > view_min + shift_left_end) or (xp1 > view_min + view_len and xp2 < view_min + view_len + shift_right_end):
        
            retained_lines.append([xp1, xp2])
                
            selected += 1
            
    
    for ci, cj in rejected_pairs:
        
        xp1 = ci
#        pos_rect = patches.Rectangle((xp1, y1), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        xp2 = cj
#        pos_rect = patches.Rectangle((xp2, y2), kmer_width, seq_height, linewidth=0, edgecolor=None, facecolor='r')
#        ax.add_patch(pos_rect)
        
        if (xp1 < view_min and xp2 > view_min + shift_left_end) or (xp1 > view_min + view_len and xp2 < view_min + view_len + shift_right_end):
        
            removed_lines.append([xp1, xp2])
            selected += 1
    
    
    print("selected {} retained, {} removed pairs ({} total)".format(len(retained_lines), len(removed_lines), len(retained_lines)+len(removed_lines)), file=sys.stderr)
    
    shift = np.median([xs[1] - xs[0] for xs in retained_lines + removed_lines])
    if override_shift is not None:
        shift = override_shift
    
    print("shifting bottom row by {}".format(shift), file=sys.stderr)
    
    all_xs = []
    for xs in retained_lines:
        xs[1] -= shift
        all_xs.extend(xs)
    for xs in removed_lines:
        xs[1] -= shift
        all_xs.extend(xs)
    
    seq_rect1 = patches.Rectangle((x1, y1), len1, seq_height, linewidth=0, edgecolor=None, facecolor='k')
    seq_rect2 = patches.Rectangle((x2 - shift, y2), len2, seq_height, linewidth=0, edgecolor=None, facecolor='k')
    ax.add_patch(seq_rect1)
    ax.add_patch(seq_rect2)
    
    if hard_xlim is not None:
        filter_lines(retained_lines, hard_xlim[0], hard_xlim[1])
        filter_lines(removed_lines, hard_xlim[0], hard_xlim[1])
        print("reduced to lines in hard view limit: {} retained and {} removed pairs ({} total)".format(len(retained_lines), len(removed_lines), len(retained_lines) + len(removed_lines)), file=sys.stderr)
    
    if len(retained_lines) + len(removed_lines) > max_lines:
        sub_samp_retained = max_lines * len(retained_lines) // (len(retained_lines) + len(removed_lines))
        sub_samp_removed = max_lines * len(removed_lines) // (len(retained_lines) + len(removed_lines))
        
        if sub_samp_retained < len(retained_lines):
            retained_lines = sorted(random.sample(retained_lines, sub_samp_retained))
        if sub_samp_removed < len(removed_lines):
            removed_lines = sorted(random.sample(removed_lines, sub_samp_removed))
            
        print("downsampled to {} non-unique pairs".format(len(removed_lines) + len(retained_lines)), file=sys.stderr)
        
        
    if filter_nonunique_self_identities:
        self1_pairs = [xs for xs in self1_pairs if xs[0] != xs[1]]
        self2_pairs = [xs for xs in self2_pairs if xs[0] != xs[1]]
    
    print("ranges: {}-{} {}-{}".format(min(xs[0] for xs in retained_lines),
                                       max(xs[0] for xs in retained_lines),
                                       min(xs[1] for xs in retained_lines),
                                       max(xs[1] for xs in retained_lines)))
        
    plo = np.percentile(all_xs, view_percentile)
    phi = np.percentile(all_xs, 100 - view_percentile)
    print("percentile range {} : {}".format(plo, phi), file=sys.stderr)
    
    if hard_xlim is  None:    
        min_x = plo - (phi - plo) * 0.05
        max_x = phi + (phi - plo) * 0.05
    else:
        min_x, max_x = hard_xlim
        
        
    unique_lines = filter_pairs(unique_pairs, min_x, max_x, shift, 0)
    
    print("selected {} unique pairs".format(len(unique_lines)), file=sys.stderr)
    
    if len(unique_lines) > max_unique_lines:
        unique_lines = sorted(random.sample(unique_lines, max_unique_lines))
        print("downsampled to {} unique pairs".format(len(unique_lines)), file=sys.stderr)
    
    self1_lines = filter_pairs(self1_pairs, min_x, max_x, 0, 0)
    self2_lines = filter_pairs(self2_pairs, min_x, max_x, shift, shift)
    self1_unique_lines = filter_pairs(self1_unique_pairs, min_x, max_x, 0, 0)
    self2_unique_lines = filter_pairs(self2_unique_pairs, min_x, max_x, shift, shift)
    
    if len(self1_lines) + len(self2_lines) > max_self_nonunique_lines:
        sub_samp1 = max_self_nonunique_lines * len(self1_lines) // (len(self1_lines) + len(self2_lines))
        sub_samp2 = max_self_nonunique_lines * len(self2_lines) // (len(self1_lines) + len(self2_lines))
        
        if sub_samp1 < len(self1_lines):
            self1_lines = sorted(random.sample(self1_lines, sub_samp1))
        if sub_samp2 < len(self2_lines):
            self2_lines = sorted(random.sample(self2_lines, sub_samp2))
            
        print("downsampled to {} non-unique self pairs".format(len(self1_lines) + len(self2_lines)), file=sys.stderr) 
    
    if len(self1_unique_lines) + len(self2_unique_lines) > max_self_unique_lines:
        sub_samp1 = max_self_unique_lines * len(self1_unique_lines) // (len(self1_unique_lines) + len(self2_unique_lines))
        sub_samp2 = max_self_unique_lines * len(self2_unique_lines) // (len(self1_unique_lines) + len(self2_unique_lines))
        
        if sub_samp1 < len(self1_unique_lines):
            self1_unique_lines = sorted(random.sample(self1_unique_lines, sub_samp1))
        if sub_samp2 < len(self2_unique_lines):
            self2_unique_lines = sorted(random.sample(self2_unique_lines, sub_samp2))
            
        print("downsampled to {} unique self pairs".format(len(self1_unique_lines) + len(self2_unique_lines)), file=sys.stderr) 
    
    
    xlims = [min_x - .05 * (max_x - min_x), max_x + 0.05 * (max_x - min_x)]
    ylims = [-.05, 1.05]
    if len(self1_pairs) + len(self1_unique_pairs) + len(self2_pairs) + len(self2_unique_pairs) > 0:
        ylims = [-1.05, 2.05]
        self_rect1 = patches.Rectangle((x1, y0), len1, seq_height, linewidth=0, edgecolor=None, facecolor='k')
        self_rect2 = patches.Rectangle((x2 - shift, y3), len2, seq_height, linewidth=0, edgecolor=None, facecolor='k')
        ax.add_patch(self_rect1)
        ax.add_patch(self_rect2)
    
    ax2 = ax.twiny()
    
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_xlim(xlims[0], xlims[1])
    ax2.set_xlim(xlims[0], xlims[1])
    
    ax.tick_params(axis = 'x', labelsize = int(height * 2))
    ax2.tick_params(axis = 'x', labelsize = int(height * 2))
    
    ax.xaxis.set_label_position('top') 
    ax.xaxis.tick_top()
    ax2.xaxis.set_label_position('bottom') 
    ax2.xaxis.tick_bottom()
    
    ax2.set_xticklabels([str(int(round(x + shift))) for x in ax.get_xticks()])
        
    ys = [y1, y2 + seq_height]
    
    scaled_ys = [height * (y - ylims[0]) / (ylims[1] - ylims[0]) for y in ys]

    plotted = 0
    if not only_unique:
        print("plotting retained lines", file=sys.stderr)
        plot_lines(ax, y1, y2 + seq_height, retained_lines, connecting_line_width, retained_col,
                   line_alpha_retained, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier, None)

            
        print("plotting removed lines", file=sys.stderr)
        plot_lines(ax, y1, y2 + seq_height, removed_lines, connecting_line_width, removed_col,
                   line_alpha_removed, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier, 0)

        
    print("plotting unique lines", file=sys.stderr)
    plot_lines(ax, y1, y2 + seq_height, unique_lines, connecting_line_width, unique_col,
               line_alpha_unique, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier, -1)
 
    
    print("plotting unique self lines", file=sys.stderr)
    plot_lines(ax, y0, y1 + seq_height, self1_unique_lines, connecting_line_width, unique_col,
               line_alpha_self_unique, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier, -1)
    plot_lines(ax, y2, y3 + seq_height, self2_unique_lines, connecting_line_width, unique_col,
               line_alpha_self_unique, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier, -1)
    
    if not only_unique:
        print("plotting non-unique self lines", file=sys.stderr)
        plot_lines(ax, y0, y1 + seq_height, self1_lines, connecting_line_width, removed_col,
                   line_alpha_self, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier)
        plot_lines(ax, y2, y3 + seq_height, self2_lines, connecting_line_width, removed_col,
                   line_alpha_self, width, height, xlims, ylims, short_edge_len, short_edge_alpha_multiplier)
    
    f.savefig(out, format = "svg", bbox_inches = "tight")
    
    plt.close(f)