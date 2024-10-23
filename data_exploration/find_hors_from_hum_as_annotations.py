#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:47:31 2023

@author: Jordan
"""

import sys
import re

if __name__ == "__main__":
    
    if len(sys.argv) != 5:
        print("usage: ./find_hors_from_hum_as_annotations.py annot.bed merge_dist min_length min_cov_fraction", file = sys.stderr)
        sys.exit(1)
    
    hum_as_bed = sys.argv[1]
    slop_distance = int(sys.argv[2])
    min_length = int(sys.argv[3])
    min_frac = float(sys.argv[4])
    
    
    curr_contig = None
    
    curr_intervals = {}
    
        
    
    for line in open(hum_as_bed):
        tokens = line.strip().split()
        contig = tokens[0]
        begin = int(tokens[1])
        end = int(tokens[2])
        
        if contig != curr_contig:
            # flush the current intervals
            for hor in curr_intervals:
                curr_begin, curr_end, cov = curr_intervals[hor]
                if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
                    print(f"{curr_contig}\t{curr_begin}\t{curr_end}\t{hor}")
            curr_intervals.clear()
            curr_contig = contig
        
        
        m = re.match("(S\d+C[\dXY/]+H\d+L)\.[\d/]+", tokens[3])

        if m is None:
            continue
        
        hor = m.group(1)
        
        if hor in curr_intervals:
            curr_begin, curr_end, cov = curr_intervals[hor]
            if curr_end + slop_distance >= begin:
                curr_intervals[hor][1] = end
                curr_intervals[hor][2] += end - begin
            else:
                if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
                    print(f"{curr_contig}\t{curr_begin}\t{curr_end}\t{hor}")
                curr_intervals[hor] = [begin, end, end - begin]
        else:
            curr_intervals[hor] = [begin, end, end - begin]
        
        
    
    for hor in curr_intervals:
        curr_begin, curr_end, cov = curr_intervals[hor]
        if curr_end - curr_begin >= min_length and cov / (curr_end - curr_begin) > min_frac:
            print(f"{curr_contig}\t{curr_begin}\t{curr_end}\t{hor}")
            
        
        
    

