#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:00:23 2022

@author: Jordan
"""

import sys
import pandas as pd
import numpy as np
import os
import subprocess
import tempfile
import collections

def tmp_file_name(suffix = None):
    tmp = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    if suffix is not None:
        tmp += "." + suffix
    return tmp

def rev_comp(seq):
    comp = {"A":"T", "a":"t", "C":"G", "c":"g", "G":"C", "g":"c", "T":"A", "t":"a", "N":"N", "n":"n"}
    return "".join(comp[c] for c in reversed(seq))

if __name__ == "__main__":
    
    assembly = sys.argv[1]
    flank_table = pd.read_table(sys.argv[2], header = 0)
    min_coverage = float(sys.argv[3])
    output_dir = sys.argv[4]
    
#    print(flank_table)
    
    if os.path.exists(output_dir):
        print(f"error, output directory {output_dir} already exists")
        exit(1)
        
    os.mkdir(output_dir)
    
    chrom_counts = flank_table.centro_chrom.value_counts()
    curr_chrom_num = collections.defaultdict(int)
    
    for i in range(flank_table.shape[0]):
        chrom = flank_table.centro_chrom.values[i]
        curr_chrom_num[chrom] += 1
        
        if np.isnan(flank_table.left_contig_coverage.values[i]):
            print("warning: skipping centromere {} on chromosome {} for failed flank projection".format(curr_chrom_num[chrom], chrom), file = sys.stderr)
            continue
        
        if (flank_table.left_contig_coverage.values[i] < min_coverage or
            flank_table.right_contig_coverage.values[i] < min_coverage or 
            flank_table.left_flank_coverage.values[i] < min_coverage or
            flank_table.right_flank_coverage.values[i] < min_coverage):
            
            print("warning: skipping centromere {} on chromosome {} for low-quality flank mapping".format(curr_chrom_num[chrom], chrom), file = sys.stderr)
            continue
        
        
        contig = flank_table.flank_contig.values[i]
        reverse = bool(flank_table.is_reverse.values[i])
        if reverse:
            left = int(flank_table.right_flank_end.values[i])
            right = int(flank_table.left_flank_begin.values[i])
        else:
            left = int(flank_table.left_flank_end.values[i])
            right = int(flank_table.right_flank_begin.values[i])
    
        
        if chrom_counts[chrom] != 1:
            fasta_file_stem = chrom + "." + str(curr_chrom_num[chrom]) + ".fasta"
            fasta_file_stem = chrom + ".fasta"
            
        final_fasta = os.path.join(output_dir, fasta_file_stem)
        
        if reverse:
            extract_fasta = tmp_file_name("fasta")
        else:
            extract_fasta = final_fasta
        
        region_file = tmp_file_name("txt")
        with open(region_file, "w") as f:
            print(f"{contig}:{left}-{right}", file = f)
        
        subprocess.check_call(f"samtools faidx -r {region_file} -o {extract_fasta} {assembly}", shell = True)
        
        os.remove(region_file)
        
        if reverse:
            # put the fasta in the forward orientation relative to the CHM flanks
            with open(final_fasta, "w") as f_out:
                length_count = collections.defaultdict(int)
                seq = ""
                with open(extract_fasta) as f_in:
                    for line in f_in:
                        if type(line) == bytes:
                            line = line.decode("utf-8")
                        line = line.strip()
                        if line.startswith(">"):
                            line += " revcomp=True"
                            print(line, file = f_out)
                        else:
                            seq += line
                            length_count[len(line)] += 1
                line_len = max(length_count, key = lambda l : length_count[l])
                seq = rev_comp(seq)
                i = 0
                while i < len(seq):
                    print(seq[i:min(i+line_len, len(seq))], file = f_out)
                    i += line_len
            os.remove(extract_fasta)
                
                        
        
        
        