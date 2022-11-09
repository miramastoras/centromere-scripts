#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:35:04 2022

@author: Jordan
"""

import sys
import re
import os
import subprocess


def get_centromere_files(sample_prefix):
    
    files = {}
    
    for hap in (1, 2):
        array_dir = sample_prefix + "_" + str(hap) + "_active_arrays"
        
        for filename in os.listdir(array_dir):
            chrom = re.search("(\S*).fa", filename).group(1)
            
            if chrom not in files:
                files[chrom] = {}
            
            files[chrom][hap] = os.path.abspath(os.path.join(array_dir, filename))
    
    return files
        

if __name__ == "__main__":
    
    child = sys.argv[1]
    parent1 = sys.argv[2]
    parent2 = sys.argv[3]
    output_dir = sys.argv[4]
    
    if os.path.isdir(output_dir):
        print("error: output directory {} already exists".format(output_dir), file = sys.stderr)
        exit()
    else:
        os.makedirs(output_dir, exist_ok = True)
    
    child_files = get_centromere_files(child)
    parent1_files = get_centromere_files(parent1)
    parent2_files = get_centromere_files(parent2)
    
    for chrom in child_files:
        for parent, parent_files in zip([parent1, parent2], [parent1_files, parent2_files]):
            if chrom in parent_files:
                for child_hap in child_files[chrom]:
                    for parent_hap in parent_files[chrom]:
                        print("# Runner # Aligning {} from child {} haplotype {} to parent {} haplotype {}".format(chrom, child, child_hap, parent, parent_hap), file = sys.stderr)
                        child_file = child_files[chrom][child_hap]
                        parent_file = parent_files[chrom][parent_hap]
                        
                        aln_output_dir = os.path.join(output_dir, "{}_{}_{}_{}_{}".format(child, child_hap, parent, parent_hap, chrom))
                        
                        subprocess.check_call(" ".join(["tandem_aligner", "--first", child_file,
                                               "--second", parent_file, "--output-dir", aln_output_dir]), shell = True)
                        print("# Runner # Completed alignment", file = sys.stderr)
                    
    
    
    
    