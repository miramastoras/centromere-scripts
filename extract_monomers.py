#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:24:59 2022

@author: Jordan
"""

import sys
import collections
import os
import subprocess
import tempfile

if __name__ == "__main__":
    
    bed = sys.argv[1]
    fasta = sys.argv[2]
    out_dir = sys.argv[3]
    
    assert(not os.path.exists(out_dir))
    
        
    families = collections.defaultdict(list)
    
    for line in open(bed):
        if type(line) == bytes:
            line = line.decode("utf-8")
        tokens = line.strip().split()
        families[tokens[3]].append((tokens[0], tokens[1], tokens[2]))
        
    
    avg_fam_size = 0.0
    for fam in families:
        avg_fam_size += len(families[fam])
    avg_fam_size /= len(families)
    
    
    os.makedirs(out_dir)
    
    for fam in families:
        family = families[fam]
        if len(family) < 0.25 * avg_fam_size:
            continue
        
        region_file = tempfile.NamedTemporaryFile() 
        for chrom, begin, end in family:
            line = "{}:{}-{}\n".format(chrom, begin, end)
            line = bytes(line, "utf-8")
            region_file.write(line)
            
        fam_fasta = os.path.join(out_dir, fa.replace("/", "-") + ".fasta")
        subprocess.check_call(["samtools", "faidx", "-r", region_file.name, "-o", fam_fasta, fasta])
            
        