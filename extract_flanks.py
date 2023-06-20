#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 14:57:52 2022

@author: Jordan
"""

import sys
import tempfile
import subprocess
import re
import os
import time
import hashlib

if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print("usage: ./extract_flanks.py ref.fasta hors.bed flank_size > flanks.fasta", file = sys.stderr)
        sys.exit(1)
    
    fasta = sys.argv[1]
    bed = sys.argv[2]
    flank_size = int(sys.argv[3])
    
    # make the index
    if not os.path.exists(fasta + ".fai"):
        subprocess.check_call(["samtools", "faidx", fasta])
    
    # extract the sequence lengths
    chrom_length = {}
    for line in open(fasta + ".fai"):
        if type(line) == bytes:
            line = line.decode("utf-8")
        tokens = line.strip().split()
        chrom_length[tokens[0]] = int(tokens[1])
    
    # construct the region file
    region_by_endpoint = {}
    region_file_name = hashlib.sha256(bytes(str(time.time()), "utf-8")).hexdigest() + ".txt"
    region_file = open(region_file_name, "w")
    for line in open(bed):
        if type(line) == bytes:
            line = line.decode("utf-8")
        
        tokens = line.strip().split()
        chrom = tokens[0]
        b = int(tokens[1]) + 1
        e = int(tokens[2])
        region_by_endpoint[(chrom, b - 1)] = (b, e)
        region_by_endpoint[(chrom, e + 1)] = (b, e)
                        
        print("{}:{}-{}".format(chrom, max(1, b - flank_size), b - 1), file = region_file)
        print("{}:{}-{}".format(chrom, e + 1, min(e + flank_size, chrom_length[chrom])), file = region_file)                      
    region_file.close()
    
    temp_fasta_name = hashlib.sha256(bytes(str(time.time()), "utf-8")).hexdigest() + ".fasta"
    try:
        subprocess.check_call(f"samtools faidx -r {region_file_name} -o {temp_fasta_name} {fasta}", shell = True)
    except subprocess.CalledProcessError as e:
        print('exit code: {}'.format(e.returncode))
        if e.output is not None:
            print('stdout: {}'.format(e.output.decode(sys.getfilesystemencoding())))
        if e.stderr is not None:
            print('stderr: {}'.format(e.stderr.decode(sys.getfilesystemencoding())))
        raise e
    
    
    region_regex = "(\S+):(\d+)-(\d+)"
    for line in open(temp_fasta_name):
        if type(line) == bytes:
            line = line.decode("utf-8")
        line = line.strip()
        if line.startswith(">"):
            m = re.search(region_regex, line[1:])
            chrom = m.group(1)
            b = int(m.group(2))
            e = int(m.group(3))
            
            if (chrom, b) in region_by_endpoint:
                region_b, region_e = region_by_endpoint[(chrom, b)]
                direction = "R"
            else:
                region_b, region_e = region_by_endpoint[(chrom, e)]
                direction = "L"
            line += "({}:{}-{}{})".format(chrom, region_b, region_e, direction)
            
        print(line)
        
    os.remove(temp_fasta_name)
    os.remove(region_file_name)
        
    