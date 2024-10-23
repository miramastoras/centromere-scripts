#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:26:09 2024

Use annotations from HumAS-HMMER to identify centromeres in an assembly
Mappings from the assembly to, e.g., CHM13 can be used to identify the chromosome
for cross-chromosome shared HORs. Outputs a BED file of the HOR arrays

@author: Jordan
"""

import subprocess
import sys
import re
import os
import argparse
import collections

# bin contigs according to which chromosome the preponderance of their mapped bases
# go to (not all contigs will get assigned)
def assign_contigs_to_chroms(paf_fp, min_prop):
    
    contig_chrom_cov = {}
    
    with open(paf_fp) as f:
        for line in f:
            tokens = line.strip().split("\t")
            seq_len = int(tokens[3]) - int(tokens[2])
            contig = tokens[0]
            chrom = tokens[5]
            
            if contig not in contig_chrom_cov:
                contig_chrom_cov[contig] = {}
            
            chrom_cov = contig_chrom_cov[contig]
            if chrom not in chrom_cov:
                chrom_cov[chrom] = 0
                
            chrom_cov[chrom] += seq_len
    
    assignments = {}
    
    debug = False
    
    for contig in sorted(contig_chrom_cov):
        total = 0
        max_chrom = None
        max_len = 0
        for chrom in contig_chrom_cov[contig]:
            chrom_len = contig_chrom_cov[contig][chrom]
            if chrom_len > max_len:
                max_len = chrom_len
                max_chrom = chrom
            total += chrom_len
        
        
        if max_len / total > min_prop:
            assignments[contig] = max_chrom
            
        
        if debug:
            print("contig: {}, assignment: {}, coverages:".format(contig, max_chrom), file = sys.stderr)
            for chrom in sorted(contig_chrom_cov[contig]):
                print("\t{}: {}, {:.5f}".format(chrom, contig_chrom_cov[contig][chrom], contig_chrom_cov[contig][chrom] / total), file = sys.stderr)


    return assignments
    

# identify the dominant censat types for each contig (if there is one) and merge
# all of the arrays of the same type into one interval
def find_active_arrays(censat_fp, min_prop):
    
    contig_arrays = {}
    
    with open(censat_fp) as f:
        for line in f:
            if line.startswith("track name"):
                continue
            tokens = line.strip().split()
            if not tokens[3].startswith("active_hor"):
                continue
            contig = tokens[0]
            begin = int(tokens[1])
            end = int(tokens[2])
            
            m = re.search("S[0-9/]+C([XYM0-9/]+)H\d+L", tokens[3])
            if m is None:
                print(tokens[3], file = sys.stderr)
            assert(m is not None)
            chrom_set = tuple(m.group(1).split("/"))
            
            if contig not in contig_arrays:
                contig_arrays[contig] = {}
            
            arrays = contig_arrays[contig]
            
            if chrom_set not in arrays:
                arrays[chrom_set] = []
                
            arrays[chrom_set].append((begin, end))
    
    debug = False
    if debug:
        for contig in sorted(contig_arrays):
            print(contig, file = sys.stderr)
            for hor_type in sorted(contig_arrays[contig]):
                print("\t" + ", ".join(str(t) for t in hor_type), file = sys.stderr)
                for b, e in contig_arrays[contig][hor_type]:
                    print("\t\t{}\t{}".format(b, e), file = sys.stderr)
    
    
    merged_arrays = []
    for contig in sorted(contig_arrays):
        arrays = contig_arrays[contig]
        total_len = 0
        max_type = None
        max_type_len = 0
        for hor_type in arrays:
            type_len = sum(e - b for b, e in arrays[hor_type])
            total_len += type_len
            if type_len > max_type_len:
                max_type_len = type_len
                max_type = hor_type
        
        if max_type_len / total_len > min_prop:
            merged_arrays.append((contig, min(b for b,e, in arrays[max_type]), 
                                  max(e for b,e, in arrays[max_type]), max_type))
        else:
            print("filter arrays on contig {} because most prevalent HOR type {} only accounts for proportion {}".format(contig, ",".join(str(v) for v in max_type), max_type_len / total_len), file = sys.stderr)
            
    
    return merged_arrays
            


# returned as (contig, begin, end, chrom) tuples
def reconcile_chrom_assignments(merged_arrays, chrom_assignments, contig_lens, shoulder):
    
    hor_types_to_arrays = {}
    for i in range(len(merged_arrays)):
        hor_type = merged_arrays[i][3]
        if hor_type not in hor_types_to_arrays:
            hor_types_to_arrays[hor_type] = []
        hor_types_to_arrays[hor_type].append(i)
    
    reconciled = []
    
    for hor_type in hor_types_to_arrays:
        if len(hor_type) == 1:
            # there is only one chromosome that this type is associated with
            if len(hor_types_to_arrays[hor_type]) == 1:
                i = hor_types_to_arrays[hor_type][0]
                chrom = "chr" + str(merged_arrays[i][3][0])
                contig = merged_arrays[i][0]
                if contig in chrom_assignments:
                    if chrom != chrom_assignments[contig]:
                        print("filtering for discordant HOR type " + chrom + " and contig chrom " + chrom_assignments[contig], file = sys.stderr)
                        print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
                        continue
                if merged_arrays[i][1] < shoulder or merged_arrays[i][2] + shoulder > contig_lens[contig]:
                    print("filtering for being at end of contig of length {} in unique type".format(contig_lens[contig]), file = sys.stderr)
                    print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
                    continue
                reconciled.append((contig, merged_arrays[i][1], merged_arrays[i][2], chrom))
            else:
                print("filtering for discontiguous array:", file = sys.stderr)
                for i in hor_types_to_arrays[hor_type]:
                    print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
            
        else:
            # we have to disambiguate using the mappings
            assignment_counts = {}
            for i in hor_types_to_arrays[hor_type]:
                contig = merged_arrays[i][0]
                if contig in chrom_assignments:
                    chrom = chrom_assignments[contig]
                    if chrom not in assignment_counts:
                        assignment_counts[chrom] = 1
                    else:
                        assignment_counts[chrom] += 1
            
            if sum(assignment_counts.values()) != len(hor_types_to_arrays[hor_type]):
                print("filtering HOR group for incomplete assignments:", file = sys.stderr)
                for i in hor_types_to_arrays[hor_type]:
                    print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
            else:
                # we were able to assign all of them, which means that we can rule out one
                # secretly being a fragmented part of the same array as another (but they can
                # still be fragmented non-secretly)
                
                for i in hor_types_to_arrays[hor_type]:
                    contig = merged_arrays[i][0]
                    chrom = chrom_assignments[contig]
                    if assignment_counts[chrom] == 1:
                        if chrom in set("chr" + str(v) for v in hor_type):
                            if merged_arrays[i][1] >= shoulder and merged_arrays[i][2] + shoulder <= contig_lens[contig]:
                                reconciled.append((contig, merged_arrays[i][1], merged_arrays[i][2], chrom))
                            else:
                                print("filtering for being at end of contig of length {} in a resolved group".format(contig_lens[contig]), file = sys.stderr)
                                print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
                        else:
                            print("filtering for discordant mapped chromosome " + chrom + ":", file = sys.stderr)
                            print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
                    else:
                        print("filtering for discontiguous array over chromosome " + chrom + " in a resolved HOR group:", file = sys.stderr)
                        print("\t" + "\t".join(str(v) for v in merged_arrays[i]), file = sys.stderr)
    
    return reconciled
                

# returns list of bools, with True if on reverse strand
def determine_strand(hor_arrays, as_hor_sf_fp, flank_size):
    
    contig_to_array = {}
    for i in range(len(hor_arrays)):
        contig_to_array[hor_arrays[i][0]] = i
        
    strand_counts = [[0, 0] for array in hor_arrays]
    
    with open(as_hor_sf_fp) as f:
        for line in f:
            if line.startswith("track name"):
                continue
            tokens = line.strip().split()
            contig = tokens[0]
            if contig in contig_to_array:
                i = contig_to_array[contig]
                mono_begin = int(tokens[1])
                mono_end = int(tokens[2])
                arr_begin = hor_arrays[i][1]
                arr_end = hor_arrays[i][2]
                
                if ((mono_begin >= arr_begin and mono_end <= arr_begin + flank_size) or 
                    (mono_begin >= arr_end - flank_size and mono_end <= arr_end)):
                    
                    # this monomer is in the shoulder of the array (and therefore unlikely to be
                    # affected by inversions)
                    
                    strand = tokens[5]
                    if strand == '+':
                        strand_counts[i][0] += 1
                    elif strand == '-':
                        strand_counts[i][1] += 1
        
    
    return [fwd < rev for fwd, rev in strand_counts]
    
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

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--censat", type=str, required = True,
                        help="cenSat annotations in BED format")
    parser.add_argument("-a", "--as_hor_sf", type=str, required = True,
                        help="alpha satellite annotations from Hum-AS_HMMER in BED format")
    parser.add_argument("-p", "--paf", type=str, required = True,
                        help="alignment of the target assembly to a haploid T2T assembly in PAF format")
    parser.add_argument("-f", "--fasta", type=str, required = True,
                        help="assembly in fasta format")
    parser.add_argument("-w", "--strand_window", type=int, default = 10000,
                        help="determine strand the arrays based on the orientation of this much sequence on each of its end")
    parser.add_argument("-m", "--min_hor_prop", type=float, default = 0.9,
                        help="minimum proportion of chromosomal HOR types to call an array for that chromosome")
    parser.add_argument("-r", "--min_chrom_prop", type=float, default = 0.7,
                        help="minimum proportion of mapped bases to assign a contig to a chromosome (for cross-chromosomal HORs)")
    parser.add_argument("-s", "--min_shoulder", type=int, default = 20000,
                        help="minimum distance from array to ends of contig (to ensure unbroken)")
    
    args = parser.parse_args()
    
    contig_lens = get_contig_lens(args.fasta)
    
    contig_chrom_assignments = assign_contigs_to_chroms(args.paf, args.min_chrom_prop)
    active_arrays = find_active_arrays(args.censat, args.min_hor_prop)
    
    unambiguous_active_arrays = reconcile_chrom_assignments(active_arrays, contig_chrom_assignments, contig_lens, args.min_shoulder)
    strands = determine_strand(unambiguous_active_arrays, args.as_hor_sf, args.strand_window)
    
    for i in range(len(unambiguous_active_arrays)):
        contig, begin, end, chrom = unambiguous_active_arrays[i]
        strand = "+"
        if strands[i]:
            strand = "-"
        
        print("\t".join(str(v) for v in [contig, begin, end, chrom, ".", strand]))
        
        
        
        
        
