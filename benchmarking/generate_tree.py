#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:10:10 2023

Simulate a tree according to a coalescent model

@author: Jordan
"""

import sys
import msprime
import re

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./generate_tree.py num_samples expected_tmrca_gens > tree.nwk", file = sys.stderr)
        exit(1);
    
    num_samples = int(sys.argv[1])
    exp_gens = float(sys.argv[2])
    
    pop_size = round(0.5 * exp_gens / (1.0 - 1.0 / num_samples))
    
    for tree in msprime.sim_ancestry(samples = num_samples, population_size = pop_size, 
                                     sequence_length = 1, ploidy = 1).trees():
        
        raw_newick = tree.as_newick()
        
        # transform continuous values to integer
        newick = ""
        prev = 0
        for m in re.finditer(r":(\d+(\.\d+))", raw_newick):
            newick += raw_newick[prev : m.start() + 1]
            newick += str(round(float(m.group(1))))
            prev = m.end()
        newick += raw_newick[prev:]
        
        print(newick)
        
            
