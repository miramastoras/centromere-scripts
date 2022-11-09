#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 17:52:39 2022

@author: Jordan
"""


import sys
import pandas as pd

inf = sys.maxsize

if __name__ == "__main__":
    
    mutation_table = sys.argv[1]
    
    tab = pd.read_table(mutation_table, header = 0)
    
    skip_substitutions = False
    
    radius = 15000
    
    columns = ["chrom", "child_hap", "child_len", "child_pos_lower", "child_pos_upper", 
               "parent", "parent_hap", "parent_len", "parent_pos_lower", "parent_pos_upper",
               "type", "length"]
    print("\t".join(columns))
    
    
    i = 0
    while i < tab.shape[0]:
        j = i + 1
        while (j < tab.shape[0] and 
               tab.chrom.values[j] == tab.chrom.values[i] and 
               tab.child_hap.values[j] == tab.child_hap.values[i] and
               (abs(tab.child_pos.values[j] - tab.child_pos.values[j - 1]) < radius or
                abs(tab.parent_pos.values[j] - tab.parent_pos.values[j - 1]) < radius)):
            j += 1
        
        net_diff = 0
        min_child_pos = inf
        max_child_pos = -inf
        min_parent_pos = inf
        max_parent_pos = -inf
        
        for k in range(i, j):
            
            min_child_pos = min(min_child_pos, tab.child_pos.values[k])
            min_parent_pos = min(min_parent_pos, tab.parent_pos.values[k])
            
            if tab.type.values[k] == "INS":
                net_diff += tab.length.values[k]
                max_child_pos = max(max_child_pos, tab.child_pos.values[k] + tab.length.values[k])
                max_parent_pos = max(max_parent_pos, tab.parent_pos.values[k])
                
            elif tab.type.values[k] == "DEL":
                net_diff -= tab.length.values[k]
                max_child_pos = max(max_child_pos, tab.child_pos.values[k])
                max_parent_pos = max(max_parent_pos, tab.parent_pos.values[k] + tab.length.values[k])
            
            
        if net_diff != 0 or not skip_substitutions:
            if net_diff < 0:
                mut_type = "DEL"
            elif net_diff > 0:
                mut_type = "INS"
            else:
                mut_type = "SUB"
            print("\t".join(str(v) for v in [tab.chrom.values[i], tab.child_hap.values[i], tab.child_len.values[i], 
                                             min_child_pos, max_child_pos, tab.parent.values[i], tab.parent_hap.values[i], 
                                             tab.parent_len.values[i], min_parent_pos, max_parent_pos, mut_type, abs(net_diff)]))
                
        i = j
    