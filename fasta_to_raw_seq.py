#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:23:18 2022

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    saw_header = False
    seq = ""
    with open(sys.argv[1]) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            if line.startswith(">"):
                assert(not saw_header)
                saw_header = True
                continue
            seq += line.strip().upper()
            
    print(seq)