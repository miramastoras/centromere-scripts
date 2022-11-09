#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:25:21 2022

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    f = None
    if len(sys.argv) == 1:
        f = sys.stdin
    elif len(sys.argv) == 2:
        f = open(sys.argv[1])
    else:
        print("error: invalid input", file = sys.stderr)
        sys.exit()
        
    for line in f:
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith(">"):
            print(line.strip())
        else:
            print(line.strip().upper())