#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:57:04 2024

Produce a Bandage-readable CSV for coloring the path(s) of a GFA file

@author: Jordan
"""

import sys

hex_alphabet = "0123456789abcdef"
from_hex = {hex_alphabet[i]:i for i in range(len(hex_alphabet))}

def to_hex(val):
    assert(val <= 255 and val >= 0)
    return hex_alphabet[val // 16] + hex_alphabet[val % 16]

def hex_color(r, g, b):
    return "#" + to_hex(r) + to_hex(g) + to_hex(b)

def color_interp(dist, length, col):
    
    base = 30
    
    interp = round((255 - base) * dist / length)
    
    if col == "r":
        return hex_color(base + interp, interp, interp)
    elif col == "g":
        return hex_color(interp, base + interp, interp)
    else:
        return hex_color(interp, interp, base + interp)

def hex_to_int(hex_code):
    p = len(hex_code) - 1
    v = 0
    place = 1
    while p >= 0:
        v += place * from_hex[hex_code[p]]
        p -= 1
        place *= 16
    return v
    

def parse_hex(color):
    assert(len(color) == 7)
    return hex_to_int(color[1:3]), hex_to_int(color[3:5]), hex_to_int(color[5:7])

def gradient(dist, length, grad_colors):
    
    assert(len(grad_colors) >= 2)
    
    pos = dist * (len(grad_colors) - 1) / length
    unit = int(pos)
    blend = pos - unit
    if unit == len(grad_colors) - 1:
        return grad_colors[-1]
    
    r1, g1, b1 = parse_hex(grad_colors[unit])
    r2, g2, b2 = parse_hex(grad_colors[unit + 1])
    
    r = round(r1 * (1.0 - blend) + r2 * blend)
    g = round(g1 * (1.0 - blend) + g2 * blend)
    b = round(b1 * (1.0 - blend) + b2 * blend)
    
    return hex_color(r, g, b)
    
    

if __name__ == "__main__":
    
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("usage: ./color_gfa_paths.py graph.gfa [which_path] > colors.csv")
        exit(1)
    
    gfa_fp = sys.argv[1]
    only_path = ""
    if len(sys.argv) >= 3:
        only_path = sys.argv[2]
    
    node_lengths = {}
    node_depth = {}
    restricted_node_depth = {}
    paths = []
    
    for line in open(gfa_fp):
        
        if line.startswith("S"):
            
            _, node_id, seq = line.strip().split()
            node_depth[node_id] = 0
            restricted_node_depth[node_id] = 0
            
            node_lengths[node_id] = len(seq)
        elif line.startswith("P"):
            _, name, steps, __ = line.strip().split()
            path = []
            for step in steps.split(","):
                node = step[:-1]
                if name == only_path:
                    restricted_node_depth[node] += 1
                node_depth[node] += 1
                path.append(node)
            paths.append((name, path))
    
    
    if len(only_path) != 0 and only_path not in [name for name, path in paths]:
        print("path " + only_path + " is not found in GFA", file = sys.stderr)
        exit(1)
        
    #grad_colors = ["#841e62", "#f3aa20", "#346b6d"]
    grad_colors = ["#6b3ec7", "#ea1a7f", "#fec603", "#a8f387", "#16d6fa"]
    null_color = "#a0a0a0"
    
    colored = set()
    print("Name,Color,Depth")
    for i in range(len(paths)):
        path_name, path = paths[i]
        if len(only_path) != 0 and path_name != only_path:
            continue
        col = "rgb"[i % 3]
        path_len = sum(node_lengths[node] for node in path)
        prefix_len = 0
        for node in path:
            if node not in colored:
                midpoint = round(prefix_len + node_lengths[node] / 2)
                if only_path:
                    node_col = gradient(midpoint, path_len, grad_colors)
                else:
                    node_col = color_interp(midpoint, path_len, col)
                if len(only_path) == 0:
                    print("{},{},{}".format(node, node_col, node_depth[node]))
                else:
                    print("{},{},\"{},{}\"".format(node, node_col, restricted_node_depth[node], node_depth[node]))
                colored.add(node)
            prefix_len += node_lengths[node]

    for node in node_depth:
        if node not in colored:
            if len(only_path) == 0:
                print("{},{},{}".format(node, null_color, node_depth[node]))
            else:
                print("{},{},\"{},{}\"".format(node, null_color, restricted_node_depth[node], node_depth[node]))
        
        
        
