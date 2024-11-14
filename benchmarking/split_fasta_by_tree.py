#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 09:40:22 2024

@author: Jordan
"""
import sys
import skbio
import os
import argparse
import heapq


# convert length 0 edges into polytomies
def polytomize(tree):
    
    assert(tree is not None)
    
    queue = [tree]
    
    while queue:
        node = queue.pop()
        queue.extend(node.children)
        if node.length == 0.0 and not node.is_tip() and node is not tree:
            node.parent.children.remove(node)
            for child in node.children:
                child.parent = node.parent
                node.parent.children.append(child)
            tree.remove(node)
    tree.prune()
    tree.assign_ids()

# return root of subtrees that divide into most divergent groups with at least
# the target number of subtrees
def choose_subtrees(tree, target_num):
    assert(target_num >= 1)
    assert(type(target_num) == int)
    
    chosen = []
    heap = [(0.0, tree.id)]
    while len(heap) + len(chosen) < target_num and len(heap) != 0:
        
        # gather equidistant nodes
        nodes = []
        length, node_id = heapq.heappop(heap)
        nodes.append(tree.find_by_id(node_id))
        while len(heap) != 0 and heap[0][0] == length:
            nodes.append(tree.find_by_id(heapq.heappop(heap)[0]))
        
        for node in nodes:
            if node.is_tip():
                # we have to keep this one, so we do it as an isolated node
                chosen.append(node)
            else:
                # queue up children
                for child in node.children:
                    heapq.heappush(heap, (length + child.length, child.id))
    
    for length, node_id in heap:
        chosen.append(tree.find_by_id(node_id))
    
    return chosen


def parse_fasta(fa):
    parsed = []
    seq = ""
    name = ""
    annotations = ""
    with open(fa) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            line = line.strip()
            if line.startswith(">"):
                if len(seq) != 0 or len(name) != 0:
                    parsed.append((name, annotations, seq))
                divider = next((i for i, char in enumerate(line) if char.isspace()), len(line))
                name = line[1:divider]
                if divider < len(line):
                    annotations = line[divider + 1:]
                else:
                    annotations = ""
                seq = ""
            else:
                seq += line
    if len(seq) != 0 or len(name) != 0:
        parsed.append((name, annotations, seq))
    return parsed

# write a single fasta sequence to an open file
def write_fasta_seq(f, name, annotations, seq, line_len = 80):
    
    header = ">" + name
    if len(annotations) != 0:
        header += " " + annotations
    
    print(header, file = f)
    for i in range(0, len(seq), line_len):
        print(seq[i:i + line_len], file = f)

# write the subtree Newicks and FASTAs
def write_output(output_dir, subtrees, fasta_contents):
    
    name_to_seq = {record[0]: i for i, record in enumerate(fasta_contents)}
    
    for i, subtree in enumerate(subtrees):
    
        with open(os.path.join(output_dir, "subgroup_{}_seqs.fasta".format(i)), "w") as f:
            for node in subtree.tips():
                if node.name not in name_to_seq:
                    continue
                name, annotations, seq = fasta_contents[name_to_seq[node.name]]
                write_fasta_seq(f, name, annotations, seq)
        
        with open(os.path.join(output_dir, "subgroup_{}_tree.nwk".format(i)), "w") as f:
            print(subtree, file = f)
    

if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", type=str, required = True,
                        help="Tree to divide sequences by in Newick format")
    parser.add_argument("-f", "--fasta", type=str, required = True,
                        help="Sequences to divide in FASTA format")
    parser.add_argument("-n", "--num_subgroups", type=int, required = True,
                        help="The minimum number of subgroups to split into (true number may be higher if tree is polytomic)")
    parser.add_argument("-o", "--output_dir", type=str, required = True,
                        help="Directory to place output into (cannot already exist)")
    
    
    args = parser.parse_args()
    
    if os.path.exists(args.output_dir):
        print("error, output directory {} already exists".format(args.output_dir))
        sys.exit(1)
        
    os.mkdir(args.output_dir)
    
    tree = skbio.TreeNode.read(open(args.tree, "r"))
    
    fasta_contents = parse_fasta(args.fasta)
    
    tree = tree.shear([name for name, annotations, seq in fasta_contents])
    
    polytomize(tree)
    
    subtrees = choose_subtrees(tree, args.num_subgroups)
    
    write_output(args.output_dir, subtrees, fasta_contents)
