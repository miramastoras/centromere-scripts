#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:31:44 2024

@author: Jordan
"""

import sys
import skbio
import collections

# hack an enum for a subtree's phylogenetic consistency status
UNOBSERVED = 0
HOMOGENOUS0 = 1
HOMOGENOUS1 = 2
DISJOINT = 3
CONTAINED0 = 4
CONTAINED1 = 5
DISCORDANT = 6

harmonic_sum_memo = [float("nan"), 1.0]
def harmonic_number(n):
    while len(harmonic_sum_memo) <= n:
        harmonic_sum_memo.append(harmonic_sum_memo[-1] + 1.0 / len(harmonic_sum_memo))
    return harmonic_sum_memo[n]

def is_phylo_consistent(tree, sample_to_row, column):
    
    cons_states = [UNOBSERVED for i in range(tree.count())]
    
    for node in tree.postorder():
        par_id = None
        if node.parent:
            par_id = node.parent.id
        # print("at node {}, name {}, parent {}".format(node.id, node.name, par_id))
        if node.is_tip():
            if node.name not in sample_to_row:
                print("sample missing from matrix: {}".format(node.name), file = sys.stderr)
                sys.exit(1)
            i = sample_to_row[node.name]
            # print("leaf allele is {}".format(column[i]))
            if column[i] != None:
                if column[i]:
                    cons_states[node.id] = HOMOGENOUS1
                else:
                    cons_states[node.id] = HOMOGENOUS0
        else:
            child_states = collections.Counter(cons_states[child.id] for child in node)
            if len(child_states) == 1 and UNOBSERVED in child_states:
                # still no state observed
                cons_states[node.id] = UNOBSERVED
            else:
                # remove unobserved 
                del child_states[UNOBSERVED]
                if len(child_states) == 1:
                    state = next(iter(child_states))
                    if state == HOMOGENOUS0 or state == HOMOGENOUS1 or child_states[state] == 1:  
                        cons_states[node.id] = state
                    else:
                        cons_states[node.id] = DISCORDANT
                elif len(child_states) > 2:
                    # any combination of 3 observed states is inconsistent
                    cons_states[node.id] = DISCORDANT
                else:
                    if HOMOGENOUS0 in child_states and HOMOGENOUS1 in child_states:
                        # can be organized into two subtrees with distinct alleles
                        cons_states[node.id] = DISJOINT
                    elif DISJOINT in child_states and HOMOGENOUS0 in child_states and child_states[DISJOINT] == 1:
                        # can be organized into subtree of allele 1 appearing within allele 0 background
                        cons_states[node.id] = CONTAINED1
                    elif DISJOINT in child_states and HOMOGENOUS1 in child_states and child_states[DISJOINT] == 1:
                        # can be organized into subtree of allele 0 appearing within allele 1 background
                        cons_states[node.id] = CONTAINED0
                    elif CONTAINED0 in child_states and HOMOGENOUS1 in child_states and child_states[CONTAINED0] == 1:
                        # can still be organized into subtree of allele 1 appearing within allele 0 background
                        cons_states[node.id] = CONTAINED0
                    elif CONTAINED1 in child_states and HOMOGENOUS0 in child_states and child_states[CONTAINED1] == 1:
                        # can still be organized into subtree of allele 0 appearing within allele 1 background
                        cons_states[node.id] = CONTAINED1
                    else:
                        cons_states[node.id] = DISCORDANT
        # print("choose state {}".format(["UNOBSERVED", "HOMOGENOUS0", "HOMOGENOUS1", "DISJOINT", "CONTAINED0", "CONTAINED1", "DISCORDANT"][cons_states[node.id]]))
    
    return cons_states[tree.id] != DISCORDANT
    

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./var_phylo_consistency.py tree.nwk var_mat.tsv", file = sys.stderr)
        sys.exit(1)
    
    tree_fp = sys.argv[1]
    mat_fp = sys.argv[2]
    
    tree = skbio.TreeNode.read(open(tree_fp, "r"))
    tree.assign_ids()
    
    mat = []
    sample_to_row = {}
    allele_0 = []
    for line in open(mat_fp, "r"):
        tokens = line.strip().split()
        sample_to_row[tokens[0]] = len(mat)
        row = []
        for i in range(1, len(tokens)):
            if len(mat) == 0:
                allele_0.append(None)
            token = tokens[i]
            allele = []
            if token != "?":
                allele = []
                for subtoken in token.split(","):
                    if allele_0[i - 1] is None:
                        allele_0[i - 1] = subtoken
                
                    allele.append(subtoken != allele_0[i - 1])
            row.append(allele)            
        
        mat.append(row)
    
    
    num_consistent = 0
    num_possibly_inconsistent = 0
    num_singleton = 0
    expected_num_singletons = 0.0
    for j in range(len(mat[0])):
        column = []
        num_0 = 0
        num_1 = 0
        for i in range(len(mat)):
            if len(mat[i][j]) != 1:
                column.append(None)
            else:
                column.append(mat[i][j][0])
            for allele in mat[i][j]:
                if allele:
                    num_1 += 1
                else:
                    num_0 += 1
        
        expected_num_singletons += 1.0 / harmonic_number(num_0 + num_1 - 1)
        
        if num_0 == 1 or num_1 == 1:
            num_singleton += 1
        else:
            if any(allele is not None for allele in column):
                num_possibly_inconsistent += 1
                if is_phylo_consistent(tree, sample_to_row, column):
                    num_consistent += 1
    

    print("Expected number of singleton sites: {}".format(expected_num_singletons))
    print("{} of {} sites ({}%) are singletons".format(num_singleton, len(mat[0]), 100.0 * num_singleton / len(mat[0])))
    print("{} of {} possibly inconsistent sites ({}%) are phylogenetically consistent".format(num_consistent, num_possibly_inconsistent, 
                                                                                              100.0 * num_consistent / num_possibly_inconsistent))
    
    
    
    
    
    
    
    
    
    