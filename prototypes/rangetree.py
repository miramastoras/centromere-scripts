#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 22:01:13 2022

@author: Jordan
"""

import bisect
import sys
inf = sys.maxsize

class KeyWrapper:
    def __init__(self, iterable, key):
        self.it = iterable
        self.key = key

    def __getitem__(self, i):
        return self.key(self.it[i])

    def __len__(self):
        return len(self.it)

def bisect_key(v):
    if type(v) == int:
        return v
    else:
        return v.value[1]

class RangeVectorEntry:
    
    def __init__(self):
        self.value = None
        self.prefix_sum = None
        self.left_pos = None
        self.right_pos = None

class RangeTreeNode:
    def __init__(self):
        self.min = None
        self.max = None
        self.vector = None
        
    def __str__(self):
        return "xm: {}, xM: {}\n{}".format(self.min, self.max, ", ".join("(v: {}, ps: {}, l:{}, r:{})".format(e.value, e.prefix_sum, e.left_pos, e.right_pos) for e in self.vector) if self.vector is not None else None)

        
    def index_interval(self, x_begin, x_end, y_begin, y_end):
        
        wrapped = KeyWrapper(self.vector, bisect_key)
        
        l = bisect.bisect_left(wrapped, y_begin)
        r = bisect.bisect_left(wrapped, y_end)
        return (l, r)
                    
        

class RangeTree2D:
    
    # data comes in format (x, y, weight)
    def __init__(self, data, debug = False):
        
        self.debug = debug
        self.leaves_reached = 0
        
        # figure out the height of the internal nodes
        n = 1
        h = 0
        while n < len(data):
            h += 1
            n *= 2
            
        num_internal = n - 1
        
        self.tree = [RangeTreeNode() for i in range(num_internal + len(data))]
        
        
        # determine the x-sorted intervals for each node
        if debug:
            print("setting up x ranges in tree")
        
        x_sorted = sorted(data)
        
        for i in range(len(x_sorted)):
            self.tree[i + num_internal].max = x_sorted[i]
            self.tree[i + num_internal].min = x_sorted[i]
            
        for i in range(num_internal - 1, -1, -1):
            l = 2 * i + 1
            if l < len(self.tree):
                
                if self.tree[l].min is not None:
                    self.tree[i].min = self.tree[l].min
                
                r = l + 1
                if r < len(self.tree) and self.tree[r].max is not None:
                    self.tree[i].max = self.tree[r].max
                elif self.tree[l].max is not None:
                    self.tree[i].max = self.tree[l].max
        
        if debug:
            print("setting up y vectors")
        # determine the y-sorted subset and create fractional cascading vectors
        
        y_sorted = sorted(data, key = lambda xyw: (xyw[1], xyw[0], xyw[2]))
        
        if len(self.tree) != 0:
            self.__init__helper(y_sorted, 0)
        
    def __init__helper(self, y_sorted, i):
        
        node = self.tree[i]
        
        l = 2 * i + 1
        r = l + 1
        
        node.vector = [RangeVectorEntry() for j in range(len(y_sorted))]
        
        y_sorted_left = []
        y_sorted_right = []
        
        prefix_sum = 0
        
        for j in range(len(y_sorted)):
            
            # add an entry to the range vector
            
            node.vector[j].value = y_sorted[j][:2]
            
            prefix_sum += y_sorted[j][2]
            node.vector[j].prefix_sum = prefix_sum
            
            node.vector[j].left_pos = len(y_sorted_left)
            node.vector[j].right_pos = len(y_sorted_right)
            
            # divvy the points among the two children
            
            if l < len(self.tree) and self.tree[l].max is not None and y_sorted[j][:2] <= self.tree[l].max:
                y_sorted_left.append(y_sorted[j])
            elif r < len(self.tree) and self.tree[r].min is not None:
                y_sorted_right.append(y_sorted[j])
        
        # recurse into the children

        del y_sorted
        if l < len(self.tree):
            self.__init__helper(y_sorted_left, l)
        if r < len(self.tree):
            self.__init__helper(y_sorted_right, r)
            
        if self.debug and l >= len(self.tree) and r >= len(self.tree):
            self.leaves_reached += 1
            if self.leaves_reached % 100000 == 0:
                print("completed construction for {} of {} leaves".format(self.leaves_reached, len(self.tree[0].vector)))
    
    def range_query_helper(self, i, x_begin, x_end, idx_1, idx_2, points):
                
        if i >= len(self.tree) or idx_1 >= idx_2:
            # there is no leaf here the y interval is exhausted
            return
        
        node = self.tree[i]
        
        if node.max is None or node.max < (x_begin, -inf) or node.min >= (x_end, -inf):
            # this interval is disjoint with the query
            return
        
        if node.min >= (x_begin, -inf) and node.max < (x_end, -inf):
            for j in range(idx_1, idx_2):
                entry = node.vector[j]
                points.append((entry.value[0], entry.value[1]))
        else:
            idx_1_l = 0
            idx_1_r = 0
            idx_2_l = 0
            idx_2_r = 0
            if idx_1 < len(node.vector):
                idx_1_l = node.vector[idx_1].left_pos
                idx_1_r = node.vector[idx_1].right_pos
            else:
                if 2 * i + 1 < len(self.tree):
                    idx_1_l = len(self.tree[2 * i + 1].vector)
                if 2 * i + 2 < len(self.tree):
                    idx_1_r = len(self.tree[2 * i + 2].vector)
            if idx_2 < len(node.vector):
                idx_2_l = node.vector[idx_2].left_pos
                idx_2_r = node.vector[idx_2].right_pos
            else:
                if 2 * i + 1 < len(self.tree):
                    idx_2_l = len(self.tree[2 * i + 1].vector)
                if 2 * i + 2 < len(self.tree):
                    idx_2_r = len(self.tree[2 * i + 2].vector)
                    
            self.range_query_helper(2 * i + 1, x_begin, x_end, 
                                    idx_1_l, idx_2_l, points)
            self.range_query_helper(2 * i + 2, x_begin, x_end,
                                    idx_1_r, idx_2_r, points)
    
    def range_weight_query_helper(self, i, x_begin, x_end, idx_1, idx_2):
                
        if i >= len(self.tree) or idx_1 >= idx_2:
            # there is no leaf here the y interval is exhausted
            return 0
        
        node = self.tree[i]
        
        # TODO: i could probably improve these bounds by using y_begin/end instead of inf
        
        if node.max is None or node.max < (x_begin, -inf) or node.min >= (x_end, -inf):
            # this interval is disjoint with the query
            return 0
        
        weight_sum = 0
        
        if node.min >= (x_begin, -inf) and node.max < (x_end, -inf):
            
            if idx_2 > 0:
                weight_sum = node.vector[idx_2 - 1].prefix_sum
            if idx_1 > 0:
                weight_sum -= node.vector[idx_1 - 1].prefix_sum
            
        else:
            idx_1_l = 0
            idx_1_r = 0
            idx_2_l = 0
            idx_2_r = 0
            if idx_1 < len(node.vector):
                idx_1_l = node.vector[idx_1].left_pos
                idx_1_r = node.vector[idx_1].right_pos
            else:
                if 2 * i + 1 < len(self.tree):
                    idx_1_l = len(self.tree[2 * i + 1].vector)
                if 2 * i + 2 < len(self.tree):
                    idx_1_r = len(self.tree[2 * i + 2].vector)
            if idx_2 < len(node.vector):
                idx_2_l = node.vector[idx_2].left_pos
                idx_2_r = node.vector[idx_2].right_pos
            else:
                if 2 * i + 1 < len(self.tree):
                    idx_2_l = len(self.tree[2 * i + 1].vector)
                if 2 * i + 2 < len(self.tree):
                    idx_2_r = len(self.tree[2 * i + 2].vector)
                
                
            
            
            weight_sum = self.range_weight_query_helper(2 * i + 1, x_begin, x_end,
                                                        idx_1_l, idx_2_l)
            weight_sum += self.range_weight_query_helper(2 * i + 2, x_begin, x_end,
                                                        idx_1_r, idx_2_r)
            
        return weight_sum
        
        
        
    def range_query(self, x_begin, x_end, y_begin, y_end):
        
        points = []
                
        idx_begin, idx_end = self.tree[0].index_interval(x_begin, x_end, y_begin, y_end)
        
        self.range_query_helper(0, x_begin, x_end, idx_begin, idx_end, points)
        
        return points
    
    def range_weight_query(self, x_begin, x_end, y_begin, y_end):
                
        idx_begin, idx_end = self.tree[0].index_interval(x_begin, x_end, y_begin, y_end)
        
        
        return self.range_weight_query_helper(0, x_begin, x_end, idx_begin, idx_end)
         
        
        
if __name__ == "__main__":
    
    data = [(1, 5, 1), (1, 6, 2), (2, 5, 4), (2, 6, 8)]
    
    tree = RangeTree2D(data)
    
    for x1 in range(1, 4):
        for x2 in range(x1, 4):
            for y1 in range(5, 8):
                for y2 in range(y1, 8):
                    w = tree.range_weight_query(x1, x2, y1, y2)
                    p = tree.range_query(x1, x2, y1, y2)
                    print("x {}:{}, y {}:{}, weight {}, num {}".format(x1, x2, y1, y2, w, len(p)))
    
    





      
        
