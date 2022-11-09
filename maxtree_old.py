#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 09:51:46 2022

@author: Jordan
"""

class Node():
    def __init__(self, key, value, info = None):
        self.key = key
        self.value = value
        self.info = info
        self.parent = None
        self.left = None
        self.right = None
        self.color = 1
        self.left_max = None
        self.right_max = None
        

    def subtree_max(self):
        m = self.value
        if self.left_max is not None:
            m = max(m, self.left_max)
        if self.right_max is not None:
            m = max(m, self.right_max)
        return m
    
    def set_left_max(self):
        if self.left is None:
            self.left_max = None
        else:
            self.left_max = self.left.subtree_max()
            
    def set_right_max(self):
        if self.right is None:
            self.right_max = None
        else:
            self.right_max = self.right.subtree_max()


class MaxTree():
    
    def __init__(self):
        self.root = None
        
    def size(self):
        return self.size_helper(self.root)
        
    def size_helper(self, node):
        if node is None:
            return 0
        return self.size_helper(node.left) + self.size_helper(node.right) + 1

    # Search the tree
    def search_tree_helper(self, node, key):
        if node is None or key == node.key:
            return node

        if key < node.key:
            return self.search_tree_helper(node.left, key)
        else:
            return self.search_tree_helper(node.right, key)

    # Balancing the tree after deletion
    def delete_fix(self, x):
        if x is None:
            return
        print("delete fix x {}, parent {}".format(x, x.parent))
        while x != self.root and x.color == 0:
            print("\tx {}, parent {}".format(x, x.parent))
            if x.parent is not None:
                if x == x.parent.left:
                    s = x.parent.right
                    if s is not None:
                        if s.color == 1:
                            s.color = 0
                            x.parent.color = 1
                            self.left_rotate(x.parent)
                            s = x.parent.right
                            
                        if s.left is not None and s.right is not None:
                            if s.left.color == 0 and s.right.color == 0:
                                s.color = 1
                                x = x.parent
                            else:
                                if s.right.color == 0:
                                    s.left.color = 0
                                    s.color = 1
                                    self.right_rotate(s)
                                    s = x.parent.right
            
                                s.color = x.parent.color
                                x.parent.color = 0
                                s.right.color = 0
                                self.left_rotate(x.parent)
                                x = self.root
                else:
                    s = x.parent.left
                    if s is not None:
                        if s.color == 1:
                            s.color = 0
                            x.parent.color = 1
                            self.right_rotate(x.parent)
                            s = x.parent.left
        
                        if s.left is not None and s.right is not None:
                            if s.left.color == 0 and s.right.color == 0:
                                s.color = 1
                                x = x.parent
                            else:
                                if s.left.color == 0:
                                    s.right.color = 0
                                    s.color = 1
                                    self.left_rotate(s)
                                    s = x.parent.left
            
                                s.color = x.parent.color
                                x.parent.color = 0
                                s.left.color = 0
                                self.right_rotate(x.parent)
                                x = self.root
        x.color = 0

    # replace u with v
    def __rb_transplant(self, u, v):
        if u.parent == None:
            self.root = v
        elif u == u.parent.left:
            u.parent.left = v
            u.parent.set_left_max()
        else:
            u.parent.right = v
            u.parent.set_right_max()
        
        if v is not None:
            v.parent = u.parent

    # Node deletion
    def delete_node(self, z):

        y = z
        y_original_color = y.color
        if z.left is None:
            x = z.right
            self.__rb_transplant(z, z.right)
        elif z.right is None:
            x = z.left
            self.__rb_transplant(z, z.left)
        else:
            y = self.minimum(z.right)
            y_original_color = y.color
            x = y.right
            if y.parent == z:
                if x is not None:
                    x.parent = y
            else:
                self.__rb_transplant(y, y.right)
                y.right = z.right
                y.right.parent = y
                y.set_right_max()

            self.__rb_transplant(z, y)
            y.left = z.left
            y.left.parent = y
            y.color = z.color
            y.set_left_max()
            
        if y_original_color == 0:
            self.delete_fix(x)

    # Balance the tree after insertion
    def fix_insert(self, k):
        while k.parent is not None and k.parent.color == 1:
            if k.parent.parent is not None:
                if k.parent == k.parent.parent.right:
                    u = k.parent.parent.left
                    if u is not None and u.color == 1:
                        u.color = 0
                        k.parent.color = 0
                        k.parent.parent.color = 1
                        k = k.parent.parent
                    else:
                        if k == k.parent.left:
                            k = k.parent
                            self.right_rotate(k)
                        k.parent.color = 0
                        k.parent.parent.color = 1
                        self.left_rotate(k.parent.parent)
                else:
                    u = k.parent.parent.right
    
                    if u is not None and u.color == 1:
                        u.color = 0
                        k.parent.color = 0
                        k.parent.parent.color = 1
                        k = k.parent.parent
                    else:
                        if k == k.parent.right:
                            k = k.parent
                            self.left_rotate(k)
                        k.parent.color = 0
                        k.parent.parent.color = 1
                        self.right_rotate(k.parent.parent)
                if k == self.root:
                    break
        self.root.color = 0

    def insert_node_helper(self, key, value, info = None):
        node = Node(key, value, info)

        y = None
        x = self.root

        while x is not None:
            y = x
            if node.key < x.key:
                if x.left_max is None or x.left_max < value:
                    x.left_max = value
                x = x.left
            else:
                if x.right_max is None or x.right_max < value:
                    x.right_max = value
                x = x.right

        node.parent = y
        if y == None:
            self.root = node
        elif node.key < y.key:
            y.left = node
        else:
            y.right = node

        if node.parent == None:
            node.color = 0
            return node

        if node.parent.parent == None:
            return node

        self.fix_insert(node)
        
        return node

    # Printing the tree
    def __print_helper(self, node, indent):
        if node is not None:
            if node.parent is not None:
                if node is node.parent.left:
                    left_indent = indent[:-1] + "     |"
                    right_indent = indent + "    |"
                else:
                    right_indent = indent[:-1] + "     |"
                    left_indent = indent + "    |"
            else:
                left_indent = "    |"
                right_indent = "    |"
                
            self.__print_helper(node.right, right_indent)
            l = node.left_max if node.left_max is not None else "."
            r = node.right_max if node.right_max is not None else "."
            print("{}----{}:{}({} l:{} r:{})".format(indent, node.key, node.value, node.color, l, r))
            self.__print_helper(node.left, left_indent)


    def left_rotate(self, x):
        y = x.right
        x.right = y.left
        if y.left is not None:
            y.left.parent = x

        y.parent = x.parent
        if x.parent == None:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y
        y.left = x
        x.parent = y
        
        
        x.right_max = y.left_max
        y.left_max = x.subtree_max()

    def right_rotate(self, x):
        y = x.left
        x.left = y.right
        if y.right is not None:
            y.right.parent = x

        y.parent = x.parent
        if x.parent == None:
            self.root = y
        elif x == x.parent.right:
            x.parent.right = y
        else:
            x.parent.left = y
        y.right = x
        x.parent = y
        
        x.left_max = y.right_max
        y.right_max = x.subtree_max()

    def print_tree(self):
        self.__print_helper(self.root, "")

    def minimum(self, node):
        while node.left is not None:
            node = node.left
        return node

    def maximum(self, node):
        while node.right is not None:
            node = node.right
        return node

    def successor(self, x):
        if x is None:
            return self.minimum(self.root)
        if x.right is not None:
            return self.minimum(x.right)

        y = x.parent
        while y is not None and x == y.right:
            x = y
            y = y.parent
        return y

    def predecessor(self,  x):
        if x is None:
            return self.maximum(self.root)
        if x.left is not None:
            return self.maximum(x.left)

        y = x.parent
        while y is not None and x == y.left:
            x = y
            y = y.parent

        return y

    def get_root(self):
        return self.root
    
    # the node with the greatest key value that is <= the query key
    def lower_bound(self, key):
        bound = None
        node = self.root
        while node != None:
            if node.key <= key and (bound is None or bound.key <= node.key):
                bound = node
            if key < node.key:
                node = node.left
            else:
                node = node.right
        return bound
    
    # the node with the least key value that is >= the query key
    def upper_bound(self, key):
        bound = None
        node = self.root
        while node != None:
            if node.key >= key and (bound is None or bound.key >= node.key):
                bound = node
            if key < node.key:
                node = node.left
            else:
                node = node.right
        return bound
    
    # The max value of any key higher than the node (not including the node itself)
    def ray_max_upper(self, node):
        m = node.right_max
        while node is not None:
            parent = node.parent
            if parent is not None:
                if node is parent.left:
                    if m is None or m < parent.value:
                        m = parent.value
                    if parent.right_max is not None and m < parent.right_max:
                        m = parent.right_max
            node = parent
        return m
    
    # The max value of any key lower than the node (not including the node itself)
    def ray_max_lower(self, node):
        m = node.left_max
        while node is not None:
            parent = node.parent
            if parent is not None:
                if node is parent.right:
                    if m is None or m < parent.value:
                        m = parent.value
                    if parent.left_max is not None and m < parent.left_max:
                        m = parent.left_max
            node = parent
        return m

    # return the node associated with a key or None if it is missing
    def find(self, key):
        return self.search_tree_helper(self.root, key)

    # delete the node associated with a key, raises error if key is not present
    def delete(self, key):
        node = self.find(key)
        if node is None:
            raise ValueError("Deleting key {} that is not in tree".format(key))
        self.delete_node(node)
        
    # add a new node
    def insert(self, key, value, info = None):
        return self.insert_node_helper(key, value, info)


if __name__ == "__main__":
    bst = MaxTree()

    print("insert 55")
    bst.insert(55, 4)
    bst.print_tree()
    print("insert 40")
    bst.insert(40, 5)
    bst.print_tree()
    print("insert 65")
    bst.insert(65, 6)
    bst.print_tree()
    print("insert 60")
    bst.insert(60, 7)
    bst.print_tree()
    print("insert 75")
    bst.insert(75, 8)
    bst.print_tree()
    print("insert 57")
    bst.insert(57, 9)
    bst.print_tree()


    print("\nAfter deleting an element")
    bst.delete(55)
    bst.print_tree()
    print()
    
    for i in range(38, 85, 7):
        ub = bst.upper_bound(i)
        lb = bst.lower_bound(i)
        if ub is not None:
            u = str(ub.key)
        else:
            u = "None"
        if lb is not None:
            l = str(lb.key)
        else:
            l = "None"
            
        print("{}: lb {}, ub {}".format(i, l, u))
        
    n = bst.minimum(bst.get_root())
    while n is not None:
        mu = bst.ray_max_upper(n)
        ml = bst.ray_max_lower(n)
        if mu is None:
            u = "None"
        else:
            u = str(mu)
        if ml is None:
            l = "None"
        else:
            l = str(ml)
            
        print("{}: lower max {}, upper max {}".format(n.key, l, u))
        
        n = bst.successor(n)
    