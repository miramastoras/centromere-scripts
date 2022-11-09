#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 17:49:16 2022

@author: Jordan
"""

# Implementing Red-Black Tree in Python


import sys


# Node creation
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
            

class RedBlackTree():
    def __init__(self):
        self.root = None

    # Balancing the tree after deletion
    def delete_fix(self, x):
        while x is not None and x != self.root and x.color == 0:
            # x must have a parent because it is not the root
#            print("fixing at {}:{}".format(x.key, x.value))
            if x == x.parent.left:
                s = x.parent.right
                if s is not None:
                    if s.color == 1:
                        s.color = 0
                        x.parent.color = 1
                        self.left_rotate(x.parent)
                        s = x.parent.right
    
                if s is not None:
                    
                    if s.left is not None and s.right is not None and s.left.color == 0 and s.right.color == 0:
                        s.color = 1
                        x = x.parent
                    else:
                        if s.right is not None and s.right.color == 0:
                            if s.left is not None:
                                s.left.color = 0
                            s.color = 1
                            self.right_rotate(s)
                            s = x.parent.right
    
                        if s is not None and x.parent is not None:
                            s.color = x.parent.color
                        if x.parent is not None:
                            x.parent.color = 0
                        if s is not None and s.right is not None:
                            s.right.color = 0
                        self.left_rotate(x.parent)
                        x = self.root
                else:
                    x = x.parent
            else:
                s = x.parent.left
                if s is not None:
                    if s.color == 1:
                        s.color = 0
                        x.parent.color = 1
                        self.right_rotate(x.parent)
                        s = x.parent.left
                
                if s is not None:
                    if s.left is not None and s.right is not None and s.left.color == 0 and s.right.color == 0:
                        s.color = 1
                        x = x.parent
                    else:
                        if s.left is not None and s.left.color == 0:
                            if s.right is not None:
                                s.right.color = 0
                            s.color = 1
                            self.left_rotate(s)
                            s = x.parent.left
    
                        if s is not None and x.parent is not None:
                            s.color = x.parent.color
                        if x.parent is not None:
                            x.parent.color = 0
                        if s is not None and s.left is not None:
                            s.left.color = 0
                        self.right_rotate(x.parent)
                        x = self.root
                else:
                    x = x.parent
        if x is not None:
            x.color = 0

    def __rb_transplant(self, u, v):
        if u.parent is None:
            self.root = v
        elif u == u.parent.left:
            u.parent.left = v
            u.parent.set_left_max()
        else:
            u.parent.right = v
            u.parent.set_right_max()
        v.parent = u.parent

    # Node deletion
    def delete_node_helper(self, node, key):
        z = None
        while node is not None:
            if node.key == key:
                z = node

            if node.key <= key:
                node = node.right
            else:
                node = node.left

        if z is None:
            print("Cannot find key in the tree")
            return
        
#        print("find {}:{}".format(z.key, z.value))

#        print("before initial delete")
#        self.prettyPrint()
        y = z
        y_original_color = y.color
        if z.left is None and z.right is None:
#            print("case 0")
            x = None
            if z.parent is not None:
                if z == z.parent.left:
                    z.parent.left = None
                    z.parent.left_max = None
                else:
                    z.parent.right = None
                    z.parent.right_max = None
            else:
                self.root = None
        elif z.left is None and z.right is not None:
#            print("case 1")
            x = z.right
            self.__rb_transplant(z, z.right)
        elif z.right is None and z.left is not None:
#            print("case 2")
            x = z.left
            self.__rb_transplant(z, z.left)
        else:
#            print("case 3")
            y = self.minimum(z.right)
            y_original_color = y.color
            x = y.right
            
            # handle left of z
            y.left = z.left
            y.left_max = z.left_max
            if y.left is not None:
                y.left.parent = y
            
            # handle right of z (if necessary)
            if y.parent != z:
                
                if x is not None:
                    x.parent = y.parent
                    
                if y == y.parent.left:
                    y.parent.left = x
                    y.parent.set_left_max()
                else:
                    y.parent.right = x
                    y.parent.set_right_max()
                
                y.right = z.right
                if y.right is not None:
                    y.right.parent = y
                
            # handle parent of z
            y.parent = z.parent
            if y.parent is None:
                self.root = y
            else:
                if z == y.parent.left:
                    y.parent.left = y
                    y.parent.set_left_max()
                else:
                    y.parent.right = y
                    y.parent.set_right_max()
                
            y.color = z.color
            
#        print("after initial delete")
#        self.prettyPrint()
            
        if y_original_color == 0:
            self.delete_fix(x)

    # Balance the tree after insertion
    def fix_insert(self, k):
        while k.parent is not None and k.parent.color == 1 and k.parent.parent is not None:
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

    # Printing the tree
    def __print_helper(self, node, indent, last):
        if node is not None:
            sys.stdout.write(indent)
            if last:
                sys.stdout.write("R----")
                indent += "     "
            else:
                sys.stdout.write("L----")
                indent += "|    "

            s_color = "RED" if node.color == 1 else "BLACK"
            print(str(node.key) + "(" + s_color + ")")
            self.__print_helper(node.left, indent, False)

    def minimum(self, node):
        while node.left is not None:
            node = node.left
        return node

    def maximum(self, node):
        while node.right is not None:
            node = node.right
        return node

    def height(self):
        return self.height_helper(self.root)

    def height_helper(self, x):
        if x is None:
            return 0
        else:
            return 1 + max(self.height_helper(x.left), self.height_helper(x.right))

    def successor(self, x):
        if x is None:
            return None
        if x.right is not None:
            return self.minimum(x.right)

        s = None
        while x.parent is not None:
            if x.parent.left == x:
                s = x.parent
                break
            x = x.parent
        return s

    def predecessor(self,  x):
        if x is None:
            return None
        if x.left is not None:
            return self.maximum(x.left)

        s = None
        while x.parent is not None:
            if x.parent.right == x:
                s = x.parent
                break
            x = x.parent
        return s

    def left_rotate(self, x):
        if x is None:
            return
        y = x.right
        if y is None:
            return 
        x.right = y.left
        if y.left is not None:
            y.left.parent = x

        if y is not None:
            y.parent = x.parent
        if x.parent is None:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y
        y.left = x
        x.parent = y
        
        x.set_right_max()
        y.set_left_max()

    def right_rotate(self, x):
        if x is None:
            return
        y = x.left
        if y is None:
            return
        x.left = y.right
        if y.right is not None:
            y.right.parent = x

        if y is not None:
            y.parent = x.parent
        if x.parent is None:
            self.root = y
        elif x == x.parent.right:
            x.parent.right = y
        else:
            x.parent.left = y
        y.right = x
        x.parent = y
        
        x.set_left_max()
        y.set_right_max()

    def insert(self, key, value, info = None):
        
        node = Node(key, value, info)

        y = None
        x = self.root

        while x is not None:
            y = x
            if node.key < x.key:
                x = x.left
                if y.left_max is None:
                    y.left_max = node.value
                else:
                    y.left_max = max(y.left_max, node.value)
            else:
                x = x.right
                if y.right_max is None:
                    y.right_max = node.value
                else:
                    y.right_max = max(y.left_max, node.value)

        node.parent = y
        if y is None:
            self.root = node
        elif node.key < y.key:
            y.left = node
        else:
            y.right = node

        if node.parent is None:
            node.color = 0
            return

        if node.parent.parent is None:
            return

        self.fix_insert(node)
        
        return node

    def lower_bound(self, key):
        bound = None
        node = self.root
        while node is not None:
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
        while node is not None:
            if node.key >= key and (bound is None or bound.key >= node.key):
                bound = node
            if key < node.key:
                node = node.left
            else:
                node = node.right
        return bound
    
    def find(self, key):
        x = self.root
        
        while x is not None and x.key != key:
            if key < x.key:
                x = x.left
            else:
                x = x.right
                
        return x
    
    def size_helper(self, node):
        if node is None:
            return 0
        return self.size_helper(node.left) + self.size_helper(node.right) + 1
    
    def size(self):
        return self.size_helper(self.root)

    def get_root(self):
        return self.root

    def delete(self, key):
        self.delete_node_helper(self.root, key)

    def print_tree(self):
        self.__print_helper(self.root, "", True)
    
    def __prettyPrint_helper(self, node, indent):
        branch_len = 20
        if node is not None:
            if node.parent is not None:
                if node is node.parent.left:
                    left_indent = indent[:-1] + " " + " " * branch_len + "|"
                    right_indent = indent + " " * branch_len + "|"
                else:
                    right_indent = indent[:-1] + " " + " " * branch_len + "|"
                    left_indent = indent + " " * branch_len + "|"
            else:
                left_indent = " " * branch_len + "|"
                right_indent = " " * branch_len + "|"
                
            self.__prettyPrint_helper(node.right, right_indent)
            print("{}{}{}:{} ".format(indent, "-" * branch_len, node.key, node.value)) #node.parent.key if node.parent is not None else None
            self.__prettyPrint_helper(node.left, left_indent)
    
    def prettyPrint(self):
        self.__prettyPrint_helper(self.root, "")


if __name__ == "__main__":
    bst = RedBlackTree()

    bst.insert(55, 1)
    bst.insert(40, 2)
    bst.insert(65, 3)
    bst.insert(60, 4)
    bst.insert(75, 5)
    bst.insert(57, 6)

    bst.print_tree()

    print("\nAfter deleting an element")
    bst.delete(40)
    bst.print_tree()