import sys

class Node:
    def  __init__(self, key, value):
        self.key = key
        self.value = value
        self.parent = None
        self.left = None
        self.right = None
        self.bf = 0

class AVLTree:

    def __init__(self):
        self.root = None
    
    def __searchTreeHelper(self, node, key):
        if node is None or key == node.key:
            return node

        if key < node.key:
            return self.__searchTreeHelper(node.left, key)
        else:
            return self.__searchTreeHelper(node.right, key)

    def __deleteNodeHelper(self, node, key):
        # search the key
        if node is None: 
            return node
        elif key < node.key:
            node.left = self.__deleteNodeHelper(node.left, key)
            if node.left is not None:
                node.left.parent = node
        elif key > node.key: 
            node.right = self.__deleteNodeHelper(node.right, key)
            if node.right is not None:
                node.right.parent = node
        else:
            # the key has been found, now delete it

            # case 1: node is a leaf node
            if node.left is None and node.right is None:
                node = None

            # case 2: node has only one child
            elif node.left is None:
                temp = node
                node = node.right

            elif node.right is None:
                temp = node
                node = node.left

            # case 3: has both children
            else:
                temp = self.minimum(node.right)
                node.key = temp.key
                node.value = temp.value
                node.right = self.__deleteNodeHelper(node.right, temp.key)
                if node.right is not None:
                    node.right.parent = node

            # Write the update balance logic here 
            # YOUR CODE HERE
        return node

    # update the balance factor the node
    def __updateBalance(self, node):
        if node.bf < -1 or node.bf > 1:
            self.__rebalance(node)
            return;

        if node.parent is not None:
            if node == node.parent.left:
                node.parent.bf -= 1

            if node == node.parent.right:
                node.parent.bf += 1

            if node.parent.bf != 0:
                self.__updateBalance(node.parent)

     # rebalance the tree
    def __rebalance(self, node):
        if node.bf > 0:
            if node.right.bf < 0:
                self.rightRotate(node.right)
                self.leftRotate(node)
            else:
                self.leftRotate(node)
        elif node.bf < 0:
            if node.left.bf > 0:
                self.leftRotate(node.left)
                self.rightRotate(node)
            else:
                self.rightRotate(node)
                
                
    # Printing the tree
    def __print_helper(self, node, indent):
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
                
            self.__print_helper(node.right, right_indent)
            print("{}{} {}:{} ".format(indent, "-" * branch_len, node.key, node.value)) #node.parent.key if node.parent is not None else None
            self.__print_helper(node.left, left_indent)
        
    def size_helper(self, node):
        if node is None:
            return 0
        return self.size_helper(node.left) + self.size_helper(node.right) + 1

    # rotate left at node x
    def leftRotate(self, x):
        y = x.right
        x.right = y.left
        if y.left is not None:
            y.left.parent = x

        y.parent = x.parent;
        if x.parent is None:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
            if y is not None:
                y.parent = x.parent
        else:
            x.parent.right = y
            if y is not None:
                y.parent = x.parent
        y.left = x
        x.parent = y

        # update the balance factor
        x.bf = x.bf - 1 - max(0, y.bf)
        y.bf = y.bf - 1 + min(0, x.bf)

    # rotate right at node x
    def rightRotate(self, x):
        y = x.left
        x.left = y.right;
        if y.right is not None:
            y.right.parent = x
        
        y.parent = x.parent;
        if x.parent is None:
            self.root = y
        elif x == x.parent.right:
            x.parent.right = y
            if y is not None:
                y.parent = x.parent
        else:
            x.parent.left = y
            if y is not None:
                y.parent = x.parent
        
        y.right = x
        x.parent = y

        # update the balance factor
        x.bf = x.bf + 1 - min(0, y.bf)
        y.bf = y.bf + 1 + max(0, x.bf)

    # search the tree for the key k
    # and return the corresponding node
    def find(self, k):
        return self.__searchTreeHelper(self.root, k)
    
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

    # find the node with the minimum key
    def minimum(self, node):
        if node is None:
            return None
        while node.left is not None:
            node = node.left
        return node

    # find the node with the maximum key
    def maximum(self, node):
        if node is None:
            return None
        while node.right is not None:
            node = node.right
        return node

    # find the successor of a given node
    def successor(self, x:
        if debug:
            print("\tgetting successor to {}: {}".format(x.key, x.value))
        # if the right subtree is not null,
        # the successor is the leftmost node in the
        # right subtree
        if x.right is not None:
            if debug:
                print("\tfinding minimum of right tree")
            return self.minimum(x.right)

        s = None
        while x.parent is not None:
            if x.parent.left == x:
                s = x.parent
                break
            x = x.parent
        return s

    # find the predecessor of a given node
    def predecessor(self, x):
        # if the left subtree is not null,
        # the predecessor is the rightmost node in the 
        # left subtree
        if x.left is not None:
            return self.maximum(x.left)

        s = None
        while x.parent is not None:
            if x.parent.right == x:
                s = x.parent
                break
            x = x.parent
        return s

    # insert the key to the tree in its appropriate position
    def insert(self, key, value):
        # PART 1: Ordinary BST insert
        node =  Node(key, value)
        y = None
        x = self.root

        while x is not None:
            y = x
            if node.key < x.key:
                x = x.left
            else:
                x = x.right

        # y is parent of x
        node.parent = y
        if y is None:
            self.root = node
        elif node.key < y.key:
            y.left = node
            if node is not None:
                node.parent = y
        else:
            y.right = node
            if node is not None:
                node.parent = y

        # PART 2: re-balance the node if necessary
        self.__updateBalance(node)
        return node


    # delete the node from the tree
    def delete(self, key):
        self.root = self.__deleteNodeHelper(self.root, key)

    # print the tree structure on the screen
    def prettyPrint(self):
        self.__print_helper(self.root, "")
        
        
    def size(self):
        return self.size_helper(self.root)
