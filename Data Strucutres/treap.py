"""This module provides an implementation of treaps based on the
the 1989 paper "Randomized Search Trees", by Aragon, Cecilia R.,
and Raimund G. Seidel.

Imports:
    queue: The module implements three types of queue, which differ only in
           the order in which the entries are retrieved.

    random: This module implements pseudo-random number generators for
        various distributions. Used to randomly generate a priority for the
        treap nodes.

"""
import random


class treap_node:
    """An object that represents a node in a treap.

    In this implementation, each treap_node corresponds to a
    node in the k_order and as a consequence a node in the graph.

    Attributes:
        rank: Position of the node in an in-order walk of the treap it is in.
        node: Integer identifying the node, which makes it easier to
            map it to a node in the graph.
        left: Left child treap_node.
        right: Right child treap_node.
        parent: Parent treap_node.
        size_left: Number of nodes in the left subtree rooted at the
            node +1 for the node itself.
        size_right: Number of nodes in the right subtree rooted at the
            node +1 for the node itself.
        priority: Random integer.
    """

    # Since we only expect specific instance attributes we use __slots__
    # to save memory and accelerate attribute access.
    __slots__ = ('rank', 'node', 'left', 'right', 'parent',
                 'size_left', 'size_right', 'priority')

    def __init__(self, node):
        """Initializes a treap_node with the aforementioned attributes."""

        self.rank = None
        self.node = node
        self.left = None
        self.right = None
        self.parent = None
        self.size_left = 1
        self.size_right = 1
        self.priority = int(random.random()*100_000_000)

    def __repr__(self):
        return (f'node: {self.node} | '
                f'prio: {self.priority} | '
                f'rank: {self.rank} | '
                f'parent: {getattr(self.parent, "node", None)} | '
                f'size_left: {self.size_left} | '
                f'size_right: {self.size_right} | '
                f'left: {getattr(self.left, "node", None)} | '
                f'right: {getattr(self.right, "node", None)}')

    def rotate_right(self):
        """Rotates the given treap_node to the right."""

        placeholder = self.left
        if placeholder.right is not None:
            self.size_left = placeholder.size_right
            placeholder.right.parent = self
        else:
            self.size_left = 1
        placeholder.size_right += self.size_right
        if self.parent is not None:
            if self.parent.right == self:
                self.parent.right = self.left
            else:
                self.parent.left = self.left
        self.left = placeholder.right
        placeholder.parent = self.parent
        self.parent = placeholder
        placeholder.right = self
        return placeholder

    def rotate_left(self):
        """Rotates the given treap_node to the left."""

        placeholder = self.right
        if placeholder.left is not None:
            self.size_right = placeholder.size_left
            placeholder.left.parent = self
        else:
            self.size_right = 1
        placeholder.size_left += self.size_left
        if self.parent is not None:
            if self.parent.right == self:
                self.parent.right = self.right
            else:
                self.parent.left = self.right
        self.right = placeholder.left
        placeholder.parent = self.parent
        self.parent = placeholder
        placeholder.left = self
        return placeholder


class Treap:
    """A treap is a randomly built binary search tree.

    Treaps use randomly chosen priorities to balance themselves.
    The treap_nodes are arranged in heap order regarding their priorities.
    Their keys obey the binary-search-tree property.

    That means for each treap_node x we have the following conditions regarding
    the treap_nodes y in its left subtree, the treap_nodes z in its right subtree,
    and its parent treap_node p:

    - x.key >= y.key
    - x.key <= z.key
    - x.priority <= p.priority

    In this implementation, the keys of the treap_nodes are their rank
    in an in-order walk of the tree.

    Attributes:
        root: Root treap_node of the treap.
    """

    def __init__(self, root=None):
        """Initializes the Treap object with a root."""
        self.root = root

    def treap_insert(self, root, node):
        """Inserts a treap_node into the treap and rebalances it.

        This function is an implementation of the TREAP-INSERT
        algorithm in the 1989 paper. First, the treap_node is inserted
        into the treap with respect to its rank and then rotations
        are performed to validate the tree in regard to the priorities
        of the treap_nodes.
        This implementation is iterative rather than recursive, as in the
        original paper. This way we don't have to change the maximum
        recursion depth and since, as of Python 3.6.5, there is no support
        for tail recursion elimination this approach is more efficient.
        (http://neopythonic.blogspot.com/2009/04/tail-recursion-elimination.html)

        Args:
            self: Treap the node is inserted in.
            root: Root node of the given treap.
            node: Treap_node to be inserted.

        Returns:
            Inserted Treap node.
        """

        if self.root is None:
            self.root = node
            return node
        c_node = self.root
        while c_node is not None:
            # Traverse the tree based on the rank of the node until
            # the node is a leaf.
            if node.rank > c_node.rank:
                if c_node.right is None:
                    c_node.size_right += 1
                    c_node.right = node
                    node.parent = c_node
                    break
                c_node.size_right += 1
                c_node = c_node.right
            else:
                if c_node.left is None:
                    c_node.size_left += 1
                    c_node.left = node
                    node.parent = c_node
                    break
                c_node.size_left += 1
                c_node = c_node.left
        c_parent = node.parent
        while (c_parent is not None and
                c_parent.priority < node.priority):
            # Perform rotations until the priorities are in the right place.
            if c_parent.right == node:
                higher_node = c_parent.rotate_left()
                c_parent = higher_node.parent
            else:
                higher_node = c_parent.rotate_right()
                c_parent = higher_node.parent
        if c_parent is None:
            self.root = higher_node

    def treap_remove(self, node):
        """Removes a given treap_node from the treap.

        This function implements the Treap-Delete function given
        in the 1989 paper. The implementation below is
        iterative rather than recursive, like it is in the paper,
        since as of Python 3.6.5 there is no support for tail recursion
        elimination.
        (http://neopythonic.blogspot.com/2009/04/tail-recursion-elimination.html)

        Args:
            self: Treap the treap_node is removed from.
            node: Treap_node that is removed.

        Returns:
            Treap_node that was removed.
        """

        # To remove a treap_node from the treap, the node is rotated downwards
        # based on the priorities of its children until it is a leaf, where we
        # can simply remove it.
        while True:
            if node is self.root:
                # For the case that the node is the root of the Treap there are 4 cases:
                # 1. No Children 2. Only has left child 3. Only has right child
                # 4. Has left and right child.
                if node.left is None and node.right is None:
                    self.root = None
                    return node
                if node.right is None:
                    node.rotate_right()
                    self.root = node.parent
                elif node.left is None:
                    node.rotate_left()
                    self.root = node.parent
                else:
                    if node.left.priority > node.right.priority:
                        node.rotate_right()
                        self.root = node.parent
                        continue
                    else:
                        node.rotate_left()
                        self.root = node.parent
                        continue
            if node.left is None and node.right is None:
                # When the node is a leaf, we remove it and adjust all necessary sizes.
                c_parent = node.parent
                c_node = node
                while c_parent is not None:
                    if c_parent.right == c_node:
                        c_parent.size_right -= 1
                        c_node = c_parent
                        c_parent = c_parent.parent
                    else:
                        c_parent.size_left -= 1
                        c_node = c_parent
                        c_parent = c_parent.parent
                if node.parent.left == node:
                    node.parent.left = None
                    node.parent = None
                    return node
                else:
                    node.parent.right = None
                    node.parent = None
                    return node
            elif node.left is None:
                node.rotate_left()
            elif node.right is None:
                node.rotate_right()
            else:
                if node.left.priority > node.right.priority:
                    node.rotate_right()
                else:
                    node.rotate_left()

    def observation_6(self, node, prev_node):
        """Removes and re-inserts a treap_node as the successor of prev_node in self.

        As mentioned in Observation 6.1 in the 2016 paper
        it is possible that after the insertion of an edge
        the positions of two nodes in the k_order may be changed.
        More precisely, a node node_1 that came before a node node_2
        in the k_order, now comes after node_2, with node_1 != node_2.

        E.g. Before the insertion k_order[2] looks like this
        [1, 2, 3], and after order_insert it may look like this [2, 1, 3].

        To make the treap consistent with the k_order, every node that
        was discarded from the candidate_set at some point, is inserted
        as the successor of the last node in new_ok, at that point in time.
        So the position of the node is consistent with the new k_order.

        Sizes and treap_node positions are adjusted accordingly.

        Args:
            self: Treap the node gets removed and reinserted in.
            node: Treap_node that has its position adjusted by deletion and
                re-insertion.
            prev_node: Predecessor treap_node in the new k-order.
        """

        self.treap_remove(node)
        c_node = prev_node.right
        if c_node is None:
            prev_node.right = node
            node.parent = prev_node
        else:
            while c_node.left is not None:
                # Iterate to the position of the successor of prev_node.
                c_node = c_node.left
            c_node.left = node
            node.parent = c_node
        c_parent = node.parent
        c_node = node
        while (c_parent is not None and
               c_parent.priority < node.priority):
            # Perform rotations until the treap is valid in regards to node priorities.
            if c_parent.right == node:
                c_node = c_parent.rotate_left()
                c_parent = c_node.parent
            else:
                c_node = c_parent.rotate_right()
                c_parent = c_node.parent
        if c_parent is None:
            self.root = c_node
            return
        while c_parent is not None:
            # Adjust sizes for the rest of the nodes in the treap.
            if c_parent.right == c_node:
                c_parent.size_right += 1
                c_node = c_parent
                c_parent = c_parent.parent
            else:
                c_parent.size_left += 1
                c_node = c_parent
                c_parent = c_parent.parent

    def treap_append(self, destination_treap, node):
            """Removes a treap_node from a Treap and appends it to another one.

            Args:
                self: Treap the treap_node is removed from.
                destination_treap: Treap the treap_node is appended to.
                node: Treap_node to be removed from the given treap and
                    appended to the destination_treap.
            """

            self.treap_remove(node)
            c_node = destination_treap.root
            # Appending the node to the treap.
            if c_node is None:
                destination_treap.root = node
                return
            while c_node.right is not None:
                c_node.size_right += 1
                c_node = c_node.right
            c_node.right = node
            c_node.size_right += 1
            node.parent = c_node
            c_parent = node.parent
            # Making the treap valid in regards to node priorities.
            while (c_parent is not None and
                   c_parent.priority < node.priority):
                higher_node = c_parent.rotate_left()
                c_parent = higher_node.parent
            if c_parent is None:
                destination_treap.root = higher_node

    def treap_prepend(self, destination_treap, node):
        """Removes a treap_node from a treap and prepends it to another one.

        Args:
            self: Treap the treap_node is removed from.
            destination_treap: Treap the treap_node is prepended to.
            node: Treap_node to be removed from the given treap and
                prepended to the destination_treap.
        """

        self.treap_remove(node)
        c_node = destination_treap.root
        # Prepending the node to the treap.
        if c_node is None:
            destination_treap.root = node
            return
        while c_node.left is not None:
            c_node.size_left += 1
            c_node = c_node.left
        c_node.left = node
        c_node.size_left += 1
        node.parent = c_node
        c_parent = node.parent
        # Making the treap valid in regards to node priorities.
        while (c_parent is not None and
               c_parent.priority < node.priority):
            higher_node = c_parent.rotate_right()
            c_parent = higher_node.parent
        if c_parent is None:
            destination_treap.root = higher_node

    def rank(self, node_number, treap_map):
        """Returns the rank of a given graph node.

        Based on the pseudocode provided in Introduction
        to Algorithms, Chapter 14.
        The rank represents the node's position in k_order[k],
        where k is the core number of the node. It is important to
        note that the rank of the treap node represents the
        position of the dbll node in its current doubly linked list.

        Args:
            self: Treap in which the rank is determined.
            node_number: Integer identifying the node in the graph.
            mapp: Dictionary that maps the node_number to the
                respective treap_node.

        Returns:
            Rank of the given node.
        """

        node = treap_map[node_number]
        rank = node.size_left
        c_node = node
        while c_node.parent is not None:
            if c_node == c_node.parent.right:
                rank += c_node.parent.size_left
            c_node = c_node.parent
        return rank
