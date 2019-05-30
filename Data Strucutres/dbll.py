"""This module provides an implementation of a doubly linked
list data structure.
As a reference the code fragments in Data Structures and
Algorithms in Python, M.T. Goodrich, R. Tamassia, and M.H. Goldwasser from
2013 were used.
"""


class dbll_node:
    """A node object that can be part of a doubly linked list.

    Attributes:
        element: Integer value the node holds.
        prev: Previous node in the doubly linked list.
        next: Next node in the doubly linked list.
    """

    # Since we only expect specific instance attributes we use __slots__
    # to save memory and accelerate attribute access.
    __slots__ = 'element', 'prev', 'next'

    def __init__(self, element, prev, next):
        """Initializes a dbll_node with the aforementioned attributes."""
        self.element = element
        self.prev = prev
        self.next = next

    def __next__(self):
        return self.next

    def __repr__(self):
        return '{}'.format(self.element)


class DoublyLinkedList:
    """A doubly linked list consisting of a trailer, header and dbll_node instances.

    In this implementation two sentinel nodes, _trailer and _header
    are used to make the implementation a little bit less complex
    when regarding edge cases.

    Attributes:
        _header: Sentinel node at the head of the doubly linked list.
        _trailer: Sentinel node at the beginning of the doubly linked list.
        _size: Amount of nodes in the doubly linked list aside from
            the two sentinal nodes.
    """

    def __init__(self):
        """Initializes a doubly linked list with the aforementioned attributes."""

        self._header = dbll_node(None, None, None)
        self._trailer = dbll_node(None, None, None)
        self._header.next = self._trailer
        self._trailer.next = self._header
        self._size = 0

    def __repr__(self):
        if self._trailer.next is self._header:
            return 'Empty dbll'
        return f'{[nodes for nodes in self if nodes is not self._header]}'

    def __iter__(self):
        current = self._trailer.next
        while current.element is not None:
            yield current
            current = current.next

    def __len__(self):
        return self._size

    def dbll_core_append(self, new_node):
        """Creates a new dbll_node and appends it to the given dbll.

        This function is used in the core_decomp function in 'order.py'.
        Since the nodes in the graph are given as integers we need to
        create a new dbll_node for each given graph node to construct
        the initial k_order.

        Args:
            self: Doubly linked list to be appended to.
            new_node: Integer used to create a new dbll_node.

        Returns:
            Return the newly created dbll_node object.
        """

        if self._size == 0:
            new = dbll_node(new_node, self._trailer, self._header)
            self._trailer.next = new
            self._header.prev = new
            self._size += 1
        else:
            new = dbll_node(new_node, self._header.prev, self._header)
            self._header.prev.next = new
            self._header.prev = new
            self._size += 1
        return new

    def dbll_append(self, new_node):
        """Appends a dbll_node to the given dbll and updates its size."""

        if self._size == 0:
            new_node.next = self._header
            new_node.prev = self._trailer
            self._trailer.next = new_node
            self._header.prev = new_node
        else:
            new_node.next = self._header
            new_node.prev = self._header.prev
            self._header.prev.next = new_node
            self._header.prev = new_node
        self._size += 1
        return new_node

    def dbll_prepend(self, new_node):
        """Prepends a dbll_node to the given dbll and updates its size."""

        new_node.next = self._trailer.next
        new_node.prev = self._trailer
        self._trailer.next.prev = new_node
        self._trailer.next = new_node
        self._size += 1
        return new_node

    def dbll_remove(self, node):
        """Removes a dbll_node from the given dbll and updates its size."""

        node.prev.next = node.next
        node.next.prev = node.prev
        node.prev = node.next = None
        # self._size -= 1
        return node

    def dbll_terminate(self, node_1, node_2, destination_dbll):
        """Appends all remaining dbbl nodes in self to destination_dbll."""

        if destination_dbll._size == 0:
            destination_dbll._trailer.next = node_1
            node_1.prev = destination_dbll._trailer
            destination_dbll._header.prev = node_2
            node_2.next = destination_dbll._header
        else:
            node_1.prev = destination_dbll._header.prev
            destination_dbll._header.prev.next = node_1
            node_2.next = destination_dbll._header
            destination_dbll._header.prev = node_2
        # destination_dbll._size += self._size
        # self._size = 0
        destination_dbll._size = 1

    def dbll_block_append(self, node_1, node_2, destination_dbll, size_change):
        """Removes the dbll nodes in [node_1, node_2] from self and appends them to destination_dbll"""

        if destination_dbll._size == 0:
            node_1.prev.next = node_2.next
            node_2.next.prev = node_1.prev
            destination_dbll._trailer.next = node_1
            node_1.prev = destination_dbll._trailer
            destination_dbll._header.prev = node_2
            node_2.next = destination_dbll._header
        else:
            node_1.prev.next = node_2.next
            node_2.next.prev = node_1.prev
            node_1.prev = destination_dbll._header.prev
            destination_dbll._header.prev.next = node_1
            node_2.next = destination_dbll._header
            destination_dbll._header.prev = node_2
        # self._size -= size_change
        # destination_dbll._size += size_change
        destination_dbll._size = 1

    def dbll_obs_6(self, append_to, node):
        node.next = append_to.next
        node.prev = append_to
        append_to.next.prev = node
        append_to.next = node
        # self._size += 1
