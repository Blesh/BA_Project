"""This script provides an implementation of the algorithms specified in the
2016 paper 'A Fast Order-Based Approach for Core Maintenance' by
Yikai Zhang, Jeffrey Xu Yu, Ying Zhang, and Lu Qin.

Imports:
    networkx: Provides functionality for the creation, manipulation, and
            modification of graphs.

    heapq: This module provides an implementation of the heap queue algorithm,
        also known as the priority queue algorithm. Used for the
        implementation of the min heap in order_insert.

    queue: The module implements three types of queues, which differ only in
        the order in which the entries are retrieved. In this script the
        standard FIFO queue is used.

    dbll: Since as of Python 3.6.5, python doesn't have a built-in doubly
        linked list data structure, the dbll module provides an
        implementation of a doubly linked class in order to create the
        k-order defined in the paper.

    treap: This module provides an implementation of treaps based on the
        the 1989 paper "Randomized Search Trees", by Aragon, Cecilia R.,
        and Raimund G. Seidel.

    timeit: This module provides a simple way to time Python code.

    collections: The order in which the nodes are processed plays a huge role
                in the order-based algorithms. Thus from collections the OrderDict is used.
                The major advantage of an ordered dict in this case, is that
                the order is maintained by itself and membership can be checked
                in constant time.
"""
import networkx as nx
import heapq
import queue
import dbll
import treap
import random
import timeit
import collections


"""The entries in the list 'k_order' are doubly linked lists that hold
dbll nodes representing their graph node counterparts.
The nodes have core numbers respectively to the index of the
doubly linked list they are in.
E.g.: The doubly linked list at k_order[1] holds the dbll nodes that
correspond to the graph nodes with core number 1.
"""
k_order = []

"""Treap_order holds the treaps representing the doubly linked lists
in the k_order.
The treap at treap_order[x] holds the treap that represents the
doubly linked list at k_order[x].
"""
treap_order = []

"""The keys in the treap_map dictionary are the values of all dbll nodes
in the k_order and the values are their corresponding treap_node objects.
If there is a dbll node x in k_order[k], the treap_map dictionary maps
this dbll node to its respective treap node in treap_order[k].
The dictionary provides a fast reference to each treap node so that
we can calculate their ranks efficiently.
"""
treap_map = {}

"""The keys in the dbll_map dictionary are integers that represent
each node in the graph. Those keys are mapped to their respective
dbll node in the k_order.
This way we have a fast reference to each dbll node in the k_order,
which is needed in order_remove.
"""
dbll_map = {}


# Source:
# https://stackoverflow.com/questions/1581895/how-check-if-a-task-is-already-in-python-queue/1581937#1581937
class Membership_queue(queue.Queue):
    """Subclass extending FIFO queue so that one can check for membership.

    The hook methods _init and _put are overwritten, so that with each queue
    there is a corresponding set. This way duplicate items are not allowed.
    Could also be achieved with an external set, but this keeps the code more
    readable.
    """

    def _init(self, maxsize):
        """Inits Membership_queue with a set."""
        self.maxsize = maxsize
        self.queue = set()

    def _put(self, item):
        """Adds an item to the set and the queue."""
        self.queue.add(item)

    def _get(self):
        """Pops an item in FIFO fashion."""
        return self.queue.pop()


def core_decomp(graph):
    """Provides the k-core decomposition of an input graph and constructs the k_order.

    The implementation is based on the algorithm provided in the 2002
    paper "An O(m) Algorithm for Cores Decomposition of Networks.", by
    Vladimir Batagelj, and Matjaz Zaversnik. Some modifications have
    been made, based on Algorithm 1 in the 2016 paper to tailor to the
    Order-Based approach.

    The nodes of the graph are inserted into the list 'nodes' in
    ascending order of their degrees. After that, they are iterated over and
    for each current node, the core number is set equal to their current
    degree, and each neighbor of that node gets its degree decreased by 1.

    While iterating over each node, the remaining degree and candidate degree
    are initialized for that vertex.

    For each currently processed vertex a dbll node is created and appended to
    k_order[k], so that the index of the doubly linked list in the k_order, represents
    the core number of the nodes in it.
    Furthermore, there is an entry created in dbll_map, mapping the current vertex to
    the newly created dbll node.

    Args:
        graph: Graph object created with the networkx library.
    """

    # Initialization
    copy_graph = graph.copy()
    max_node = max(graph)
    bins = []
    maximum_degree = 0
    # Init necessary lists with enough space for all vertices in the graph.
    nodes = [-1]*(max_node + 2)
    position = [0]*(max_node + 1)
    core = [-1]*(max_node + 1)
    for node in graph:
        # Compute the degree of each vertex.
        deg = 0
        for neighbor in graph[node]:
            deg += 1
        core[node] = deg
        if deg > maximum_degree:
            maximum_degree = deg
    for i in range(maximum_degree + 1):
        bins.append(0)  # Init bins.
    for node in graph:
        bins[core[node]] += 1  # Count how many nodes are in each bin.
    start = 1
    for i in range(maximum_degree + 1):
        # Determine at what position each bin starts.
        current_bin = bins[i]
        bins[i] = start
        start += current_bin
    for node in graph:
        # Order the nodes in the graph based on their degree.
        position[node] = bins[core[node]]
        nodes[position[node]] = node
        bins[core[node]] += 1
    for i in range(maximum_degree, 1, -1):
        # Reset the starting position of each bin.
        bins[i] = bins[i-1]
    bins[0] = 1
    k = 1
    graph_nodes = (node for node in nodes if node != -1)
    c_node = next(graph_nodes)
    try:
        while True:
            k_order.append(dbll.DoublyLinkedList())
            while core[c_node] < k:
                # Init remaining and candidate degree.
                graph.nodes[c_node]['remaining_degree'] = len(copy_graph[c_node])
                graph.nodes[c_node]['candidate_degree'] = 0
                # Create a dbll node for the current graph node and append it to the k_order.
                new_node = k_order[k-1].dbll_core_append(c_node)
                dbll_map[c_node] = new_node
                for neighbor in graph[c_node]:
                    if core[neighbor] > core[c_node]:
                        # Move the neighbor node one bin to the left by swapping its
                        # position with the first vertex in the current bin and shift the starting
                        # position of the bin (l.184-193).
                        neighbor_degree = core[neighbor]
                        position_neighbor = position[neighbor]
                        current_bin = bins[neighbor_degree]
                        first_vertex = nodes[current_bin]
                        if neighbor != first_vertex:
                            position[neighbor] = current_bin
                            position[first_vertex] = position_neighbor
                            nodes[position_neighbor] = first_vertex
                            nodes[current_bin] = neighbor
                        bins[neighbor_degree] += 1
                        core[neighbor] -= 1
                copy_graph.remove_node(c_node)
                c_node = next(graph_nodes)
            k += 1
    except StopIteration:
        return core


def treap_setup(k_order):
    """Creates a treap_order corresponding to a given k_order.

    For each doubly linked list and their respective dbll nodes
    in k_order a treap is constructed and inserted into treap_order.

    While iterating over each dbll node a corresponding treap_node is created and
    an entry in treap_map is created, which maps the dbll node values to the treap node.
    The rank of the treap_node is set to the position of the currently
    processed dbll node in its doubly linked list.
    Finally, the treap node is inserted in treap_order[k], where k is the index of the
    doubly linked list in k_order.
    After the function is done, each index in k_order and treap_order
    represent their respective treap and doubly linked list.

    E.g.: Assume we have the doubly linked list k_order[1] = [11, 23, 5, 71].
    Then treap_setup creates a treap in treap_order[1] with rank(11) = 1,
    rank(23) = 2, rank(5) = 3, and rank(71) = 4.

    Args:
        k_order: List of doubly linked lists representing the order
                the graph nodes were processed in core_decomp.
    """

    for i in range(len(k_order)):
        rank = 0
        treap_order.append(treap.Treap())
        for nodes in k_order[i]:
            # Create treap_node corresponding to current dbll node.
            place_holder = treap.treap_node(nodes.element)
            # Map the value of the dbll node to the treap node object.
            treap_map[nodes.element] = place_holder
            # Set the treap node rank to the current position of the dbll in k_order[i].
            place_holder.rank = rank
            treap_order[i].treap_insert(treap_order[i].root, place_holder)
            rank += 1


def compute_mcd_init(graph):
    """Computes the maximum-core degree of all nodes in the graph.

    The maximum-core degree(mcd) of a node is defined as the
    number of its neighbors with equal or higher core
    number. The mcd is saved as a dict entry in the dictionary of
    each node.

    Args:
        graph: Graph object created with the networkx library.
    """

    for node in graph:
        mcd_count = 0
        for neighbor in graph[node]:
            if core[neighbor] >= core[node]:
                mcd_count += 1
        graph.nodes[node]['mcd'] = mcd_count


def compute_mcd(node):
    """Computes the maximum-core degree of a single node in the graph.

    Since order_remove uses the maximum-core degree of the nodes, this
    function is needed to update the maximum-core degree of certain nodes after
    order_insert and order_remove.

    Args:
        node: Integer representing a node in the graph.
    """

    mcd_count = 0
    for neighbor in graph[node]:
        if core[neighbor] >= core[node]:
            mcd_count += 1
    graph.nodes[node]['mcd'] = mcd_count


def order_insert(graph, u, v):
    """Inserts an edge between two given graph nodes and updates the necessary core numbers.

    Adds an edge between graph node u and graph node v. The graph nodes, which
    need their core numbers increased as a result of the insertion are updated.
    The graph nodes x are classified based on their remaining_degree and
    candidate_degree:

    - Branch 1(l.337): x.remaining_degree + x.candidate_degree > k.
    - Branch 2(l.346): x.remaining_degree + x.candidate_degree <= k and x.candidate_degree = 0.
    - Branch 3(l.573): x.remaining_degree + x.candidate_degree <= k and x.candidate_degree not 0.

    Graph nodes that have their core number increased have their corresponding
    dbll nodes prepended to k_order[k+1]. Accordingly the corresponding treap
    nodes are prepended to treap_order[k+1].

    Args:
        graph: Graph object created with the networkx library.
        u, v: Integers representing the graph nodes between which the edge is to be inserted.
    """

    # Initialization.
    B = []  # Stack B to maintain the next vertex to jump to in Branch 2.
    b_set = set()  # Used to lazily delete elements in B.
    obs_6 = {}
    candidates = collections.OrderedDict()
    new_ok = dbll.DoublyLinkedList()
    graph.add_edge(u, v)
    if core[u] < core[v]:
        graph.nodes[u]['remaining_degree'] += 1
        graph.nodes[u]['mcd'] += 1
        k = core[u]
        root = u
    elif core[u] > core[v]:
        graph.nodes[v]['remaining_degree'] += 1
        graph.nodes[v]['mcd'] += 1
        k = core[v]
        root = v
    else:
        k = core[u]
        graph.nodes[v]['mcd'] += 1
        graph.nodes[u]['mcd'] += 1
        if treap_order[k].rank(u, treap_map) < treap_order[k].rank(v, treap_map):
            root = u
            graph.nodes[u]['remaining_degree'] += 1
        else:
            root = v
            graph.nodes[v]['remaining_degree'] += 1

    # Core phase.
    if graph.nodes[root]['remaining_degree'] <= k:
        return
    # heapq.heappush(B, (treap_order[k].rank(root, treap_map), root))
    c_node = k_order[k]._trailer.next
    if dbll_map[root].prev is not k_order[k]._trailer:
        k_order[k].dbll_block_append(c_node, dbll_map[root].prev, new_ok, treap_order[k].rank(root, treap_map))
        c_node = k_order[k]._trailer.next
    while c_node is not k_order[k]._header:
        while len(B) != 0 and B[0][1] in b_set:
            # Lazy removal of nodes in B.
            heapq.heappop(B)
        b_set.add(c_node.element)
        # Branch 1.
        if (graph.nodes[c_node.element]['remaining_degree'] + graph.nodes[c_node.element]['candidate_degree'] > k):
            # Add the current node to candidates.
            candidates[k_order[k].dbll_remove(c_node)] = None
            for neighbor in graph[c_node.element]:
                if (core[neighbor] == k and (treap_order[k].rank(c_node.element, treap_map) <
                                             treap_order[k].rank(neighbor, treap_map))):
                    graph.nodes[neighbor]['candidate_degree'] += 1
                    # Since the neighbors canidate_degee is greater than 0 it gets pushed on
                    # the heap B.
                    heapq.heappush(B, (treap_order[k].rank(neighbor, treap_map), neighbor))
            c_node = k_order[k]._trailer.next
        # Branch 2
        elif graph.nodes[c_node.element]['candidate_degree'] == 0:
            try:
                j, jump_node = heapq.heappop(B)
                # All nodes from c_node to jump_node.prev are moved to new_ok
                k_order[k].dbll_block_append(c_node, dbll_map[jump_node].prev, new_ok, j)
                c_node = k_order[k]._trailer.next
            except IndexError:
                # B is empty so all remaining nodes in k_order[k] are appended to new_ok
                # and the algorithm leaves the loop.
                k_order[k].dbll_terminate(c_node, k_order[k]._header.prev, new_ok)
                break
        # Branch 3.
        else:
            new_ok.dbll_append(k_order[k].dbll_remove(c_node))
            graph.nodes[c_node.element]['remaining_degree'] += graph.nodes[c_node.element]['candidate_degree']
            graph.nodes[c_node.element]['candidate_degree'] = 0
            remove_candidates(graph, new_ok, b_set, c_node, k, obs_6, candidates)
            c_node = k_order[k]._trailer.next

    # Ending phase.
    for node in obs_6:
        # Adjust the position of the nodes in obs_6.
        treap_order[k].observation_6(treap_map[node.element], treap_map[obs_6[node].element])
    for node in reversed(candidates):  # reversed() to maintain the order the candidates had in k_order[k].
        graph.nodes[node.element]['candidate_degree'] = 0
        core[node.element] += 1
        compute_mcd(node.element)
        for neighbor in graph[node.element]:
            if (core[neighbor] == core[node.element] or core[neighbor] == core[node.element]-1):
                compute_mcd(neighbor)
        k_order[k+1].dbll_prepend(node)
        treap_order[k].treap_prepend(treap_order[k+1], treap_map[node.element])
    k_order[k] = new_ok


def remove_candidates(graph, new_ok, b_set, node, k, obs_6, candidates):
    """Adjusts all remaining and candidate_degrees that are affected by 'node' being in Branch 3.

    When a node in order_insert is processed and ends up in the third branch,
    that means that the sum of their candidate_degree and remaining_degree is
    not bigger than k, but their candidate degree is not 0. So there is a
    neighbor node before that node in the k_order that is a possible candidate.

    Since the input node cannot be in the (k+1)-core, the neighbor loses a potential
    neighbor in the (k+1)-core and therefore has its remaining degree decreased by one.
    That can result in the neighbor losing his place in candidates and may cascade onto
    to the neighbor of the neighbor and so forth.

    Args:
        graph: Graph object created with the networkx library.
        candidates: OrderedDict of graph nodes that are curently candidates.
        new_ok: Doubly linked list that replaces k_order[k] after oder_insert.
        b_set: Set of graph nodes to lazily delete them from heap B.
        node: Currently processed dbll node in oder_insert.
        k: Minimum core number of the two nodes between which the edge was inserted.
        obs_6: Dictionary that maps the nodes that are removed from candidates
            to the currently last node in new_ok.
    """

    q = Membership_queue()
    for neighbor in graph[node.element]:
        if dbll_map[neighbor] in candidates:
            graph.nodes[neighbor]['remaining_degree'] -= 1
            if (graph.nodes[neighbor]['remaining_degree'] +
                    graph.nodes[neighbor]['candidate_degree'] <= k):
                q.put(dbll_map[neighbor])
    while not q.empty():
        c_node = q.get()
        graph.nodes[c_node.element]['remaining_degree'] += graph.nodes[c_node.element]['candidate_degree']
        graph.nodes[c_node.element]['candidate_degree'] = 0
        del candidates[c_node]
        obs_6[c_node] = new_ok._header.prev  # Create obs_6 entry to later adjust treap_node positions.
        new_ok.dbll_append(c_node)
        k_neighbor = (x for x in graph[c_node.element] if core[x] == k)
        for neighbor in k_neighbor:
            neighbor_node = dbll_map[neighbor]
            if (treap_order[k].rank(node.element, treap_map) < treap_order[k].rank(neighbor, treap_map)):
                graph.nodes[neighbor]['candidate_degree'] -= 1
                if graph.nodes[neighbor]['candidate_degree'] == 0:
                    b_set.add(neighbor)  # Lazy removal of 'neighbor' from B.
            elif (neighbor_node in candidates and (treap_order[k].rank(c_node.element, treap_map) <
                                                   treap_order[k].rank(neighbor, treap_map))):
                graph.nodes[neighbor]['candidate_degree'] -= 1
                if (graph.nodes[neighbor]['candidate_degree'] +
                        graph.nodes[neighbor]['remaining_degree'] <= k and
                        (neighbor_node not in q.queue)):
                    q.put(neighbor_node)
            elif neighbor_node in candidates:
                graph.nodes[neighbor]['remaining_degree'] -= 1
                if (graph.nodes[neighbor]['candidate_degree'] +
                        graph.nodes[neighbor]['remaining_degree'] <= k and
                        (neighbor_node not in q.queue)):
                    q.put(neighbor_node)


def order_remove(graph, u, v):
    """Removes the edge between the two given graph nodes and updates the necessary core numbers.

    The edge between the two graph nodes u and v is removed. The core numbers
    of affected nodes are updated and the remaining_degree and candidate_degree
    of each one of those nodes is adjusted to the new k_order.
    Dbll nodes that have their core number decreased get appended to
    k_order[k-1]. Accordingly, the corresponding treap nodes are
    appended to treap_order[k-1].

    Args:
        graph: Graph object created with the networkx library.
        u, v: Integers representing the graph nodes between which the edge is to be removed.
    """

    # Initialization.
    candidates = collections.OrderedDict()
    propagation_queue = Membership_queue()
    graph.remove_edge(u, v)
    if core[u] < core[v]:
        graph.node[u]['remaining_degree'] -= 1
        graph.nodes[u]['mcd'] -= 1
        k = core[u]
        propagation_queue.put(u)
    elif core[u] > core[v]:
        graph.node[v]['remaining_degree'] -= 1
        graph.nodes[v]['mcd'] -= 1
        k = core[v]
        propagation_queue.put(v)
    else:
        graph.nodes[u]['mcd'] -= 1
        graph.nodes[v]['mcd'] -= 1
        k = core[u]
        propagation_queue.put(v)
        propagation_queue.put(u)
        if treap_order[k].rank(u, treap_map) < treap_order[k].rank(v, treap_map):
            graph.nodes[u]['remaining_degree'] -= 1
        else:
            graph.nodes[v]['remaining_degree'] -= 1

    # Core phase.
    # Updating core number of affected nodes based on mcd values.
    # Same approach as in traversal_remove.
    while not propagation_queue.empty():
        c_node = propagation_queue.get()
        if core[c_node] == k and graph.nodes[c_node]['mcd'] < k:
            candidates[c_node] = None
            core[c_node] -= 1
            for neighbor in graph[c_node]:
                if core[neighbor] == k:
                    graph.nodes[neighbor]['mcd'] -= 1
                    if graph.nodes[neighbor]['mcd'] < k:
                        propagation_queue.put(neighbor)

    # Ending phase.
    # Update k_order and treap_order.
    # Create a copy of candidates since we cannot iterate over a dict and alter it at the same time.
    copy_candidates = {**candidates}
    for node in candidates:
        # compute_mcd(node)
        graph.nodes[node]['remaining_degree'] = 0
        for neighbor in graph[node]:
            if (core[neighbor] == k and (treap_order[k].rank(neighbor, treap_map) <
                                         treap_order[k].rank(node, treap_map))):
                graph.nodes[neighbor]['remaining_degree'] -= 1
            if core[neighbor] >= k or (neighbor in copy_candidates):
                graph.nodes[node]['remaining_degree'] += 1
            if (core[neighbor] == core[node] or core[neighbor] == core[node]+1):
                compute_mcd(neighbor)
        del copy_candidates[node]
        k_order[k-1].dbll_append(k_order[k].dbll_remove(dbll_map[node]))
        treap_order[k].treap_append(treap_order[k-1], treap_map[node])


# Insert the graph definition here:

def file_to_graph(file):
    """Creates a networkx graph through the given file.

    The file has to consist of two columns, where two integers in one line
    represent an edge in the graph.

    E.g.: The file has to look like this:
                        1   2
                        5   6
                        6   10
                        ...

    Args:
        file: Text file containing integers in the format mentioned above.
    """

    # graph = DefaultGraph()  # Uncomment for traversal.
    graph = nx.Graph()  # Uncomment for order and color.

    with open(file) as f:
        for line in f:
            node_1, node_2 = tuple(line.split())
            graph.add_edge(int(node_1), int(node_2))
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph


# Initialization
core = core_decomp(graph)
treap_setup(k_order)
compute_mcd_init(graph)


# After the graph definition and initialization the code to test the script
# can be inserted here:
