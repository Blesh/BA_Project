"""This script provides an implementation of the algorithms specified in
'Incremental k-core decomposition: algorithms and evaluation.',
by Sariyüce  et al., published in
the VLDB Journal — The International Journal on Very Large Data Bases in 2016.

Imports:
    networkx: Provides functionality for the creation, manipulation, and
            modification of graphs.

    timeit: This module provides a simple way to time Python code.

    random: This module implements pseudo-random number generators for
            various distributions. Used to randomly choose nodes for our
            time performance evaluation.

    collections: In the paper the attributes of graph nodes are lazily initialized.
                We use the defaultdict class in collections to mimic this
                approach.
"""
import networkx as nx
import timeit
import random
from collections import defaultdict


class DefaultGraph(nx.Graph):
    """This class initializes inserted graph nodes with a defaultdict for their node attributes.

    We change the 'add_edges_from' and 'add_edge' functions of the Graph class
    in the networkx library. So whenever either of those functions is used
    to add edges, the nodes added to the Graph instance will have a defaultdict
    initialized for their attributes instead of a normal dictionary.

    The traversal algorithm uses lazy init for its node attributes, so this approach
    avoids multiple setdefault() calls or other tedious workarounds in the algorithm.
    """

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        u, v = u_of_edge, v_of_edge
        # add nodes.
        if u not in self._node:
            self._adj[u] = self.adjlist_inner_dict_factory()
            self._node[u] = defaultdict(int)  # defaultdict initialized.
        if v not in self._node:
            self._adj[v] = self.adjlist_inner_dict_factory()
            self._node[v] = defaultdict(int)  # defaultdict initialized.
        # add the edge.
        datadict = self._adj[u].get(v, self.edge_attr_dict_factory())
        datadict.update(attr)
        self._adj[u][v] = datadict
        self._adj[v][u] = datadict

    def add_edges_from(self, ebunch_to_add, **attr):
        for e in ebunch_to_add:
            ne = len(e)
            if ne == 3:
                u, v, dd = e
            elif ne == 2:
                u, v = e
                dd = {}
            else:
                raise nx.NetworkXError(
                    "Edge tuple %s must be a 2-tuple or 3-tuple." % (e,))
            if u not in self._node:
                self._adj[u] = self.adjlist_inner_dict_factory()
                self._node[u] = defaultdict(int)  # defaultdict initialized.
            if v not in self._node:
                self._adj[v] = self.adjlist_inner_dict_factory()
                self._node[v] = defaultdict(int)  # defaultdict initialized.
            datadict = self._adj[u].get(v, self.edge_attr_dict_factory())
            datadict.update(attr)
            datadict.update(dd)
            self._adj[u][v] = datadict
            self._adj[v][u] = datadict


def init_rcd(graph, n):
    """Initializes the residential core degrees of all nodes in the graph.

    The core degrees of all nodes in the graph are calculated, from the maximum-core degree
    up to the n-core degree. The integers representing those values are stored
    in a dictionary with the key being a node in the graph, and the value
    being a list of its residential core degrees.

    E.g.: graph.nodes[12]['rcd'] = [5, 4, 2]. So the node 12 has mcd 5, pcd 4,
        and 3-core degree 2.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer representing up to which residential core degree for each node
        is to be calculated.
    """

    for node in graph:
        graph.nodes[node]['rcd'] = [0]*n  # Inits a list with suitable size for each node.
    for i in range(n):
        for node in graph:
            if i == 0:  # Computation of the maximum-core degree.
                for neighbor in graph[node]:
                    if core[neighbor] >= core[node]:
                        graph.nodes[node]['rcd'][0] += 1
            else:
                for neighbor in graph[node]:
                    # Computation of n-core degree.
                    if ((core[neighbor] == core[node] and
                            graph.nodes[neighbor]['rcd'][i-1] > core[node]) or
                            core[node] < core[neighbor]):
                        graph.nodes[node]['rcd'][i] += 1


def compute_rcd(graph, node, n):
    """Computes the n-core degree for a given graph node.

    Args:
        graph: Graph object created with the networkx library.
        node: Integer representing a graph node.
        n: Integer representing which n-core degree should be calculated.
    """

    graph.nodes[node]['rcd'][n] = 0
    if n == 0:
        # Computation of the maximum-core degree.
        for neighbor in graph[node]:
            if core[neighbor] >= core[node]:
                graph.nodes[node]['rcd'][n] += 1
    else:
        # Computation of the n-core degree.
        for neighbor in graph[node]:
            if ((core[node] == core[neighbor] and
                    graph.nodes[neighbor]['rcd'][n-1] > core[node]) or
                    core[node] < core[neighbor]):
                graph.nodes[node]['rcd'][n] += 1


def traversal_insert(graph, n, u, v):
    """Inserts an edge between two given graph nodes and updates the necessary core numbers.

    An edge between the nodes u and v is inserted. For every node in graph, the core number
    is increased if necessary due to the insertion of the edge.
    To determine which nodes need to be updated the information gained
    by the residential core degrees up to the n-core degree is used.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer defining up to which n-hop neighborhood should be considered
           for updating the core numbers.
        u, v: Integers representing the graph nodes between which the edge is to be inserted.
    """

    # Initialization.
    root = u
    if core[v] < core[u]:
        root = v
    graph.add_edge(u, v)
    # Adjust residential core degrees after the edge of the insertion.
    multihop_prepare_rcd_insertion(graph, n, u, v)
    s = []  # Stack s.
    k = core[root]
    changed = set()  # Stores all the nodes that have their core number changed.
    reset_attr = set()  # Stores all the nodes that need their node attributes reset at the end.
    candidates = set()  # Stores all nodes that might need their core number changed.
    s.append(root)
    graph.nodes[root]['cd'] = graph.nodes[root]['rcd'][n-1]
    graph.nodes[root]['visited'] = 1
    candidates.add(root)
    pruned_vertices = set()

    # Core phase.
    while s:
        node = s.pop()
        reset_attr.add(node)
        if node == root:
            root_micd = graph.nodes[node]['micd']
            graph.nodes[node]['micd'] = 0
        if graph.nodes[node]['cd'] > k:  # and (graph.nodes[node]['micd'] < k or k == 0):
            # Node has enough suitable neighbors to increase its core number.
            # candidates.add(node)
            if graph.nodes[node]['micd'] == k:
                for neighbor in graph[node]:
                    if (core[neighbor] == k and graph.nodes[neighbor]['rcd'][0] > k and
                            not graph.nodes[neighbor]['evicted'] and not graph.nodes[neighbor]['visited']):
                        reset_attr.add(neighbor)
                        graph.nodes[neighbor]['pruned'].append(node)
                        pruned_vertices.add(neighbor)
            else:
                for neighbor in graph[node]:
                    if (core[neighbor] == k and graph.nodes[neighbor]['rcd'][n-2] > k and
                            not graph.nodes[neighbor]['visited']):
                        s.append(neighbor)
                        candidates.add(neighbor)
                        graph.nodes[neighbor]['visited'] = 1
                        # Due to propagate_eviction the 'cd' of the neighbor might have already
                        # been decreased so addition must be used to take that into account.
                        graph.nodes[neighbor]['cd'] += graph.nodes[neighbor]['rcd'][n-1]
        else:  # Node does not have enough suitable neighbors and therefore is evicted.
            if not graph.nodes[node]['evicted']:
                propagate_eviction(graph, k, node, candidates, reset_attr)
    for node in pruned_vertices:
        if not graph.nodes[node]['visited']:
            for z_nodes in graph.nodes[node]['pruned']:
                graph.nodes[z_nodes]['cd'] -= 1
                if graph.nodes[z_nodes]['cd'] == k:
                    propagate_eviction(graph, k, z_nodes, candidates, reset_attr)

    # Ending phase.
    graph.nodes[root]['micd'] = root_micd
    for node in reset_attr:
        if node in candidates and not graph.nodes[node]['evicted']:
            changed.add(node)
            for neighbor in graph[node]:
                if core[neighbor] == k+1:
                    graph.nodes[node]['micd'] -= 1
                if core[neighbor] == k:
                    graph.nodes[neighbor]['micd'] += 1
            core[node] += 1
        # Every node attribute besides 'rcd' has to be reset, since they are not
        # initialized at each function call. So without resetting further function calls
        # might result in errors.
        rcd = graph.nodes[node]['rcd']
        micd = graph.nodes[node]['micd']
        graph.nodes[node].clear()
        graph.nodes[node]['pruned'] = []
        graph.nodes[node]['rcd'] = rcd
        graph.nodes[node]['micd'] = micd
    # Adjust residential core degrees of the updated vertices and its neighbors.
    multihop_recompute_rcd_insertion(graph, n, changed)


def propagate_eviction(graph, k, node, candidates, reset_attr):
    """Evicts the given node and additionally any other node that needs to be evicted as a result.

    The given graph node is evicted since it does not have enough suitable neighbors
    to be in the (k+1)-core. Every neighbor has its cd reduced due to them losing a potential
    neighbor.That means the neighbors of the input node might not be in
    the (k+1)-core themselves. Resulting in the neighbor being evicted as well.
    This eviction mechanism propagates recursively until no other node needs to be evicted.

    Args:
        graph: Graph object created with the networkx library.
        k: Minimum core number of the two nodes between which the edge was inserted.
        node: Integer representing the graph node to be evicted.
        reset_attr: List used to keep track of the nodes, that after
                    the algorithm terminates need their node attributes reset.
    """

    graph.nodes[node]['evicted'] = 1
    for neighbor in graph[node]:
        if core[neighbor] == k:
            graph.nodes[neighbor]['cd'] -= 1
            reset_attr.add(neighbor)
            if (graph.nodes[neighbor]['cd'] == k and not graph.nodes[neighbor]['evicted']):
                propagate_eviction(graph, k, neighbor, candidates, reset_attr)


def traversal_remove(graph, n, u, v):
    """Removes the edge between the two given graph nodes and updates the necessary core numbers.

    The edge (u, v) in the graph is removed. The core number of the nodes in the graph
    is updated based on their maximum-core degree.
    If a node's maximum core degree is smaller than k after the deletion of
    the edge, its core number is decreased by 1, since the node doesn't have enough
    neighbors to be in the k-core.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer defining up to which n-hop neighborhood should be considered
        while running the algorithm.
        u, v: Integers representing the nodes between which the edge is removed.
    """

    # Initialization.
    root = u
    if core[v] < core[u]:
        root = v
    graph.remove_edge(u, v)
    changed = set()  # Stores all the nodes that have their core number changed.
    reset_attr = set()  # Stores all the nodes that need their node attributes reset at the end.
    # Adjust residential core degrees after the edge of the insertion.
    multihop_prepare_rcd_removal(graph, n, u, v)
    k = core[root]

    # Core Phase.
    if core[u] != core[v]:
        graph.nodes[root]['visited'] = 1
        reset_attr.add(root)
        graph.nodes[root]['cd'] = graph.nodes[root]['rcd'][0]
        if graph.nodes[root]['cd'] < k:  # Node does not have enough suitable neighbors to be in the k-core.
            propagate_dismissal(graph, k, root, changed, reset_attr)
    else:
        graph.nodes[u]['visited'] = 1
        reset_attr.add(u)
        graph.nodes[u]['cd'] = graph.nodes[u]['rcd'][0]
        if graph.nodes[u]['cd'] < k:  # Node does not have enough suitable neighbors to be in the k-core.
            propagate_dismissal(graph, k, u, changed, reset_attr)
        graph.nodes[v]['visited'] = 1
        reset_attr.add(v)
        graph.nodes[v]['cd'] = graph.nodes[v]['rcd'][0]
        # Node does not have enough suitable neighbors to be in the k-core and
        # was not already dismissed by the potential propagate_dismissal call node u.
        if not graph.nodes[v]['dismissed'] and graph.nodes[v]['cd'] < k:
            propagate_dismissal(graph, k, v, changed, reset_attr)
    
    # End Phase.
    for node in reset_attr:
        # Every node attribute besides 'rcd' has to be reset, since they are not
        # initialized at each function call. So without resetting further function calls
        # might result in errors.
        micd = graph.nodes[node]['micd']
        rcd = graph.nodes[node]['rcd']
        graph.nodes[node].clear()
        graph.nodes[node]['pruned'] = []
        graph.nodes[node]['micd'] = micd
        graph.nodes[node]['rcd'] = rcd
    for node in changed:
        for neighbor in graph[node]:
            if core[neighbor] == k:
                graph.nodes[node]['micd'] += 1
            if core[neighbor] == k-1:
                graph.nodes[neighbor]['micd'] -= 1
    multihop_recompute_rcd_removal(graph, n, changed)


def propagate_dismissal(graph, k, node, changed, reset_attr):
    """Reduces the core number of all nodes that cannot be in the k-core after the edge removal.

    The given node is dismissed and has its core number decreased. This dismissal results
    in every neighbor of node losing a neighbor in the k-core themselves. As a result, those
    neighbors might no longer be in the k-core and they are dismissed. This process
    propagates until no further node needs to be dismissed.
    Every node that has its core number decremented is added to the set changed, which
    is needed to recompute the rcd values in the 'End Phase' of traversal_remove.

    Args:
        graph: Graph object created with the networkx library..
        k: Minimum core number of the two nodes between which the edge was deleted.
        node: Integer representing the graph node on which the
            function is called.
        changed: Set of nodes that have there core number reduced already.
        reset_attr: List used to keep track of the nodes, that after
                    the algorithm terminates need their node attributes reset.
    """

    graph.nodes[node]['dismissed'] = 1
    core[node] -= 1
    changed.add(node)
    for neighbor in graph[node]:
        if core[neighbor] == k:
            if not graph.nodes[neighbor]['visited']:
                # The cd of the neighbor might already have been altered in a previous
                # propgate_dismissal call, so addition must be used to take that into account.
                graph.nodes[neighbor]['cd'] += graph.nodes[neighbor]['rcd'][0]
                graph.nodes[neighbor]['visited'] = 1
            reset_attr.add(neighbor)
            # cd of every neighbor is reduced due to the dismissal of the input node.
            graph.nodes[neighbor]['cd'] -= 1
            if graph.nodes[neighbor]['cd'] < k and not graph.nodes[neighbor]['dismissed']:
                propagate_dismissal(graph, k, neighbor, changed, reset_attr)


def multihop_prepare_rcd_insertion(graph, n, u, v):
    """Prepares the rcd values needed in traversal_insert.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer representing up to which n-hop neighborhood should.
        be considered while calculating rcd values.
        u, v: Integers representing graph nodes between which the edge is inserted.
    """

    root = u
    if core[v] < core[u]:
        root = v
    k = core[root]
    frontiers = [set() for index in range(n)]
    if core[u] != core[v]:  # If their core number is equal only rcd changes originating from
        for h in range(n):  # the root occure.
            graph.nodes[root]['rcd'][h] += 1  # The root gained a neighbor with higher core number.
            if h < n-1 and graph.nodes[root]['rcd'][h] == k+1:
                # If the n-core degree of the root increases to k+1, means that it was k before.
                # So by Definition 7 in the paper this causes a change in the rcd values
                # of its neighbors with equal core number.
                frontiers[h+1].add(root)
            if h > 0:
                for node in frontiers[h]:
                    for neighbor in graph[node]:
                        if core[neighbor] == k:
                            graph.nodes[neighbor]['rcd'][h] += 1
                            if h < n-1 and graph.nodes[neighbor]['rcd'][h] == k+1:
                                # If the neighbors rcd value itself is now k+1, its neighbors
                                # rcd values have to be adjusted as well.
                                frontiers[h+1].add(neighbor)
    else:
        for h in range(n):
            # The maximum-core degree has to be adjusted separately since it does not depend
            # on the residential core degree of its neighbors.
            if h == 0:
                graph.nodes[u]['rcd'][h] += 1
                if graph.nodes[u]['rcd'][h] == k+1:
                    frontiers[h+1].add(u)
                graph.nodes[v]['rcd'][h] += 1
                if graph.nodes[v]['rcd'][h] == k+1:
                    frontiers[h+1].add(v)
            else:
                # At each iteration the h-core degree of u and v has to be adjusted since
                # they have equal core number at the point of insertion.
                if graph.nodes[v]['rcd'][h-1] > k:
                    graph.nodes[u]['rcd'][h] += 1
                    if h < n-1 and graph.nodes[u]['rcd'][h] == k+1:
                        frontiers[h+1].add(u)
                if graph.nodes[u]['rcd'][h-1] > k:
                    graph.nodes[v]['rcd'][h] += 1
                    if h < n-1 and graph.nodes[v]['rcd'][h] == k+1:
                        frontiers[h+1].add(v)
                for node in frontiers[h]:
                    # This loop works just like it did in the first case. One thing to consider
                    # is that the residential core degrees of u and v are adjusted before and therefore
                    # have to be excluded.
                    for neighbor in graph[node]:
                        if (not (node == u and neighbor == v) and
                                not (node == v and neighbor == u) and
                                core[neighbor] == k):
                            graph.nodes[neighbor]['rcd'][h] += 1
                            if h < n-1 and graph.nodes[neighbor]['rcd'][h] == k+1:
                                frontiers[h+1].add(neighbor)


def multihop_recompute_rcd_insertion(graph, n, changed):
    """Recomputes the necessary rcd values after traversal_insert.

    The nodes in 'changed' and their neighbors with
    core(neighbor) = core[node] - 1 have their rcd values recomputed
    up to the n-core degree.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer representing up to which n-hop neighborhood should
        be considered while calculating rcd values.
        changed: Set of nodes that had their core number changed in
                traversal_remove.
    """

    for node in changed:
        graph.nodes[node]['visited'] = 1
    for h in range(n):
        updated = set()    # Set of nodes that need their h-core degree recomputed.
        for node in changed:
            # This loop determines the neighbors of the nodes in changed that
            # need their h-core degree recomputed.
            for neighbor in graph[node]:
                if (not graph.nodes[neighbor]['visited'] and
                        (core[neighbor] == core[node] or core[neighbor] == core[node] - 1)):
                    updated.add(neighbor)
                    graph.nodes[neighbor]['visited'] = 1
        for node in updated:
            changed.add(node)
        for node in changed:
            # Computation of the h-core degree at each loop iteration.
            compute_rcd(graph, node, h)
    for node in changed:
        # Node attribute 'visited' has to be reset, since they could result in errors
        # in traversal_remove and traversal_insert otherwise.
        graph.nodes[node]['visited'] = 0


def multihop_prepare_rcd_removal(graph, n, u, v):
    """Prepares the rcd values that are needed in traversal_remove.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer representing up to which n-hop neighborhood should
        be considered while calculating rcd values.
        u, v: Integers representing graph nodes between which the edge is removed.
    """

    root = u
    if core[v] < core[u]:
        root = v
    k = core[root]
    frontiers = [set() for index in range(n)]
    if core[u] != core[v]:  # If their core number is equal only rcd changes originating from
        for h in range(n):  # the root occur.
            graph.nodes[root]['rcd'][h] -= 1
            if h < n-1 and graph.nodes[root]['rcd'][h] == k:
                # If the n-core degree of the root decreases to k, means that it was k+1 before.
                # So by Definition 7 in the paper this causes a change in the rcd values
                # of its neighbors with equal core number.
                frontiers[h+1].add(root)
            if h > 0:
                for node in frontiers[h]:
                    for neighbor in graph[node]:
                        if core[neighbor] == k:
                            graph.nodes[neighbor]['rcd'][h] -= 1
                            if h < n-1 and graph.nodes[neighbor]['rcd'][h] == k:
                                # If the neighbors rcd value itself is now k, its neighbors
                                # rcd values have to be adjusted as well.
                                frontiers[h+1].add(neighbor)
    else:
        # Create a copy of the residential-core degrees of u and v.
        old_rcd_u = [*graph.nodes[u]['rcd']]
        old_rcd_v = [*graph.nodes[v]['rcd']]
        for h in range(n):
            # The maximum-core degree has to be adjusted separately since it does not depend
            # on the residential core degree of its neighbors.
            if h == 0:
                graph.nodes[u]['rcd'][h] -= 1
                if graph.nodes[u]['rcd'][h] == k:
                    frontiers[h+1].add(u)
                graph.nodes[v]['rcd'][h] -= 1
                if graph.nodes[v]['rcd'][h] == k:
                    frontiers[h+1].add(v)
            else:
                # At each iteration the h-core degree of u and v has to be adjusted since
                # they have equal core number at the point of insertion.
                if old_rcd_v[h-1] > k:
                    graph.nodes[u]['rcd'][h] -= 1
                    if h < n-1 and graph.nodes[u]['rcd'][h] == k:
                        frontiers[h+1].add(u)
                if old_rcd_u[h-1] > k:
                    graph.nodes[v]['rcd'][h] -= 1
                    if h < n-1 and graph.nodes[v]['rcd'][h] == k:
                        frontiers[h+1].add(v)
                for node in frontiers[h]:
                    # This loop works just like it did in the first case. One thing to consider
                    # is that the residential core degrees of u and v are adjusted before and therefore
                    # have to be excluded.
                    for neighbor in graph[node]:
                        if (not (node == u and neighbor == v) and
                                not (node == v and neighbor == u) and
                                core[neighbor] == k):
                            graph.nodes[neighbor]['rcd'][h] -= 1
                            if h < n-1 and graph.nodes[neighbor]['rcd'][h] == k:
                                frontiers[h+1].add(neighbor)


def multihop_recompute_rcd_removal(graph, n, changed):
    """Recomputes the necessary rcd values after traversal_remove.

    The nodes in 'changed' and their neighbors with
    core(neighbor) = core[node]+1 have their rcd values recomputed
    up to the n-core degree.

    Args:
        graph: Graph object created with the networkx library.
        n: Integer representing up to which n-hop neighborhood should
        be considered while calculating rcd values.
        changed: Set of nodes that had their core number changed in
                traversal_remove.
    """

    for node in changed:
        graph.nodes[node]['visited'] = 1
    for h in range(n):
        updated = set()  # Set of nodes that need their h-core degree recomputed.
        for node in changed:
            # This loop determines the neighbors of the nodes in changed that
            # need their h-core degree recomputed.
            for neighbor in graph[node]:
                if (not graph.nodes[neighbor]['visited'] and
                        (core[neighbor] == core[node] or core[neighbor] == core[node] + 1)):
                    updated.add(neighbor)
                    graph.nodes[neighbor]['visited'] = 1
        for node in updated:
            changed.add(node)
        for node in changed:
            # Computation of the h-core degree at each loop iteration.
            compute_rcd(graph, node, h)
    # Node attribute 'visited' has to be reset, since they could result in errors
    # in traversal_remove and traversal_insert otherwise.
    for node in changed:
        graph.nodes[node]['visited'] = 0


def micd_init(graph):
    for node in graph:
        graph.nodes[node]['pruned'] = []
        for neighbor in graph[node]:
            if core[neighbor] > core[node]:
                graph.nodes[node]['micd'] += 1


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

    graph = DefaultGraph()  # Uncomment for traversal.
    # graph = nx.Graph()  # Uncomment for order and color.

    with open(file) as f:
        for line in f:
            node_1, node_2 = tuple(line.split())
            graph.add_edge(int(node_1), int(node_2))
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph

# Initialization.

core = nx.core_number(graph)
init_rcd(graph, 2)
micd_init(graph)

# After the graph definition and initialization the code to test the script
# can be inserted below:
