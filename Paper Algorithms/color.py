"""This script provides an implementation of the algorithms specified in the
2012 paper 'Efficient Core Maintenance in Large Dynamic Graphs' by
Rong-Hua Li and Jeffrey Xu Yu.

Imports:
    networkx: Provides functionality for the creation, manipulation, and
            modification of graphs.

    timeit: This module provides a simple way to time Python code.

    random: This module implements pseudo-random number generators for
            various distributions. Used to randomly choose nodes for our
            time performance evaluation.
"""
import networkx as nx
import timeit
import random


def color_insert(graph, u, v):
    """Inserts an edge between two given graph nodes and updates the necessary core numbers.

    Updates the core number of the nodes in a graph after the insertion of
    an edge between the two vertices u and v. The core numbers have to be stored
    in a dictionary or list named 'core'. The graph has to be an instance of Graph(),
    from the networkx library.

    Args:
        graph: Graph object created with the networkx library.
        u, v: Integers representing the nodes between which
            the edge is to be inserted.
    """

    nx.set_node_attributes(graph, 0, 'visited')
    nx.set_node_attributes(graph, 0, 'color')
    candidates = set()
    graph.add_edge(u, v)
    if core[u] > core[v]:
        k = core[v]
        xy_prune_color(graph, v, v, k, candidates)
        recolor_insert(graph, k, candidates)
        update_insert(graph, k, candidates)
    else:
        k = core[u]
        xy_prune_color(graph, u, u, k, candidates)
        recolor_insert(graph, k, candidates)
        update_insert(graph, k, candidates)


def xy_prune_color(graph, node, root, k, candidates):
    """Prunes nodes that cannot have their core number changed.

    Fills the set 'candidates' with every node in the graph,
    that may have its core number increased after the edge insertion from
    color_insert.

    Args:
        graph: Graph object created with the networkx library.
        node: Integer representing the currently processed node.
        root: Integer representing the node from which the pruning process originates.
        k: Minimum core number of the two nodes between which an edge was inserted.
        candidates: Set of integers representing the nodes that possibly
                    have their core number changed.
    """

    graph.nodes[node]['visited'] = 1
    graph.nodes[node]['y_node'] = 0
    graph.nodes[node]['x_node'] = 0
    for neighbor in graph[node]:
        if core[neighbor] >= k:  # Computation of x_node.
            graph.nodes[node]['x_node'] += 1
        if node is not root and core[neighbor] > k:  # Computation of y_node.
            graph.nodes[node]['y_node'] += 1
    if graph.nodes[node]['x_node'] > k:  # X-pruning.
        if graph.nodes[node]['y_node'] < k or k == 0:  # Y-pruning.
            for neighbor in graph[node]:
                if (graph.nodes[neighbor]['visited'] == 0 and core[neighbor] == k):
                    xy_prune_color(graph, neighbor, node, k, candidates)
        if graph.nodes[node]['color'] == 0:
            candidates.add(node)
            graph.nodes[node]['color'] = 1


def recolor_insert(graph, k, candidates):
    """Determines the candidates that definitely need their core number updated.

    Recursively determines which nodes in 'candidates' cannot increase their core
    number after the edge insertion into the graph.

    Args:
        graph: Graph object created with the networkx library.
        k: Minimum core number of the two nodes between which the edge was inserted.
        candidates: Set of integers representing the nodes that possibly
                    have their core number changed.
    """

    flag = 0
    for node in candidates:
        if graph.nodes[node]['color'] == 1:
            graph.nodes[node]['x_node'] = 0
            for neighbor in graph[node]:  # Computation of modified x_node.
                if (graph.nodes[neighbor]['color'] == 1 or core[neighbor] > k):
                    graph.nodes[node]['x_node'] += 1
            if graph.nodes[node]['x_node'] <= k:
                graph.nodes[node]['color'] = 0
                flag = 1
    if flag == 1:
        recolor_insert(graph, k, candidates)


def update_insert(graph, k, candidates):
    """Increments the core number of every colored node in the set 'candidates'."""
    for node in candidates:
        if graph.nodes[node]['color'] == 1:
            core[node] += 1


def color_remove(graph, u, v):
    """Removes the edge between the two given graph nodes and updates the necessary core numbers.

    Updates the core number of the nodes in a graph after the removal
    of an edge between the two vertices u and v. The core numbers have to be stored
    in a dictionary or list named 'core'. The graph has to be an instance of Graph(),
    from the networkx library.

    Args:
        graph: Graph object created with the networkx library.
        u, v: Integers representing the graph nodes between which
            the edge is to be removed.
    """

    nx.set_node_attributes(graph, 0, 'visited')
    nx.set_node_attributes(graph, 0, 'color')
    candidates = set()
    graph.remove_edge(u, v)
    graph.nodes[u]['x_node'] = len([neighbor for neighbor in graph[u] if core[neighbor] >= core[u]])
    graph.nodes[v]['x_node'] = len([neighbor for neighbor in graph[v] if core[neighbor] >= core[v]])
    if core[u] > core[v]:
        k = core[v]
        if graph.nodes[v]['x_node'] < k:  # X-pruning.
            y_prune_color(graph, v, k, candidates)
            recolor_remove(graph, k, candidates)
            update_remove(graph, k, candidates)
    if core[u] < core[v]:
        k = core[u]
        if graph.nodes[u]['x_node'] < k:  # X-pruning.
            y_prune_color(graph, u, k, candidates)
            recolor_remove(graph, k, candidates)
            update_remove(graph, k, candidates)
    if core[v] == core[u]:
        k = core[u]
        if (graph.nodes[u]['x_node'] < k and graph.nodes[v]['x_node'] < k):  # X-pruning.
            y_prune_color(graph, u, k, candidates)
            if graph.nodes[v]['color'] == 0:
                # If the subcores of u and v are disconnected y_prune_color has to be
                # called on both of them.
                nx.set_node_attributes(graph, 0, 'visited')
                y_prune_color(graph, v, k, candidates)
                recolor_remove(graph, k, candidates)
                update_remove(graph, k, candidates)
            else:
                recolor_remove(graph, k, candidates)
                update_remove(graph, k, candidates)
        if graph.nodes[u]['x_node'] < k and graph.nodes[v]['x_node'] >= k:  # X-pruning.
            y_prune_color(graph, u, k, candidates)
            recolor_remove(graph, k, candidates)
            update_remove(graph, k, candidates)
        if (graph.nodes[u]['x_node'] >= k and graph.nodes[v]['x_node'] < k):  # X-pruning.
            y_prune_color(graph, v, k, candidates)
            recolor_remove(graph, k, candidates)
            update_remove(graph, k, candidates)


def y_prune_color(graph, node, k, candidates):
    """Prunes nodes that cannot have their core number increased after the edge removal.

    Fills the set 'candidates' with nodes that may have their core number decreased and
    applies Y-pruning to prune the vertices that cannot have their core number decreased.

    Args:
        graph: Graph object created with the networkx library.
        node: Integer representing the currently processed node.
        k: Minimum core number of the two nodes between which the edge was deleted.
        candidates: Set of integers representing the nodes that possibly
                    have their core number changed.
    """

    graph.nodes[node]['visited'] = 1
    if graph.nodes[node]['color'] == 0:
        candidates.add(node)
        graph.nodes[node]['color'] = 1
    # Computation of y_node.
    graph.nodes[node]['y_node'] = len([neighbor for neighbor in graph[node] if core[neighbor] > k])
    if graph.nodes[node]['y_node'] < k:  # Y-pruning.
        for neighbor in graph[node]:
            if (graph.nodes[neighbor]['visited'] == 0 and core[neighbor] == k):
                y_prune_color(graph, neighbor, k, candidates)


def recolor_remove(graph, k, candidates):
    """Determines the candidates that definitely need their core number updated.

    Recursively determines which nodes in 'candidates' cannot decrease their core
    number after the edge removal in colore_remove.

    Args:
        graph: Graph object created with the networkx library.
        k: Minimum core number of the two nodes between which the edge was deleted.
        candidates: Set of integers representing the nodes that possibly
                    have their core number changed.
    """

    flag = 0
    for node in candidates:
        if graph.nodes[node]['color'] == 1:
            graph.nodes[node]['x_node'] = 0
            for neighbor in graph[node]:  # Computation of modified x_node.
                if (graph.nodes[neighbor]['color'] == 1 or core[neighbor] > k):
                    graph.nodes[node]['x_node'] += 1
            if graph.nodes[node]['x_node'] < k:
                graph.nodes[node]['color'] = 0
                flag = 1
    if flag == 1:
        recolor_remove(graph, k, candidates)


def update_remove(graph, k, candidates):
    """Decrements the core number of every node in the set 'candidates'."""
    for node in candidates:
        if graph.nodes[node]['color'] == 0:
            core[node] = core[node] - 1


# Insert the graph definition here:

# Initialization.
core = nx.core_number(graph)

# After the graph definition and initialization the code to test the script
# can be inserted here:
