# ------------------------------------------------------------------------------------
"""
This is just some library syntax so you don't have to look that up.
The networkx library implementation is based on adjacency list representation
and is implemented with python dictionaries.
So it looks something like this:

        {node_1 {neighbor_1{}, neighbor_2{}}, node_2{neighbor_1{}, neighbor_2{}}}

We also take advantage of the fact that each node has a dictionary that
holds its node attributes. So that is where we store node properties like
remaining degree, maximum core degree, visited.
To access the node attributes you will see the following in the code a lot:

        graph.nodes[<node>][<node_attribute>]

The following syntax, which you will see a lot as well, allows us to iterate over
the neighbors of a node:

        graph[<node>]

Maybe an important thing to note is that in all the implemented algorithms the
graph nodes are assumed to be simple integers.
"""
# ------------------------------------------------------------------------------------
"""
This function creates a networkx graph based on a given text file.
If you want to run the algorithms on a graph given through one of the text
files just insert this to the end of the respective script, where it says
to insert the graph definition, specify the path and the graph will be created.
Depending on what algorithm you run, you have to uncomment
the right graph initialization at the beginning of the function.

Otherwise, you can just create your own graph like in the example graphs
further down.
"""

'''
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
    # graph = nx.Graph()  # Uncomment for order and color.

    with open(file) as f:
        for line in f:
            node_1, node_2 = tuple(line.split())
            graph.add_edge(int(node_1), int(node_2))
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph

graph = file_to_graph('<path_to_txt_file>')
'''
# ------------------------------------------------------------------------------------
"""
Following code tests whether the specified function updated the core numbers correctly.
If you would like to test insertion and removal functions simultaneously
you can just add the removal function below the insertion one in 'Test insertion'.
"""
'''
# Test insertion:
# When testing traversal_insert function parameter n has to added.
graph_nodes = list(nx.nodes(graph))
node_1 = random.choice(list(graph_nodes))
node_2 = random.choice(list(graph_nodes))
if node_1 != node_2 and (not (node_1 in graph[node_2])):

    <function_to_test>(graph, node_1, node_2)

core_2 = nx.core_number(graph)
result = []
for node in graph:
    result.append(core[node] == core_2[node])
print(all(result))


# Test removal:
# When testing traversal_remove function parameter n has to added.
graph_edges = list(graph.edges())
random_edge = random.choice(graph_edges)
graph_edges.remove(random_edge)
node_1, node_2 = random_edge

<function_to_test>(graph, node_1, node_2)

core_2 = nx.core_number(graph)
result = []
for node in graph:
    result.append(core[node] == core_2[node])
print(all(result))
'''
# ------------------------------------------------------------------------------------
# Writes the time taken of any function to a text file.
'''
with open("<filename>.txt", 'w') as f:
    f.write("\tInsertion\t\tDeletion \n")
    for _ in range(<how often>):
        node_1 = random.choice(list(nx.nodes(graph)))
        node_2 = random.choice(list(nx.nodes(graph)))
        if node_1 != node_2 and (not (node_1 in graph[node_2])):
            insert = timeit.timeit("""<function_1>(<parameters_1>)""", globals=globals(), number=1)
            remove = timeit.timeit("""<function_2>(<parameters_2>)""", globals=globals(), number=1)
            f.write(f'{_}\t{insert:.6f}\t\t\t{remove:.6f} \n')
'''
# ------------------------------------------------------------------------------------
# Figure 1 from the 2012 paper.
'''
# graph = DefaultGraph()  # when used in traversal.py
graph = nx.Graph()
edges = [(1, 2), (1, 3), (2, 3), (2, 4), (2, 9), (3, 4), (3, 5), (3, 7),
         (3, 6), (4, 7), (4, 6), (4, 5), (5, 6), (5, 7), (5, 8), (6, 7),
         (6, 8), (7, 9), (9, 10), (10, 11), (10, 18), (11, 12),
         (11, 13), (12, 13), (13, 14), (13, 17), (14, 15), (14, 16),
         (14, 17), (15, 16), (15, 17), (16, 17), (18, 16), (18, 17)]
graph.add_edges_from(edges)
'''
# ------------------------------------------------------------------------------------
# Figure 3 from the 2013 paper.
'''
# graph = DefaultGraph()  # when used in traversal.py
graph = nx.Graph()
edges = [(1, 2), (1, 3), (1, 9), (2, 3), (2, 9), (2, 4), (2, 10), (2, 11),
         (2, 12), (3, 9), (3, 4), (3, 5), (3, 6), (3, 10), (3, 11), (3, 12),
         (4, 5), (4, 6), (4, 9), (4, 10), (4, 11), (4, 12),
         (5, 6), (5, 9), (6, 9), (7, 8), (7, 9), (8, 9),
         (10, 11), (10, 12), (11, 12)]
graph.add_edges_from(edges)
'''
# ------------------------------------------------------------------------------------
# Figure 2 from the 2013 Paper
'''
# graph = DefaultGraph()  # when used in traversal.py
graph = nx.Graph()
edges = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5),
         (5, 6), (5, 14), (6, 7), (7, 8), (8, 9), (8, 15), (9, 10), (9, 15),
         (10, 11), (10, 16), (11, 12), (11, 18), (11, 19), (12, 16), (12, 13),
         (13, 15), (13, 14), (14, 15), (17, 18), (17, 19), (17, 20),
         (18, 19), (18, 20), (19, 20)]
graph.add_edges_from(edges)
'''
# ------------------------------------------------------------------------------------
# Figure 3 from the 2016 paper.
# Shortened the bottom half, but still illustrates what it is
# supposed to show.
'''
# graph = DefaultGraph()  # when used in traversal.py
graph = nx.Graph()
edges = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5),
         (5, 6), (5, 7), (5, 8), (5, 9), (5, 13), (6, 7), (7, 10),
         (8, 9), (8, 14), (10, 11), (10, 12), (10, 13), (11, 12),
         (11, 13), (12, 13), (14, 15), (14, 20), (15, 16), (16, 17),
         (17, 18), (18, 19), (20, 21), (21, 22), (22, 23), (23, 24)]
graph.add_edges_from(edges)
'''
# ------------------------------------------------------------------------------------
