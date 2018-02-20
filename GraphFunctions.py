"""Graphs.py

Various simple functions for graph input.

Each function's input graph G should be represented in such a way that "for v in G"
loops through the vertices, and "G[v]" produces a list of the neighbors of v;
for instance, G may be a dictionary mapping each vertex to its neighbor set.

D. Eppstein, April 2004.
"""
from typing import Set

from Models import Graph, DirectedGraph, Vertex, UndirectedGraph


def is_undirected(graph: DirectedGraph):
    """Check that the graph represents a simple undirected graph."""
    for v in graph:
        if v in graph[v]:
            return False
        for w in graph[v]:
            if v not in graph[w]:
                return False
    return True


def max_out_degree(graph: DirectedGraph):
    """Return the maximum vertex out-degree of graph."""
    return max([len(graph.outgoing_links_of(v)) for v in graph])


def min_out_degree(graph: DirectedGraph):
    """Return the minimum vertex out-degree of graph."""
    return min([len(graph.outgoing_links_of(v)) for v in graph])


def max_in_degree(graph: DirectedGraph):
    """Return the maximum vertex in-degree of graph."""
    return max([len(graph.ingoing_links_of(v)) for v in graph])


def min_in_degree(graph: DirectedGraph):
    """Return the minimum vertex in-degree of graph."""
    return min([len(graph.ingoing_links_of(v)) for v in graph])


def copy_graph(graph: DirectedGraph):
    """
    Make a copy of a graph G and return the copy.
    Any information stored in edges G[v][w] is discarded.
    
    Most of the time, copy.deepcopy will be preferable to this function;
    however, unlike deepcopy, this function can change the data type
    of the adjacency list of the given graph.

    The second argument should be a callable that turns a sequence
    of neighbors into an appropriate representation of the adjacency list.
    Note that, while Set, list, and tuple are appropriate values for
    adjacency_list_type, dict is not -- use Util.map_to_constant instead.
    """
    return DirectedGraph({v: set(graph[v]) for v in graph})


def induced_subgraph(vertex_subset: Set[Vertex], original_graph: Graph):
    """
    The subgraph consisting of all edges between pairs of vertices in the vertex subset.
    """
    if original_graph is DirectedGraph:
        graph_builder_class = DirectedGraph
    else:
        graph_builder_class = UndirectedGraph
    return graph_builder_class({x: set(y for y in original_graph[x] if y in vertex_subset)
                                for x in original_graph if x in vertex_subset})


def union(*graphs: DirectedGraph):
    """Return a graph having all edges from the argument graphs."""
    output = {}
    for graph in graphs:
        for vertex in graph:
            output.setdefault(vertex, set()).update(list(graph[vertex]))
    return DirectedGraph(output)


def check_independent_set(vertex_subset: Set[Vertex], graph: DirectedGraph):
    """
    True if V is an independent set of vertices in G, False otherwise.
    """
    for vertex in vertex_subset:
        for other_vertex in graph[vertex]:
            if other_vertex in vertex_subset:
                return False
    return True

