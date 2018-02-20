"""
The Graph class is intended to generify how graphs are used
in the PADS classes. The Tree class is similar.
The reasoning for these: some of the classes and functions are
expecting a dictionary of children - for example, when given
some G, G["x"] = ["a", "b", "c"]. Some other times a mapping
of parents is expected, which could be G["a"]="x", G["b"]="x"...
I decided that I don't like this. Plus, some people want to
work with object oriented stuff, and not arrays/lists/dicts.
"""
from collections import defaultdict
from typing import TypeVar, Generic, Optional, Set, Dict, DefaultDict

Vertex = TypeVar("Vertex")


class Tree(Generic(Vertex)):
    def __init__(self, children_of: Dict[Vertex, Set[Vertex]] = None, parent_of: Dict[Vertex, Vertex] = None):
        if children_of is not None and parent_of is not None:
            raise ValueError("Cannot create a tree from multiple sources.")
        if children_of is None and parent_of is None:
            raise ValueError("Cannot create a tree from zero sources.")

        self._children_of: DefaultDict[Vertex, Set[Vertex]] = defaultdict(list)
        self._parent_of: DefaultDict[Vertex, Vertex] = defaultdict()

        if children_of:
            self._children_of = children_of
            for vertex in children_of.keys():
                for child in children_of[vertex]:
                    self._parent_of[child] = vertex

        elif parent_of:
            self._parent_of = parent_of
            for parent, child in parent_of.keys():
                self._children_of[parent].append(child)

    def set_relation(self, parent: Vertex, child: Vertex):
        """If child already had a parent, this will remove the previous link"""
        if child in self._parent_of:
            self._children_of[self._parent_of[child]].remove(child)
        self._parent_of[child] = parent
        self._children_of[parent].add(child)

    def remove_vertex(self, vertex: Vertex):
        """Will be slightly faster if you remove subtrees top-to-bottom"""
        if vertex in self._parent_of:
            self._children_of[self._parent_of[vertex]].remove(vertex)
            del self._parent_of[vertex]
        if vertex in self._children_of:
            for child in self._children_of[vertex]:
                self._parent_of[child] = None
            del self._children_of[vertex]

    def parent_of(self, vertex: Vertex) -> Optional[Vertex]:
        return self._parent_of[vertex]

    def children_of(self, vertex: Vertex) -> Set[Vertex]:
        return self._children_of[vertex]


class Graph:
    def outgoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        pass

    def ingoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        pass

    def all_vertices(self) -> Set[Vertex]:
        pass

    def __iter__(self) -> Set[Vertex]:
        return self.all_vertices()

    def __getitem__(self, item: Vertex) -> Set[Vertex]:
        return self.outgoing_links_of(item)


class DirectedGraph(Graph):
    def __init__(self, outgoing_links_of: Dict[Vertex, Set[Vertex]]):
        self._outgoing_links_of = outgoing_links_of
        self._ingoing_links_of: DefaultDict[Vertex, Set[Vertex]] = defaultdict(list)
        self._all_vertices = set(outgoing_links_of.keys())
        for vertex in outgoing_links_of:
            for out_neighbor in outgoing_links_of[vertex]:
                self._ingoing_links_of[out_neighbor].append(vertex)

    def outgoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        return self._outgoing_links_of[vertex]

    def ingoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        return self._ingoing_links_of[vertex]

    def all_vertices(self) -> Set[Vertex]:
        return self._all_vertices


class UndirectedGraph(Graph):
    def __init__(self, neighbors_of: Dict[Vertex, Set[Vertex]]):
        self._neighbors_of = neighbors_of
        self._all_vertices = set(neighbors_of.keys())

    def neighbors_of(self, vertex: Vertex) -> Set[Vertex]:
        return self._neighbors_of[vertex]

    def outgoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        return self.outgoing_links_of(vertex)

    def ingoing_links_of(self, vertex: Vertex) -> Set[Vertex]:
        return self.outgoing_links_of(vertex)

    def all_vertices(self) -> Set[Vertex]:
        return self._all_vertices
