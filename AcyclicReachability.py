"""
Bit-parallel algorithm for testing which vertices can reach which other vertices
 in a DAG (Directed Acyclic Graph).

Usage:
    R = Reachability(G)
    ...
    R.reachable(source,destination)
returns a boolean value: True if G contains a path from the source vertex
to the destination vertex, and False otherwise. The initialization of R
performs a linear number of bit-vector operations, after which each reachability
test takes constant time to perform.

D. Eppstein, April 2009.
"""

import unittest
from typing import List, Dict, TypeVar, Generic

from PartialOrder import TopologicalOrder

T = TypeVar("T")


class Reachability(Generic[T]):
    def __init__(self, graph: Dict[T, List[T]]):
        """Initialize a reachability data structure for the given DAG."""
        self._key_of = {}
        self._can_reach = []
        sorted_vertices = TopologicalOrder(graph)
        sorted_vertices.reverse()
        for vertex in sorted_vertices:
            key = self._key_of[vertex] = len(self._can_reach)
            bits = 1 << key
            for other_vertex in graph[vertex]:
                bits |= self._can_reach[self._key_of[other_vertex]]
            self._can_reach.append(bits)

    def reachable(self, source: T, destination: T):
        """Test whether the DAG has a path from source to destination."""
        return (1 << self._key_of[destination]) & self._can_reach[self._key_of[source]] != 0


class ReachabilityTest(unittest.TestCase):
    def testReachable(self):
        G = {"A": ["C"], "B": ["C", "D"], "C": ["D", "E"], "D": [], "E": []}
        R = Reachability(G)
        for s in "ABCDE":
            for t in "ABCDE":
                self.assertEqual(R.reachable(s, t),
                                 s <= t and s + t not in ["AB", "DE"])


if __name__ == "__main__":
    unittest.main()
