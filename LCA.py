"""LCA.py

Range minimization and tree least common ancestor data structures
with linear space and preprocessing time, and constant query time,
from Bender and Farach-Colton, "The LCA Problem Revisited",
Proc. LATIN 2000 (pp.88-94), http://www.cs.sunysb.edu/~bender/pub/lca.ps

Some experimentation would be needed to determine how large a query
range needs to be to make this faster than computing the min of the range
directly, and how much input data is needed to make the linear space
version pay off compared to the much simpler LogarithmicRangeMin that
it uses as a subroutine.

D. Eppstein, November 2003.

Some experimentation has determined that LogarithmicRangeMin is actually
faster than the linear-space version, in all of my test cases. Sorry.

I've also edited the code a bunch to update it to Python 3, make it
cleaner, and choose the Log version as the default.

Itamar Curiel, February 2018
"""

import random
import unittest
from collections import defaultdict
from typing import Sized, List, Tuple, Dict

from UnionFind import UnionFind


def _xenumerate(iterable, reverse=False):
    """Return pairs (x,i) for x in iterable, where i is
    the index of x in the iterable.
    It's like enumerate(iterable), but each element is
    (value, index) and not (index, value). This is
    useful for sorting.
    """
    if reverse:
        indices = range(len(iterable) - 1, -1, -1)
    else:
        indices = range(len(iterable))
    return [(iterable[i], i) for i in indices]


class SliceableRange(Sized):
    def __init__(self, data: List):
        self._data = list(data)

    def __len__(self) -> int:
        return len(self._data)

    def get_slice_edges(self, interval: slice):
        if type(interval) != slice:
            raise ValueError("Can only access LCA object by slice notation.")
        start, stop, step = interval.indices(len(self))
        if step != 1:
            raise ValueError("Step/Stride not permitted in LCA.")
        return start, stop


class RangeMin(SliceableRange):
    """An offline structure that gives range minimums fast.
    If X is any list:
        RangeMin(X)[i:j] == min(X[i:j]).
    Initializing RangeMin(X) takes time and space linear in len(X),
    and querying the minimum of a range takes constant time per query.
    """

    def __init__(self, data: List):
        """Uses an LCA structure on a Cartesian tree for the input."""
        super().__init__(data)
        if len(self) <= 1:
            self._lca_function = lambda: min(self._data)
        else:
            big = [max(left, right)
                   for left, right in zip(
                    self._smallest_nearest_valuea(side_is_left=False),
                    self._smallest_nearest_valuea(side_is_left=True)
                )]
            parents = {i: big[i][1] for i in range(len(self)) if big[i]}
            self._lca_function = LCA(parents)

    def __getitem__(self, interval: slice):
        """When called by X[left:right], return min(X[left:right])."""
        left, right = self.get_slice_edges(interval)
        if right <= left:
            return None  # empty range has no minimum
        return self._data[self._lca_function(left, right - 1)]

    def _smallest_nearest_valuea(self, *, side_is_left) -> List[Tuple[int, int]]:
        """All nearest smaller values.
        For each x in the data, find the value smaller than x in the closest
        position to the left of x (if not reversed) or to the right of x
        (if reversed), and return list of pairs (smaller value,position).
        Due to our use of positions as a tie-breaker, values equal to x
        count as smaller on the left and larger on the right.
        """
        stack = [(min(self._data), -1)]  # protect stack top with sentinel
        output = [(0, -1) for _ in range(len(self._data))]
        for xi in _xenumerate(self._data, side_is_left):
            while stack[-1] > xi:
                stack.pop()
            output[xi[1]] = stack[-1]
            stack.append(xi)
        return output


class RestrictedRangeMin(SliceableRange):
    """Linear-space RangeMin for integer data obeying the constraint:
        abs(X[i]-X[i-1])==1.
    Meaning: each value is equal to +1 or -1 from the previous value.
    We don't actually check this constraint, so results may be incorrect
    if it is violated.  For the use of this data structure from LCA, the
    data is actually pairs rather than integers, but the minimums of
    all ranges are in the same positions as the minimums of the integers
    in the first positions of each pair, so the data structure still works.
    """

    def __init__(self, data: List):
        # Compute parameters for partition into blocks.
        # Position i in X becomes transformed into:
        #     position i & self._block_mask in block i >> self.block_len
        super().__init__(data)
        self._block_size = _log2(len(data)) // 2
        self._block_mask = (1 << self._block_size) - 1
        block_len = 1 << self._block_size

        # Partition data into blocks, find minimums within
        # each block, prefix minimums in each block,
        # and suffix minimums in each block.
        block_ids = []  # map block to block id
        ids_to_range_mins = {}  # map block index to PrecomputedRangeMin
        block_minimums = []  # map block index to min value
        self._prefix = [-999]  # map data index to prefix min of block, first value unused
        self._suffix = []  # map data index to suffix min of block
        for i in range(0, len(data), block_len):
            block = data[i:i + block_len]
            block_minimums.append(min(block))
            self._prefix += _prefix_minimums(block)
            self._suffix += _prefix_minimums(block, reverse=True)
            if len(block) < block_len:
                block_id = -1
            else:
                block_id = RestrictedRangeMin._calc_block_fingerprint(block)
            block_ids.append(block_id)
            if block_id not in ids_to_range_mins:
                ids_to_range_mins[block_id] = PrecomputedRangeMin(_xenumerate(block))
        self._blocks = [ids_to_range_mins[b] for b in block_ids]

        # Build data structure for inter-block queries
        self._block_range = LogarithmicRangeMin(block_minimums)

    def __getitem__(self, interval):
        """When called by X[left:right], return min(data[left:right])."""
        left, right = self.get_slice_edges(interval)
        first_block = left >> self._block_size
        last_block = (right - 1) >> self._block_size
        if first_block == last_block:
            i = left & self._block_mask
            position = self._blocks[first_block][i:i + right - left][1]
            return self._data[position + (first_block << self._block_size)]
        else:
            best = min(self._suffix[left], self._prefix[right])
            if last_block > first_block + 1:
                best = min(best, self._block_range[first_block + 1:last_block])
            return best

    @staticmethod
    def _calc_block_fingerprint(block):
        """Calculates value such that all blocks with the same
        pattern of increments and decrements get the same result.
        Does so with a binary number made of "flags" for these.
        """
        block_id = 0
        for i in range(1, len(block)):
            block_id <<= 1
            block_id += (block[i] > block[i - 1])
        return block_id


class PrecomputedRangeMin(SliceableRange):
    """RangeMin solved in quadratic space by precomputing all solutions."""

    def __init__(self, data):
        super().__init__(data)
        self._minimums = [_prefix_minimums(data[i:]) for i in range(len(data))]

    def __getitem__(self, interval):
        """When called by X[left:right], return min(X[left:right])."""
        left, right = self.get_slice_edges(interval)
        return self._minimums[left][right - left - 1]

    def __len__(self):
        return len(self._minimums)


class LogarithmicRangeMin(SliceableRange):
    """RangeMin in O(n log n) space and constant query time."""

    def __init__(self, data: List):
        """Computes min(data[i:i+2**j]) for each possible i,j."""
        super().__init__(data)
        self._minimums = m = [list(data)]
        for j in range(_log2(len(data))):
            m.append(list(map(min, m[-1][:-1 << j], m[-1][1 << j:])))

    def __getitem__(self, interval):
        """When called by X[left:right], return min(data[left:right])."""
        left, right = self.get_slice_edges(interval)
        j = _integer_log_table[right - left]
        row = self._minimums[j]
        return min(row[left], row[right - 2 ** j])

    def __len__(self):
        return len(self._minimums[0])


class LCA:
    """Structure for finding least common ancestors in trees.
    Tree nodes may be any hashable objects; a tree is specified
    by a dictionary mapping nodes to their parents.
    LCA(T)(x,y) finds the LCA of nodes x and y in tree T.
    """

    def __init__(self, parent_of: Dict, range_min_factory=LogarithmicRangeMin):
        """Construct LCA structure from tree parent relation.
        The 'parents' dict should map each node to its parent node."""
        children_of = defaultdict(list)
        for x in parent_of:
            children_of[parent_of[x]].append(x)
        roots = [x for x in children_of if x not in parent_of]
        if len(roots) != 1:
            raise ValueError("LCA input is made of multiple trees")

        self._representatives = {}
        """self._representatives[x] will be the index in the
        euler-traversal in which the node is last seen."""
        levels = self._euler_traverse(children_of, roots[0])
        if [x for x in parent_of if x not in self._representatives]:
            raise ValueError("LCA input is not a tree")
        self._range_min = range_min_factory(levels)

    def get_lowest_common_ancestor_of(self, *nodes):
        """Find least common ancestor of a set of nodes."""
        r = [self._representatives[x] for x in nodes]
        return self._range_min[min(r):max(r) + 1][1]

    def __call__(self, *nodes):
        """calls self.get_lowest_common_ancestor_of(...)"""
        return self.get_lowest_common_ancestor_of(*nodes)

    def _euler_traverse(self, children, root) -> List[Tuple[int, int]]:
        """Perform Euler traversal of tree,
        updating self._representatives and returning levels"""
        levels = []
        stack = [(0, root, True)]
        while stack:
            current_level, node, is_first_arrival = stack.pop()
            if is_first_arrival:
                self._representatives[node] = len(levels)
                stack.append((current_level, node, False))
                for child in children[node]:
                    stack.append((current_level + 1, child, True))
                    stack.append((current_level, node, False))
            else:
                levels.append((current_level, node))

        return levels


class OfflineLCA(defaultdict):
    """Find LCAs of all pairs in a given sequence, using Union-Find."""

    def __init__(self, parent, pairs):
        """Set up to find LCAs of pairs in tree defined by parent.
        LCA of any given pair x,y can then be found by self[x][y].
        However unlike the online LCA structure we can not find LCAs
        of pairs that are not supplied to us at startup time.
        """

        # set up dictionary where answers get stored
        # noinspection PyTypeChecker
        defaultdict.__init__(self, dict)
        for u, v in pairs:
            self[u][v] = self[v][u] = None

        # set up data structure for finding node ancestor on search path
        # self.descendants forms a collection of disjoint sets,
        #    one set for the descendants of each search path node.
        # self.ancestors maps disjoint set ids to the ancestors themselves.
        self.descendants = UnionFind()
        self.ancestors = {}

        # invert the parent relationship so we can traverse the tree
        self.children = defaultdict(list)
        for x, px in parent.items():
            self.children[px].append(x)
        root = [x for x in self.children if x not in parent]
        if len(root) != 1:
            raise ValueError("LCA input is not a tree")

        # initiate depth first traversal
        self.visited = set()
        self.traverse(root[0])

    def traverse(self, node):
        """Perform depth first traversal of tree."""
        self.ancestors[self.descendants[node]] = node
        for child in self.children[node]:
            self.traverse(child)
            self.descendants.union(child, node)
            self.ancestors[self.descendants[node]] = node
        self.visited.add(node)
        for query in self[node]:
            if query in self.visited:
                lca = self.ancestors[self.descendants[query]]
                self[node][query] = self[query][node] = lca


def _prefix_minimums(iterable, reverse=False):
    """Compute table of prefix minimums
    (or suffix minimums, if reverse=True) of a list.

    For example, given the list:
        [6, 4, 8, 6, 1, 5, 2]
    prefix minimums:
        [6, 4, 4, 4, 1, 1, 1]
    suffix minimums:
        [1, 1, 1, 1, 1, 2, 2]
    """
    current = None
    output = [None for _ in iterable]
    for x, i in _xenumerate(iterable, reverse):
        if current is None:
            current = x
        else:
            current = min(current, x)
        output[i] = current
    return output


_integer_log_table = [None, 0]


def _log2(n):
    """Make table of logs reach up to n and return floor(log_2(n))."""
    while len(_integer_log_table) <= n:
        _integer_log_table.extend([1 + _integer_log_table[-1]] * len(_integer_log_table))
    return _integer_log_table[n]


# if run as "python LCA.py", run tests on random data
# and check that RangeMin's results are correct.

class RandomRangeMinTest(unittest.TestCase):
    def testRangeMin(self):
        for trial in range(20):
            data = [random.choice(range(1000000))
                    for _ in range(random.randint(1, 100))]
            range_minimum_of = RangeMin(data)
            for sample in range(100):
                i = random.randint(0, len(data) - 1)
                j = random.randint(i + 1, len(data))
                self.assertEqual(range_minimum_of[i:j], min(data[i:j]))


class LCATest(unittest.TestCase):
    parent_of = {'b': 'a', 'c': 'a', 'd': 'a', 'e': 'b', 'f': 'b', 'g': 'f', 'h': 'g', 'i': 'g'}
    true_LCA_answers = {
        ('a', 'b'): 'a',
        ('b', 'c'): 'a',
        ('c', 'd'): 'a',
        ('d', 'e'): 'a',
        ('e', 'f'): 'b',
        ('e', 'g'): 'b',
        ('e', 'h'): 'b',
        ('c', 'i'): 'a',
        ('a', 'i'): 'a',
        ('f', 'i'): 'f',
    }

    def testLCA(self):
        lca_of = LCA(self.parent_of)
        for k, v in self.true_LCA_answers.items():
            self.assertEqual(lca_of(*k), v)

    def testLogLCA(self):
        # noinspection PyTypeChecker
        lca_of = LCA(self.parent_of, range_min_factory=LogarithmicRangeMin)
        for k, v in self.true_LCA_answers.items():
            self.assertEqual(lca_of(*k), v)

    def testOfflineLCA(self):
        lca_of = OfflineLCA(self.parent_of, self.true_LCA_answers.keys())
        for (p, q), v in self.true_LCA_answers.items():
            self.assertEqual(lca_of[p][q], v)


class testSpeeds(unittest.TestCase):
    def testSpeeds(self):
        test_sample_size = 2 ** 15

        objects = [num for num in range(test_sample_size)]
        parent_of = {}
        for i in range(1, 17):
            parent_of[i] = random.randint(0, i - 1)
        for i in range(17, len(objects)):
            parent_of[i] = random.randint(i - 8, i - 5)

        import time
        t0 = time.clock()

        lca_of = LCA(parent_of, range_min_factory=RestrictedRangeMin)
        for i in range(1, len(objects) - 1):
            lca_of(i, i + 1)
        t1 = time.clock()

        print(f"NORMAL - size {test_sample_size} took {t1 - t0} s")

        t0 = time.clock()

        # noinspection PyTypeChecker
        lca_of = LCA(parent_of, range_min_factory=LogarithmicRangeMin)
        for i in range(1, len(objects) - 1):
            lca_of(i, i + 1)
        t1 = time.clock()

        print(f"LOG - size {test_sample_size} took {t1 - t0} s")


if __name__ == "__main__":
    unittest.main()
