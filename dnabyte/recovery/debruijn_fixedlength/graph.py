from typing import Dict, List
from collections import defaultdict

class DeBruijnGraphGraph:
    """
    Weighted De Bruijn graph.

    Nodes are (k-1)-mers.
    Directed edges connect overlapping nodes and store how many reads
    support that transition.
    """

    def __init__(self, kmer_size, min_coverage, branch_ratio, logger=None):

        self.logger = logger
        
        self.logger = logger
        
        self.min_coverage = min_coverage
        self.branch_ratio = branch_ratio

        if kmer_size < 2:
            raise ValueError("kmer_size must be at least 2.")

        self.k = kmer_size

        # graph[left][right] = coverage
        self.graph = defaultdict(lambda: defaultdict(int))
        
        self.indegree = defaultdict(int)
        self.outdegree = defaultdict(int)



    def clear(self):
        """Remove all nodes and edges."""

        self.graph.clear()
        self.indegree.clear()
        self.outdegree.clear()

    def build(self, reads: List[str]):
        """
        Construct a weighted De Bruijn graph from DNA reads.
        """

        self.clear()
        self.indegree.clear()
        self.outdegree.clear()

        for read in reads:

            if len(read) < self.k:
                continue

            for i in range(len(read) - self.k + 1):

                left = read[i:i+self.k-1]
                right = read[i+1:i+self.k]

                self.graph[left][right] += 1
                self.outdegree[left] += 1
                self.indegree[right] += 1

    def prune(self, min_coverage=1):
        """
        Remove edges with coverage below min_coverage.
        """

        for left in list(self.graph):

            for right in list(self.graph[left]):

                if self.graph[left][right] < min_coverage:
                    del self.graph[left][right]

            if len(self.graph[left]) == 0:
                del self.graph[left]

        self.indegree.clear()
        self.outdegree.clear()

        for left, neighbours in self.graph.items():

            self.outdegree[left] = len(neighbours)

            for right in neighbours:
                self.indegree[right] += 1

    def start_nodes(self):
        nodes = set(self.graph)

        for n in self.graph.values():
            nodes.update(n)

        return [n for n in nodes if self.indegree[n] == 0]

    def end_nodes(self):
        """
        Return nodes with outdegree zero.
        """

        nodes = set(self.graph.keys())

        for neighbors in self.graph.values():
            nodes.update(neighbors.keys())

        return [
            node
            for node in nodes
            if self.outdegree[node] == 0
        ]

    def best_successor(self, node):

        if node not in self.graph:
            return None

        if not self.graph[node]:
            return None

        def score(neighbour):

            weight = self.graph[node][neighbour]

            future = sum(
                self.graph.get(neighbour, {}).values()
            )

            return (weight, future, neighbour)

        return max(self.graph[node], key=score)

    def __len__(self):
        return len(self.graph)

    def __str__(self):

        edges = sum(len(v) for v in self.graph.values())

        return (
            f"DeBruijnGraph("
            f"k={self.k}, "
            f"nodes={len(self.graph)}, "
            f"edges={edges})"
        )

    def consensus_not_fixed(self):

        if len(self.graph) == 0:
            return ""

        starts = self.start_nodes()

        if starts:
            start = max(
                starts,
                key=lambda n: sum(self.graph[n].values())
            )
        else:
            start = max(
                self.graph,
                key=lambda n: sum(self.graph[n].values())
            )

        sequence = start
        current = start

        visited_edges = set()

        while True:

            nxt = self.best_successor(current)

            if nxt is None:
                break

            edge = (current, nxt)

            if edge in visited_edges:
                break

            visited_edges.add(edge)

            sequence += nxt[-1]
            current = nxt

        return sequence
    
    def consensus(self, target_length):

        if target_length < self.k - 1:
            raise ValueError(
                f"target_length must be at least {self.k - 1}"
            )

        if len(self.graph) == 0:
            return ""

        starts = self.start_nodes()

        if starts:
            start_nodes = sorted(
                starts,
                key=lambda n: sum(self.graph[n].values()),
                reverse=True
            )
        else:
            start_nodes = sorted(
                self.graph,
                key=lambda n: sum(self.graph[n].values()),
                reverse=True
            )

        def find_path(node, sequence, visited_edges):

            # Exactly the requested length
            if len(sequence) == target_length:
                return sequence

            # Too long
            if len(sequence) > target_length:
                return None

            neighbours = self.graph.get(node, {})

            # Try strongest edges first
            ordered_neighbours = sorted(
                neighbours,
                key=lambda n: self.graph[node][n],
                reverse=True
            )

            for nxt in ordered_neighbours:

                edge = (node, nxt)

                if edge in visited_edges:
                    continue

                result = find_path(
                    nxt,
                    sequence + nxt[-1],
                    visited_edges | {edge}
                )

                if result is not None:
                    return result

            return None

        # Try each possible starting node
        for start in start_nodes:

            result = find_path(
                start,
                start,
                set()
            )

            if result is not None:
                return result

        # No path of exactly target_length exists
        return ""
    
    def remove_branches(self, ratio: float = 0.2):
        """
        Remove weak competing branches.

        For every node with multiple outgoing edges, keep only edges whose
        coverage is at least `ratio` times the strongest outgoing edge.

        Example
        -------
        A -> B (100)
        A -> C (18)
        A -> D (4)

        ratio = 0.2

        Keeps:
            A -> B
        Removes:
            A -> C
            A -> D
        """

        for node in list(self.graph):

            neighbours = self.graph[node]

            if len(neighbours) <= 1:
                continue

            strongest = max(neighbours.values())

            threshold = strongest * ratio

            for nxt in list(neighbours):

                if neighbours[nxt] < threshold:
                    del neighbours[nxt]

            if len(neighbours) == 0:
                del self.graph[node]

        # Recompute degrees
        self.indegree.clear()
        self.outdegree.clear()

        for left, neighbours in self.graph.items():

            self.outdegree[left] = len(neighbours)

            for right in neighbours:
                self.indegree[right] += 1
    
    def print_graph(self):

        for node in sorted(self.graph):

            print(node)

            for nxt, weight in self.graph[node].items():

                print(
                    f"   -> {nxt} ({weight})"
                )
        
    def stats(self):

        edges = sum(len(v) for v in self.graph.values())
        coverage = sum(
            sum(v.values())
            for v in self.graph.values()
        )

        return {
            "nodes": len(self.graph),
            "edges": edges,
            "coverage": coverage,
            "starts": len(self.start_nodes()),
            "ends": len(self.end_nodes())
        }
