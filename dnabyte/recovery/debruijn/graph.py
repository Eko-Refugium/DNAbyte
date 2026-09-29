"""
Weighted De Bruijn graph for DNA consensus reconstruction.
"""

from collections import defaultdict
from typing import List


class DeBruijnGraphGraph:
    """
    Weighted De Bruijn graph.

    Nodes are (k-1)-mers.
    Edges represent k-mers and store read coverage.

    expected_length is optional. When supplied, consensus()
    requires an exact-length path.
    """

    def __init__(
        self,
        kmer_size,
        min_coverage,
        branch_ratio,
        expected_length=None,
        logger=None,
    ):
        self.logger = logger

        if kmer_size < 2:
            raise ValueError("kmer_size must be at least 2.")

        if min_coverage < 1:
            raise ValueError("min_coverage must be at least 1.")

        if not 0 < branch_ratio <= 1:
            raise ValueError("branch_ratio must be between 0 and 1.")

        if expected_length is not None and expected_length < kmer_size:
            raise ValueError(
                "expected_length must be >= kmer_size."
            )

        self.k = int(kmer_size)
        self.min_coverage = int(min_coverage)
        self.branch_ratio = float(branch_ratio)
        self.expected_length = (
            int(expected_length)
            if expected_length is not None
            else None
        )

        # graph[left][right] = coverage
        self.graph = defaultdict(lambda: defaultdict(int))

        self.indegree = defaultdict(int)
        self.outdegree = defaultdict(int)

    def clear(self):
        """Remove all graph data."""

        self.graph.clear()
        self.indegree.clear()
        self.outdegree.clear()

    def _rebuild_degrees(self):
        """Recalculate weighted node degrees."""

        self.indegree.clear()
        self.outdegree.clear()

        for left, neighbours in self.graph.items():

            for right, coverage in neighbours.items():

                self.outdegree[left] += coverage
                self.indegree[right] += coverage

    def build(self, reads: List[str]):
        """
        Build the weighted De Bruijn graph.
        """

        self.clear()

        for read in reads:

            if read is None:
                continue

            read = str(read).strip().upper()

            if len(read) < self.k:
                continue

            if any(base not in "ACGT" for base in read):
                continue

            for i in range(len(read) - self.k + 1):

                left = read[i:i + self.k - 1]
                right = read[i + 1:i + self.k]

                self.graph[left][right] += 1

        self._rebuild_degrees()

        if self.logger:
            self.logger.debug(
                f"Built De Bruijn graph: "
                f"{self.node_count()} nodes, "
                f"{self.edge_count()} edges"
            )

    def prune(self, min_coverage=None):
        """
        Remove edges below the coverage threshold.
        """

        if min_coverage is None:
            min_coverage = self.min_coverage

        min_coverage = int(min_coverage)

        for left in list(self.graph):

            neighbours = self.graph[left]

            for right in list(neighbours):

                if neighbours[right] < min_coverage:
                    del neighbours[right]

            if not neighbours:
                del self.graph[left]

        self._rebuild_degrees()

    def remove_branches(self, ratio=None):
        """
        Remove weak outgoing branches.
        """

        if ratio is None:
            ratio = self.branch_ratio

        removed = 0

        for node in list(self.graph):

            neighbours = self.graph[node]

            if len(neighbours) <= 1:
                continue

            strongest = max(neighbours.values())
            threshold = strongest * ratio

            for nxt in list(neighbours):

                if neighbours[nxt] < threshold:
                    del neighbours[nxt]
                    removed += 1

            if not neighbours:
                del self.graph[node]

        self._rebuild_degrees()

        if self.logger:
            self.logger.debug(
                f"Removed {removed} weak branch edges"
            )

    def start_nodes(self):
        """Return nodes without incoming edges."""

        nodes = set(self.graph.keys())

        for neighbours in self.graph.values():
            nodes.update(neighbours.keys())

        return [
            node
            for node in nodes
            if self.indegree[node] == 0
        ]

    def end_nodes(self):
        """Return nodes without outgoing edges."""

        nodes = set(self.graph.keys())

        for neighbours in self.graph.values():
            nodes.update(neighbours.keys())

        return [
            node
            for node in nodes
            if self.outdegree[node] == 0
        ]

    def _candidate_starts(self):
        """
        Return possible starting nodes.

        Prefer true graph starts, but allow cyclic graphs.
        """

        starts = self.start_nodes()

        if starts:
            return starts

        return list(self.graph.keys())

    def _find_exact_path(self, start, target_length):
        """
        Find a path producing exactly target_length bases.

        Search is coverage-guided but length-constrained.

        Returns:
            sequence or None
        """

        sequence = start
        current = start

        # Number of added bases still required.
        remaining = target_length - len(sequence)

        if remaining < 0:
            return None

        visited_edges = set()

        def search(node, seq, remaining):

            if remaining == 0:
                return seq

            neighbours = self.graph.get(node)

            if not neighbours:
                return None

            candidates = []

            for nxt, coverage in neighbours.items():

                edge = (node, nxt)

                if edge in visited_edges:
                    continue

                candidates.append(
                    (
                        coverage,
                        sum(
                            self.graph.get(nxt, {}).values()
                        ),
                        nxt,
                    )
                )

            # Strongest-supported paths first.
            candidates.sort(
                key=lambda x: (x[0], x[1], x[2]),
                reverse=True,
            )

            for _, _, nxt in candidates:

                # Every graph edge adds exactly one base.
                if remaining < 1:
                    continue

                edge = (node, nxt)

                visited_edges.add(edge)

                result = search(
                    nxt,
                    seq + nxt[-1],
                    remaining - 1,
                )

                if result is not None:
                    return result

                visited_edges.remove(edge)

            return None

        return search(
            current,
            sequence,
            remaining,
        )

    def _best_unconstrained_path(self):
        """
        Generate the strongest path when no exact length is required.
        """

        starts = self._candidate_starts()

        if not starts:
            return ""

        start = max(
            starts,
            key=lambda node: self.outdegree[node]
        )

        sequence = start
        current = start

        visited_edges = set()

        while True:

            neighbours = self.graph.get(current)

            if not neighbours:
                break

            candidates = [
                nxt
                for nxt in neighbours
                if (current, nxt) not in visited_edges
            ]

            if not candidates:
                break

            nxt = max(
                candidates,
                key=lambda n: (
                    self.graph[current][n],
                    sum(self.graph.get(n, {}).values()),
                    n,
                )
            )

            edge = (current, nxt)

            visited_edges.add(edge)

            sequence += nxt[-1]
            current = nxt

        return sequence

    def consensus(self, expected_length=None):
        """
        Generate a consensus sequence.

        If expected_length is supplied, the returned sequence must
        have exactly that length.

        If no exact-length path exists, return "".
        """

        if not self.graph:
            return ""

        if expected_length is None:
            expected_length = self.expected_length

        if expected_length is not None:

            expected_length = int(expected_length)

            if expected_length < self.k:
                return ""

            starts = self._candidate_starts()

            # Try the strongest starts first.
            starts = sorted(
                starts,
                key=lambda node: self.outdegree[node],
                reverse=True,
            )

            for start in starts:

                sequence = self._find_exact_path(
                    start,
                    expected_length,
                )

                if sequence is not None:

                    if len(sequence) == expected_length:
                        return sequence

            if self.logger:
                self.logger.warning(
                    "No exact-length De Bruijn path found "
                    f"for length {expected_length}"
                )

            return ""

        return self._best_unconstrained_path()

    def edge_count(self):
        """Return number of unique edges."""

        return sum(
            len(neighbours)
            for neighbours in self.graph.values()
        )

    def node_count(self):
        """Return number of graph nodes."""

        nodes = set(self.graph.keys())

        for neighbours in self.graph.values():
            nodes.update(neighbours.keys())

        return len(nodes)

    def total_coverage(self):
        """Return total edge coverage."""

        return sum(
            sum(neighbours.values())
            for neighbours in self.graph.values()
        )

    def stats(self):
        """Return graph statistics."""

        return {
            "nodes": self.node_count(),
            "edges": self.edge_count(),
            "coverage": self.total_coverage(),
            "starts": len(self.start_nodes()),
            "ends": len(self.end_nodes()),
            "expected_length": self.expected_length,
        }

    def print_graph(self):

        for node in sorted(self.graph):

            print(node)

            neighbours = sorted(
                self.graph[node].items(),
                key=lambda item: (-item[1], item[0])
            )

            for nxt, coverage in neighbours:

                print(
                    f"   -> {nxt} ({coverage})"
                )

    def __len__(self):
        return self.node_count()

    def __str__(self):

        return (
            f"DeBruijnGraph("
            f"k={self.k}, "
            f"nodes={self.node_count()}, "
            f"edges={self.edge_count()}, "
            f"expected_length={self.expected_length})"
        )
