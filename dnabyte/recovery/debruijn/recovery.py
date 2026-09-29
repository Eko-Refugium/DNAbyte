"""
Weighted De Bruijn graph consensus reconstruction.
"""

from dnabyte.consensus import Consensus
from dnabyte.recovery.debruijn.graph import DeBruijnGraphGraph


class DeBruijnGraph(Consensus):

    def __init__(self, params, logger=None):

        self.logger = logger

        self.min_coverage = getattr(
            params,
            "min_coverage",
            1
        )

        self.branch_ratio = getattr(
            params,
            "branch_ratio",
            0.2
        )

        self.kmer_size = getattr(
            params,
            "kemer_size_debruijn",
            21
        )

        self.sequence_length = getattr(
            params,
            "sequence_length",
            None
        )

    def recover(self, clustered_data):

        consensus_list = []

        for cluster_id, reads in clustered_data.items():

            print(
                f"Processing cluster {cluster_id} "
                f"with {len(reads)} reads"
            )

            if not reads:
                print(f"Cluster {cluster_id} is empty")
                continue

            graph = DeBruijnGraphGraph(
                self.kmer_size,
                self.min_coverage,
                self.branch_ratio,
                expected_length=self.sequence_length,
                logger=self.logger
            )

            # Build a completely new graph for this cluster.
            graph.build(reads)

            print(
                f"Cluster {cluster_id}: "
                f"reads={len(reads)}, "
                f"graph_nodes={len(graph)}, "
                f"expected_length={self.sequence_length}"
            )

            print("Graph stats:", graph.stats())

            graph.prune(
                min_coverage=self.min_coverage
            )

            graph.remove_branches(
                ratio=self.branch_ratio
            )

            sequence = graph.consensus(
                expected_length=self.sequence_length
            )

            print(
                f"Cluster {cluster_id}: "
                f"consensus length = {len(sequence)}"
            )

            if not sequence:
                print(
                    f"Cluster {cluster_id}: "
                    f"NO CONSENSUS"
                )
                continue

            if (
                self.sequence_length is not None
                and len(sequence) != self.sequence_length
            ):
                print(
                    f"Cluster {cluster_id}: "
                    f"wrong length {len(sequence)}, "
                    f"expected {self.sequence_length}"
                )
                continue

            # IMPORTANT:
            # append, don't overwrite.
            consensus_list.append(sequence)

            print(
                f"Cluster {cluster_id}: "
                f"CONSENSUS ADDED"
            )

            print(
                f"Number of consensus sequences so far: "
                f"{len(consensus_list)}"
            )

        print(
            f"FINAL NUMBER OF CONSENSUS SEQUENCES: "
            f"{len(consensus_list)}"
        )

        return consensus_list, {
            "num_clusters": len(clustered_data),
            "num_consensus": len(consensus_list),
            }



def attributes(params):

    return {
        "kemer_size_debruijn": getattr(
            params,
            "kemer_size_debruijn",
            60
        ),

        "min_coverage": getattr(
            params,
            "min_coverage",
            5
        ),

        "branch_ratio": getattr(
            params,
            "branch_ratio",
            0.2
        ),

        "sequence_length": getattr(
            params,
            "sequence_length",
            None
        ),
    }