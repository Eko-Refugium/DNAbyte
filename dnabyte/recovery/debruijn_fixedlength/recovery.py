"""
Weighted De Bruijn graph consensus recovery.
"""

from dnabyte.consensus import Consensus
from dnabyte.recovery.debruijn.graph import DeBruijnGraphGraph


class DeBruijnGraph(Consensus):

    def __init__(self, params, logger=None):

        self.logger = logger

        self.min_coverage = getattr(params, "min_coverage", 1)
        self.branch_ratio = getattr(params, "branch_ratio", 0.2)

        # Keep the existing parameter name used by the project.
        self.kmer_size = getattr(
            params,
            "kemer_size_debruijn",
            21
        )

        # Encoding-defined rigid sequence length.
        self.sequence_length = getattr(
            params,
            "sequence_length",
            None
        )

        self.left_primer = getattr(
            params,
            "left_primer",
            None
        )

        self.right_primer = getattr(
            params,
            "right_primer",
            None
        )

        self.trim_primers = getattr(
            params,
            "trim_primers",
            False
        )

    def recover(self, clustered_data):
        """
        Recover one consensus sequence for each cluster.

        Parameters
        ----------
        clustered_data : dict
            Dictionary containing clusters of reads.

        Returns
        -------
        consensus : list[str]
            One consensus sequence for every successfully
            reconstructed cluster.

        stats : dict
            Recovery statistics.
        """

        consensus = []

        total_reads = 0
        failed_clusters = 0

        # Process every cluster independently.
        for cluster_id, reads in clustered_data.items():

            total_reads += len(reads)

            if not reads:
                failed_clusters += 1
                continue

            # Optional primer removal / read validation.
            reads = self._prepare_reads(reads)

            if not reads:
                failed_clusters += 1
                continue

            # Build a new graph for this cluster.
            graph = DeBruijnGraphGraph(
                self.kmer_size,
                self.min_coverage,
                self.branch_ratio,
                expected_length=self.sequence_length,
                logger=self.logger
            )

            graph.build(reads)

            # Remove low-coverage edges.
            graph.prune(
                min_coverage=self.min_coverage
            )

            # Remove weak branches.
            graph.remove_branches(
                ratio=self.branch_ratio
            )

            # Recover the consensus sequence.
            sequence = graph.consensus(
                expected_length=self.sequence_length
            )

            # No valid path.
            if not sequence:
                failed_clusters += 1

                if self.logger:
                    self.logger.warning(
                        "No consensus recovered for cluster %s",
                        cluster_id
                    )

                continue

            # Rigid sequence length check.
            if (
                self.sequence_length is not None
                and len(sequence) != self.sequence_length
            ):
                failed_clusters += 1

                if self.logger:
                    self.logger.warning(
                        "Cluster %s produced sequence length %d, "
                        "expected %d",
                        cluster_id,
                        len(sequence),
                        self.sequence_length
                    )

                continue

            # One consensus sequence for this cluster.
            consensus.append(sequence)

        stats = {
            "num_clusters": len(clustered_data),
            "num_consensus": len(consensus),
            "total_reads": total_reads,
            "failed_clusters": failed_clusters,
        }

        return consensus, stats

    def _prepare_reads(self, reads):
        """
        Validate reads and optionally remove primers.
        """

        prepared = []

        for read in reads:

            if not isinstance(read, str):
                continue

            read = read.upper().strip()

            if not read:
                continue

            # Validate DNA sequence.
            if any(base not in "ACGT" for base in read):
                continue

            # Optional primer trimming.
            if self.trim_primers:

                if (
                    self.left_primer
                    and read.startswith(self.left_primer.upper())
                ):
                    read = read[
                        len(self.left_primer):
                    ]

                if (
                    self.right_primer
                    and read.endswith(self.right_primer.upper())
                ):
                    read = read[
                        :-len(self.right_primer)
                    ]

            # A read must contain at least one k-mer.
            if len(read) < self.kmer_size:
                continue

            prepared.append(read)

        return prepared


def attributes(params):
    """
    Parameters exposed by the De Bruijn recovery plugin.
    """

    return {
        "kemer_size_debruijn": getattr(
            params,
            "kemer_size_debruijn",
            21
        ),

        "min_coverage": getattr(
            params,
            "min_coverage",
            1
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

        "left_primer": getattr(
            params,
            "left_primer",
            None
        ),

        "right_primer": getattr(
            params,
            "right_primer",
            None
        ),

        "trim_primers": getattr(
            params,
            "trim_primers",
            False
        ),
    }