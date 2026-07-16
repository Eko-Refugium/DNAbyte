from collections import defaultdict

from dnabyte.cluster import Cluster


class PrimerGrouper(Cluster):
    """
    Groups DNA sequences by primer regions.

    This version keeps the complete sequences.
    Primers are only used as a similarity/grouping feature.
    Nothing is removed from the reads.
    """

    def __init__(self, params=None, logger=None):
        self.params = params
        self.logger = logger

        self.primer_length_cluster = getattr(
            params,
            "primer_length_cluster",
            20
        ) if params else 20


    def cluster(self, sequences):

        if self.logger:
            self.logger.info(
                f"Primer clustering {len(sequences)} sequences"
            )

        groups = defaultdict(list)

        for seq in sequences:

            left = seq[:self.primer_length_cluster]
            right = seq[-self.primer_length_cluster:]

            key = (
                left,
                right
            )

            # keep original sequence
            groups[key].append(seq)


        # convert to ClusteredDNA-compatible format
        clustered = {
            i: seqs
            for i, seqs in enumerate(groups.values())
        }


        if self.logger:
            self.logger.info(
                f"Created {len(clustered)} clusters"
            )

        info = {}

        return clustered, info

def attributes(params):

    return {
        "primer_length_cluster": getattr(
            params,
            "primer_length_cluster",
            20
        )
    }