"""
Weighted De Bruijn graph for DNA consensus reconstruction.
"""

from collections import defaultdict
from typing import Dict, List

from dnabyte.consensus import Consensus



class Simple(Consensus):
    def __init__(self, params, logger=None):

        self.logger = logger


    def recover(self, clustered_data):
        """
        Recover consensus sequences from clustered DNA reads by majority vote per basepair.

        Input:
            ClusteredDNA

        Output:
            list of consensus sequences
            stats dictionary
        """

        consensus = []

        total_reads = 0

        for reads in clustered_data.values():

            total_reads += len(reads)

            if not reads:
                continue

            # Find the length of the shortest read
            min_length = min(len(read) for read in reads)

            # Initialize a list to hold the consensus sequence
            consensus_sequence = []

            # For each position in the reads, find the most common base
            for i in range(min_length):
                base_count = defaultdict(int)
                for read in reads:
                    base_count[read[i]] += 1
                # Find the base with the maximum count
                consensus_base = max(base_count, key=base_count.get)
                consensus_sequence.append(consensus_base)

            consensus.append(''.join(consensus_sequence))


        stats = {
            "total_reads": total_reads,
            "consensus_count": len(consensus),
        }

        return consensus, stats
    
def attributes(params):

    return {
        "kemer_size_debruijn": getattr(params, "kemer_size_debruijn", 60),
        "min_coverage": getattr(params, "min_coverage", 1),
        "branch_ratio": getattr(params, "branch_ratio", 0.2),
    }