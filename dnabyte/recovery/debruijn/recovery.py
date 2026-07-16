"""
Weighted De Bruijn graph for DNA consensus reconstruction.
"""

from collections import defaultdict
from typing import Dict, List

from dnabyte.consensus import Consensus
from dnabyte.recovery.debruijn.graph import DeBruijnGraphGraph



class DeBruijnGraph(Consensus):
    def __init__(self, params, logger=None):

        self.logger = logger

        print(params)
        
        
        self.min_coverage = getattr(params, "min_coverage")
        self.branch_ratio = getattr(params, "branch_ratio")
        
        self.kmer_size = getattr(params, "kemer_size_debruijn")

    def recover(self, clustered_data):
        """
        Recover consensus sequences from clustered DNA reads.

        Input:
            ClusteredDNA

        Output:
            list of consensus sequences
            stats dictionary
        """

        consensus = []

        total_reads = 0

        print(f"Recovering consensus sequences from {len(clustered_data)} clusters.")
        print(clustered_data)

        for reads in clustered_data.values():

            total_reads += len(reads)

            graph = DeBruijnGraphGraph(self.kmer_size, self.min_coverage, self.branch_ratio, logger=self.logger)

            graph.build(reads)

            graph.prune()

            # graph.remove_branches()

            sequence = graph.consensus()

            if sequence:
                consensus.append(sequence)


        stats = {
            "num_clusters": len(clustered_data),
            "num_consensus": len(consensus),
            "total_reads": total_reads,
        }

        print(graph.stats())

        return consensus, stats
    
def attributes(params):

    return {
        "kemer_size_debruijn": getattr(params, "kemer_size_debruijn", 85),
        "min_coverage": getattr(params, "min_coverage", 1),
        "branch_ratio": getattr(params, "branch_ratio", 0.2),
    }