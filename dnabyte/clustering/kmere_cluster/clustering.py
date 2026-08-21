from collections import defaultdict
from typing import Dict, List, Set, Tuple

from dnabyte.cluster import Cluster
print("LOADED KMERE CLUSTER MODULE")

class KmerClusterer(Cluster):

    def __init__(self, params, logger=None):
        print("LOADED KMERE CLUSTER MODULE")
        self.k = getattr(params, "kmer_size_cluster",15)
        self.threshold = getattr(params, "kmer_threshold", 0.7)
        self.logger = logger


    def cluster(self, sequences):

        groups = []

        for seq in sequences:

            kmers = self.get_kmers(seq)

            found = False

            for group in groups:

                score = self.jaccard(
                    kmers,
                    group["kmers"]
                )

                if score >= self.threshold:
                    group["seqs"].append(seq)
                    found = True
                    break

            if not found:
                groups.append({
                    "kmers": kmers,
                    "seqs": [seq]
                })

        info = {}

        return {
            i: g["seqs"]
            for i, g in enumerate(groups)
        }, info


    def get_kmers(self, seq):

        return {
            seq[i:i+self.k]
            for i in range(len(seq)-self.k+1)
        }


    def jaccard(self, a, b):

        return len(a & b) / len(a | b)



def attributes(params):

    return {
        "kmer_size_cluster": getattr(params, "kmer_size_cluster", 15),
        "kmer_threshold": getattr(params, "kmer_threshold", 0.7),
    }