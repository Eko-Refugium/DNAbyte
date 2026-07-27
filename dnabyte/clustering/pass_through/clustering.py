from collections import defaultdict
from typing import Dict, List, Set, Tuple

from dnabyte.cluster import Cluster
print("LOADED PASS THROUGH CLUSTER MODULE")

class PassThrough(Cluster):

    def __init__(self, params, logger=None):
        print("LOADED PASS THROUGH CLUSTER MODULE")
        self.logger = logger


    def cluster(self, sequences):
        """
        Clusters sequences by passing them through without any modifications.
        """
        info = {}

        return {
            i: [seq] for i, seq in enumerate(sequences)
        }, info




def attributes(params):

    return {
        
    }