from collections import defaultdict
from typing import Dict, List, Set, Tuple

from dnabyte.cluster import Cluster
class PassThrough(Cluster):

    def __init__(self, params, logger=None):
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