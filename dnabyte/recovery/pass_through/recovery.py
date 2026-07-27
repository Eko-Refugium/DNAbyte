"""
Pass-through consensus reconstruction.
"""

from collections import defaultdict
from typing import Dict, List

from dnabyte.consensus import Consensus
from dnabyte.recovery.debruijn.graph import DeBruijnGraphGraph



class PassThrough(Consensus):
    def __init__(self, params, logger=None):

        self.logger = logger
        
        
    def recover(self, clustered_data):
        readsfinal = []
        for cluster_id, reads in clustered_data.items():
            readsfinal.extend(reads)



        return readsfinal, {}
def attributes(params):

    return {
    }