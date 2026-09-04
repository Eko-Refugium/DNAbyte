from collections import defaultdict
from typing import Dict, List

from dnabyte.cluster import Cluster

class IdentityClustering(Cluster):
    """
    Simple clustering that groups identical (or near-identical) sequences together.
    Perfect for GC+ which produces fixed-length identical copies with no overlap.
    
    Groups sequences by exact match first, then optionally by edit distance
    if hamming_threshold is set.
    """

    def __init__(self, params, logger=None):
        self.hamming_threshold = getattr(params, "hamming_threshold", 0)
        self.logger = logger

    def cluster(self, sequences):
        """
        Cluster sequences by identity.
        
        Groups identical sequences together.
        If hamming_threshold > 0, also groups sequences within that edit distance.
        """
        groups_dict = defaultdict(list)
        
        # First pass: group by exact match
        for seq in sequences:
            groups_dict[seq].append(seq)
        
        # If no hamming threshold, we're done
        if self.hamming_threshold == 0:
            if self.logger:
                self.logger.info(f"Identity clustering: {len(groups_dict)} unique sequences")
            return {i: seqs for i, seqs in enumerate(groups_dict.values())}, {}
        
        # Second pass: merge groups within hamming distance
        final_groups = []
        processed = set()
        
        for representative, copies in groups_dict.items():
            if representative in processed:
                continue
            
            group = copies.copy()
            processed.add(representative)
            
            # Find other groups within hamming distance
            for other_rep, other_copies in groups_dict.items():
                if other_rep in processed:
                    continue
                if self.hamming_distance(representative, other_rep) <= self.hamming_threshold:
                    group.extend(other_copies)
                    processed.add(other_rep)
            
            final_groups.append(group)
        
        if self.logger:
            self.logger.info(
                f"Identity clustering: {len(final_groups)} groups "
                f"from {len(groups_dict)} unique sequences (hamming_threshold={self.hamming_threshold})"
            )
        
        return {i: g for i, g in enumerate(final_groups)}, {}

    @staticmethod
    def hamming_distance(seq1, seq2):
        """Calculate hamming distance between two sequences of equal length."""
        if len(seq1) != len(seq2):
            return float('inf')
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def attributes(params):
    return {
        "hamming_threshold": getattr(params, "hamming_threshold", 0),
    }
