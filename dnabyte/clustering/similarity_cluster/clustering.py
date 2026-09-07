from collections import defaultdict
from typing import Dict, List
from dnabyte.cluster import Cluster


class SimilarityClusterer(Cluster):
    """
    Clusters DNA sequences based on sequence similarity.
    
    Groups sequences that are X% identical together, regardless of encoding method.
    Works by comparing each sequence against cluster representatives using
    sequence identity/similarity metrics.
    
    This method is:
    - Encoding-agnostic (works for church, wukong, goldman, etc.)
    - Error-tolerant (groups similar sequences even with errors)
    - Scalable (uses greedy assignment rather than pairwise comparison)
    
    Parameters:
    - similarity_threshold: Minimum percent identity to group sequences (0.0-1.0)
                          e.g., 0.85 = 85% identical
    - gap_penalty: How much to penalize gaps in alignment (default: 1.0)
    """

    def __init__(self, params=None, logger=None):
        self.similarity_threshold = getattr(params, "similarity_threshold", 0.85)
        self.gap_penalty = getattr(params, "gap_penalty", 1.0)
        self.mean = getattr(params, "mean", 1)  # Expected number of copies per sequence
        self.logger = logger

    def calculate_identity(self, seq1, seq2):
        """
        Calculate sequence identity percentage (0.0 to 1.0).
        
        Handles sequences of different lengths by allowing gaps.
        Returns the identity of the shorter sequence against the longer one.
        """
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0
        
        # Use the shorter sequence length as reference
        short_len = min(len(seq1), len(seq2))
        long_len = max(len(seq1), len(seq2))
        
        # Align to the start (simple approach)
        matches = sum(1 for i in range(short_len) if seq1[i] == seq2[i])
        
        # Penalize length difference
        identity = matches / long_len
        
        return identity

    def cluster(self, sequences):
        """
        Cluster sequences using greedy assignment based on similarity.
        
        Each sequence is assigned to the first cluster whose representative
        sequence has sufficient identity. New clusters are created as needed.
        """
        if not sequences or len(sequences) == 0:
            return {}, {}

        clusters = []
        used = [False] * len(sequences)

        for i, seq in enumerate(sequences):
            if used[i]:
                continue

            # Find or create cluster for this sequence
            found_cluster = False
            
            for cluster in clusters:
                rep_seq = cluster['representative']
                identity = self.calculate_identity(seq, rep_seq)
                
                if identity >= self.similarity_threshold:
                    cluster['sequences'].append(seq)
                    found_cluster = True
                    break
            
            # Create new cluster if no match found
            if not found_cluster:
                clusters.append({
                    'representative': seq,
                    'sequences': [seq]
                })
            
            used[i] = True

        # Convert to output format
        result = {i: cluster['sequences'] for i, cluster in enumerate(clusters)}

        # Filter out clusters with fewer than mean // 2 sequences (noise/artifacts)
        min_cluster_size = self.mean // 2
        if min_cluster_size > 0:
            filtered_result = {}
            discarded_clusters = 0
            discarded_sequences = 0
            for idx, sequences in result.items():
                if len(sequences) >= min_cluster_size:
                    filtered_result[idx] = sequences
                else:
                    discarded_clusters += 1
                    discarded_sequences += len(sequences)
            result = filtered_result
        else:
            discarded_clusters = 0
            discarded_sequences = 0

        info = {
            'similarity_threshold': self.similarity_threshold,
            'num_clusters': len(result),
            'num_sequences': sum(len(seqs) for seqs in result.values()) if result else 0,
            'avg_cluster_size': sum(len(seqs) for seqs in result.values()) / len(result) if result else 0,
            'min_cluster_size': self.mean // 2,
            'discarded_clusters': discarded_clusters,
            'discarded_sequences': discarded_sequences,
        }

        if self.logger:
            self.logger.info(
                f"Similarity clustering: {len(sequences)} sequences -> {len(result)} clusters "
                f"(avg size: {info['avg_cluster_size']:.1f}, discarded: {discarded_clusters} clusters, {discarded_sequences} sequences)"
            )

        return result, info


def attributes(params):
    return {
        "similarity_threshold": getattr(params, "similarity_threshold", 0.85),
        "gap_penalty": getattr(params, "gap_penalty", 1.0),
        "mean": getattr(params, "mean", 1),
    }
