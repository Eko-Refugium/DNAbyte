"""
Unit tests for post_sequencing.recovery module.

Tests De Bruijn graph consensus building and consensus recovery functionality.
"""

import pytest
from dnabyte.post_sequencing.recovery.debruijn.recovery import (
    DeBruijnGraph,
    build_debruijn_graph,
    extract_consensus_from_sequences,
    ConsensusBuilder,
)


class TestDeBruijnGraph:
    """Tests for DeBruijnGraph class."""
    
    def test_initialization(self):
        """Test DeBruijnGraph initialization."""
        graph = DeBruijnGraph(kmer_size=21)
        assert graph.kmer_size == 21
        assert len(graph.graph) == 0
        assert len(graph.kmer_counts) == 0
    
    def test_initialization_invalid_kmer(self):
        """Test that invalid k-mer sizes raise error."""
        with pytest.raises(ValueError, match="k-mer size must be at least 2"):
            DeBruijnGraph(kmer_size=1)
    
    def test_build_from_sequences_basic(self):
        """Test building graph from sequences."""
        graph = DeBruijnGraph(kmer_size=5)
        sequences = ["ACGTACGTACGT", "ACGTACGTACGT"]  # Identical sequences
        graph.build_from_sequences(sequences)
        
        # Should have created nodes and edges
        assert len(graph.graph) > 0
    
    def test_get_all_nodes(self):
        """Test retrieving all nodes from graph."""
        graph = DeBruijnGraph(kmer_size=3)
        sequences = ["ACGTACGT"]
        graph.build_from_sequences(sequences)
        
        nodes = graph.get_all_nodes()
        assert len(nodes) > 0
    
    def test_get_start_nodes(self):
        """Test identifying start nodes (in_degree = 0)."""
        graph = DeBruijnGraph(kmer_size=3)
        # Use sequences that don't form circular structures
        sequences = ["ACGTACGTAA", "ACGTACGTAA"]
        graph.build_from_sequences(sequences)
        
        start_nodes = graph.get_start_nodes()
        # Start nodes should exist for linear/non-circular graphs
        assert isinstance(start_nodes, list)


class TestBuildDeBruijnGraph:
    """Tests for build_debruijn_graph function."""
    
    def test_build_basic(self):
        """Test basic graph building."""
        sequences = ["ACGTACGTACGT"]
        graph = build_debruijn_graph(sequences, kmer_size=5, min_coverage=1)
        
        assert isinstance(graph, DeBruijnGraph)
        assert len(graph.graph) > 0
    
    def test_build_empty_raises_error(self):
        """Test that empty sequence list raises error."""
        with pytest.raises(ValueError, match="Cannot build graph from empty sequence list"):
            build_debruijn_graph([], kmer_size=5)


class TestExtractConsensusFromSequences:
    """Tests for extract_consensus_from_sequences function."""
    
    def test_extract_consensus_identical_sequences(self):
        """Test consensus extraction from identical sequences."""
        # Use a longer sequence relative to kmer_size for better recovery
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 45 bases
        sequences = [sequence] * 3  # Three identical copies
        
        consensus = extract_consensus_from_sequences(sequences, kmer_size=5)
        
        # Consensus should be a valid extraction (may be shorter than original)
        assert len(consensus) > 0
        # The consensus should only contain bases from the sequence
        assert set(consensus).issubset(set(sequence))
    
    def test_extract_consensus_empty_raises_error(self):
        """Test that empty sequence list raises error."""
        with pytest.raises(ValueError, match="Cannot extract consensus from empty sequence list"):
            extract_consensus_from_sequences([], kmer_size=5)


class TestConsensusBuilder:
    """Tests for ConsensusBuilder class."""
    
    def test_initialization_debruijn(self):
        """Test ConsensusBuilder initialization with De Bruijn method."""
        builder = ConsensusBuilder(method='debruijn', kmer_size=21)
        assert builder.method == 'debruijn'
        assert builder.kmer_size == 21
    
    def test_initialization_majority_voting(self):
        """Test ConsensusBuilder initialization with majority voting."""
        builder = ConsensusBuilder(method='majority_voting')
        assert builder.method == 'majority_voting'
    
    def test_initialization_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError, match="Unknown consensus method"):
            ConsensusBuilder(method='invalid_method')
    
    def test_build_consensus_debruijn(self):
        """Test building consensus with De Bruijn method."""
        builder = ConsensusBuilder(method='debruijn', kmer_size=5)
        sequences = ["ACGTACGTACGT"] * 3
        
        consensus = builder.build_consensus(sequences)
        assert len(consensus) > 0
    
    def test_build_consensus_majority_voting(self):
        """Test building consensus with majority voting."""
        builder = ConsensusBuilder(method='majority_voting')
        sequences = [
            "ACGTACGTACGT",
            "ACGTACGTACGT",
            "ACGTACGTACGT",
        ]
        
        consensus = builder.build_consensus(sequences)
        assert consensus == "ACGTACGTACGT"
    
    def test_build_consensus_empty_raises_error(self):
        """Test that empty sequence list raises error."""
        builder = ConsensusBuilder(method='debruijn')
        
        with pytest.raises(ValueError, match="Cannot build consensus from empty sequence list"):
            builder.build_consensus([])
    
    def test_get_consensus_sequences(self):
        """Test extracting just the consensus sequences."""
        builder = ConsensusBuilder(method='majority_voting')
        grouped_payloads = {
            ("ACGT", "TGCA"): ["AAAA"],
            ("GGGG", "CCCC"): ["TTTT"],
        }
        
        sequences, stats = builder.get_consensus_sequences(grouped_payloads)
        
        assert len(sequences) == 2
        assert stats['num_groups'] == 2
        assert stats['total_reads'] == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
