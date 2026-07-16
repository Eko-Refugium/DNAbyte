"""
Unit tests for post_sequencing.clustering module.

Tests primer extraction and grouping functionality.
"""

import pytest
from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerExtractor, PrimerGrouper


class TestPrimerExtractor:
    """Tests for PrimerExtractor class."""
    
    def test_extract_by_length_basic(self):
        """Test basic primer extraction by fixed length."""
        sequences = [
            "ACGTACGTACGTACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNNTGCATGCATGCATGCATGCA",
            "ACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGGGGGGGGTGCATGCATGCATGCATGCA",
            "ACGTACGTACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAATGCATGCATGCATGCATGCA",
        ]
        primer_length = 20
        left, right = PrimerExtractor.extract_by_length(sequences, primer_length)
        
        assert left == "ACGTACGTACGTACGTACGT"
        assert right == "TGCATGCATGCATGCATGCA"
        assert len(left) == primer_length
        assert len(right) == primer_length
    
    def test_extract_by_length_empty_raises_error(self):
        """Test that empty sequence list raises ValueError."""
        with pytest.raises(ValueError, match="Cannot extract primers from empty"):
            PrimerExtractor.extract_by_length([], 20)
    
    def test_strip_primers_basic(self):
        """Test primer stripping from a sequence."""
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        sequence = left_primer + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" + right_primer
        
        payload = PrimerExtractor.strip_primers(sequence, left_primer, right_primer)
        assert payload == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"


class TestPrimerGrouper:
    """Tests for PrimerGrouper class."""
    
    def test_group_by_primers_basic(self):
        """Test basic grouping by primer pairs with single-copy payloads."""
        sequences = [
            "ACGTACGTACGTACGTACGTAAAAAAA1ATGCATGCATGCATGCATGCA",
            "ACGTACGTACGTACGTACGTAAAAAAA2ATGCATGCATGCATGCATGCA",
            "ACGTACGTACGTACGTACGTAAAAAAA3ATGCATGCATGCATGCATGCA",
        ]
        
        primer_length = 20
        groups = PrimerGrouper.group_by_primers(sequences, primer_length)
        
        # With different payloads (single-copy mode), each unique payload gets its own group
        # Key format: (left_primer, right_primer, payload_hash)
        assert len(groups) == 3  # 3 unique payloads -> 3 groups
        
        # Verify each group has expected structure
        for key, payloads in groups.items():
            # Each key should be 3-element tuple (left, right, payload_hash)
            assert len(key) == 3
            assert key[0] == "ACGTACGTACGTACGTACGT"
            assert key[1] == "TGCATGCATGCATGCATGCA"
            # Each group should have exactly one payload for single-copy
            assert len(payloads) == 1
    
    def test_group_by_primers_with_coverage(self):
        """Test grouping by primer pairs with identical payloads (coverage case)."""
        sequences = [
            "ACGTACGTACGTACGTACGTAAAAAAAAATGCATGCATGCATGCATGCA",  # 10 copies of same payload
            "ACGTACGTACGTACGTACGTAAAAAAAAATGCATGCATGCATGCATGCA",
            "ACGTACGTACGTACGTACGTAAAAAAAAATGCATGCATGCATGCATGCA",
        ]
        
        primer_length = 20
        groups = PrimerGrouper.group_by_primers(sequences, primer_length)
        
        # With identical payloads (coverage mode), all grouped together with single key
        assert len(groups) == 1  # 1 unique payload -> 1 group
        
        key = list(groups.keys())[0]
        assert len(key) == 3  # (left_primer, right_primer, payload_hash)
        assert len(groups[key]) == 3  # 3 sequences with same payload
    
    def test_group_by_primers_empty_raises_error(self):
        """Test that empty sequence list raises ValueError."""
        with pytest.raises(ValueError, match="Cannot group empty sequence list"):
            PrimerGrouper.group_by_primers([], 20)
    
    def test_group_stats_basic(self):
        """Test statistics calculation for grouped sequences."""
        groups = {
            ("ACGT", "TGCA"): ["AAA", "AAA", "AAA"],
            ("AAAA", "TTTT"): ["GGG", "GGG"],
        }
        
        stats = PrimerGrouper.get_group_stats(groups)
        
        assert stats['num_groups'] == 2
        assert stats['total_sequences'] == 5
        assert stats['max_group_size'] == 3
        assert stats['min_group_size'] == 2
        assert stats['mean_group_size'] == 2.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
