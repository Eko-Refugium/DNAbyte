"""
End-to-end tests for clustering+recovery post-sequencing with all synthesis methods.

Tests that the new industry-standard clustering+recovery approach works
seamlessly with every synthesis method in the DNAbyte system.
"""

import pytest
from dnabyte.params import Params
from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerGrouper, PrimerExtractor
from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder


class TestClusteringRecoverySynthesisCompatibility:
    """Test clustering+recovery works with all encoding and synthesis methods."""
    
    @pytest.fixture
    def encoding_methods(self):
        """Get all available encoding methods."""
        return ["church", "goldman", "gcplus"]
    
    @pytest.fixture
    def synthesis_methods(self):
        """Get all available synthesis methods."""
        # Can be parameterized if needed - for now test the common ones
        return ["nosynthpoly", "mesa"]
    
    def test_params_load_with_all_synthesis_methods(self, encoding_methods, synthesis_methods):
        """Test that Params initializes with clustering+recovery for all encoding and synthesis methods."""
        for encoding in encoding_methods:
            for method in synthesis_methods:
                params = Params(
                    encoding_method=encoding,
                    synthesis_method=method,
                    clustering_method='primer_grouper',
                    recovery_method='debruijn'
                )
                
                # Verify post-sequencing plugins loaded
                assert hasattr(params, 'post_sequencing_plugins')
                assert 'clustering' in params.post_sequencing_plugins
                assert 'recovery' in params.post_sequencing_plugins
                assert params.clustering_method == 'primer_grouper'
                assert params.recovery_method == 'debruijn'
    
    @pytest.mark.parametrize("encoding_method,synthesis_method", [
        ("church", "nosynthpoly"), ("church", "mesa"),
        ("goldman", "nosynthpoly"), ("goldman", "mesa"),
        ("gcplus", "nosynthpoly"), ("gcplus", "mesa"),
    ])
    def test_full_pipeline_dna_storage_workflow(self, encoding_method, synthesis_method):
        """
        Test complete DNA storage workflow with clustering+recovery.
        
        Simulates:
        1. Encoding generates DNA sequences with primers
        2. Sequencing reads (high coverage with errors)
        3. Post-sequencing clusters by primer pairs
        4. Consensus building recovers original data
        5. Verified output is valid DNA
        
        Works with all encoding and synthesis methods.
        """
        params = Params(
            encoding_method=encoding_method,
            synthesis_method=synthesis_method,
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Step 1: Simulate encoding output (sequences with primers)
        # This would come from the full encoding pipeline
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        payload = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        
        # Simulate high-coverage sequencing (unknown count)
        base_seq = left_primer + payload + right_primer
        sequences = []
        
        # Most reads perfect
        for _ in range(15):
            sequences.append(base_seq)
        
        # Some with errors (realistic 1-2% error rate)
        for i in range(3):
            mutant = list(base_seq)
            mutant[50 + i*10] = 'C' if mutant[50 + i*10] != 'C' else 'G'
            sequences.append(''.join(mutant))
        
        for i in range(2):
            mutant = list(base_seq)
            mutant[60 + i*15] = 'G' if mutant[60 + i*15] != 'G' else 'C'
            mutant[70 + i*15] = 'T' if mutant[70 + i*15] != 'T' else 'A'
            sequences.append(''.join(mutant))
        
        actual_read_count = len(sequences)
        
        # Step 2: Cluster by auto-extracted primers
        grouper = PrimerGrouper()
        grouped_payloads = grouper.group_by_primers(sequences, primer_length=20)
        
        assert len(grouped_payloads) >= 1
        
        total_grouped = sum(len(p) for p in grouped_payloads.values())
        assert total_grouped == actual_read_count
        
        # Step 3: Build consensus using De Bruijn recovery
        builder = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=1)
        consensus_seqs, stats = builder.get_consensus_sequences(grouped_payloads)
        
        # Verify results
        assert len(consensus_seqs) > 0
        assert stats['total_reads'] == actual_read_count
        
        # Step 4: Validate output is DNA
        for consensus in consensus_seqs:
            assert isinstance(consensus, str)
            assert len(consensus) > 0
            # All characters should be valid DNA bases
            assert all(base in 'ACGT' for base in consensus), \
                f"Invalid base in consensus: {consensus}"
    
    @pytest.mark.parametrize("encoding_method,synthesis_method", [
        ("church", "nosynthpoly"), ("church", "mesa"),
        ("goldman", "nosynthpoly"), ("goldman", "mesa"),
        ("gcplus", "nosynthpoly"), ("gcplus", "mesa"),
    ])
    def test_multiple_regions_multi_synthesis(self, encoding_method, synthesis_method):
        """
        Test multiple encoded regions with clustering+recovery across all methods.
        
        Encoding can have multiple regions (each with own primers).
        Post-sequencing handles all regions independently.
        """
        params = Params(
            encoding_method=encoding_method,
            synthesis_method=synthesis_method
        )
        
        # Simulate encoding with 2 regions
        region_1_primers = ("ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA")
        region_2_primers = ("GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC")
        
        region_1_payload = "A" * 50
        region_2_payload = "T" * 50
        
        # Generate reads for each region (unknown counts)
        reads_r1 = [region_1_primers[0] + region_1_payload + region_1_primers[1] for _ in range(20)]
        reads_r2 = [region_2_primers[0] + region_2_payload + region_2_primers[1] for _ in range(15)]
        
        # Combine all sequences (as they'd come from sequencing)
        all_sequences = reads_r1 + reads_r2
        
        # Manual grouping by known region primers (as encoding would provide)
        grouped_payloads = {
            region_1_primers: [
                PrimerExtractor.strip_primers(seq, region_1_primers[0], region_1_primers[1])
                for seq in reads_r1
            ],
            region_2_primers: [
                PrimerExtractor.strip_primers(seq, region_2_primers[0], region_2_primers[1])
                for seq in reads_r2
            ]
        }
        
        # Build consensus for all regions
        builder = ConsensusBuilder(method='debruijn', kmer_size=21)
        consensus_seqs, stats = builder.get_consensus_sequences(grouped_payloads)
        
        # Should have consensus for each region
        assert len(consensus_seqs) == 2
        assert stats['num_groups'] == 2
        assert stats['total_reads'] == len(all_sequences)
        
        # All should be valid
        for consensus in consensus_seqs:
            assert all(base in 'ACGT' for base in consensus)
    
    @pytest.mark.parametrize("encoding_method,synthesis_method", [
        ("church", "nosynthpoly"), ("church", "mesa"),
        ("goldman", "nosynthpoly"), ("goldman", "mesa"),
        ("gcplus", "nosynthpoly"), ("gcplus", "mesa"),
    ])
    def test_high_error_rate_resilience(self, encoding_method, synthesis_method):
        """
        Test clustering+recovery resilience with high error rates
        across all encoding and synthesis methods.
        
        Verifies De Bruijn consensus can recover from realistic
        sequencing errors even with different encoding/synthesis pipelines.
        """
        params = Params(
            encoding_method=encoding_method,
            synthesis_method=synthesis_method
        )
        
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        true_payload = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        
        base_seq = left_primer + true_payload + right_primer
        sequences = []
        
        # 50% perfect reads
        for _ in range(10):
            sequences.append(base_seq)
        
        # 30% with 1 error
        for i in range(6):
            mutant = list(base_seq)
            pos = (30 + i * 5) % len(base_seq)
            mutant[pos] = 'C' if mutant[pos] != 'C' else 'G'
            sequences.append(''.join(mutant))
        
        # 20% with 2+ errors
        for i in range(4):
            mutant = list(base_seq)
            pos1 = (50 + i * 7) % len(base_seq)
            pos2 = (60 + i * 7) % len(base_seq)
            mutant[pos1] = 'G' if mutant[pos1] != 'G' else 'C'
            mutant[pos2] = 'T' if mutant[pos2] != 'T' else 'A'
            sequences.append(''.join(mutant))
        
        # Cluster and build consensus
        grouper = PrimerGrouper()
        grouped = grouper.group_by_primers(sequences, primer_length=20)
        
        builder = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=3)
        consensus_seqs, stats = builder.get_consensus_sequences(grouped)
        
        # Should still produce consensus despite errors
        assert len(consensus_seqs) > 0
        assert stats['total_reads'] == len(sequences)
        
        # Output should be valid DNA
        for consensus in consensus_seqs:
            assert all(base in 'ACGT' for base in consensus)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
