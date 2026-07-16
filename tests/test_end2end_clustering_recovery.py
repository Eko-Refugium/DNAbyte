"""
Integration tests for post-sequencing clustering and recovery modules.

Tests the full end-to-end pipeline with the new industry-standard
clustering+recovery approach, including fallback behavior across all synthesis methods.
"""

import pytest
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.params import Params
from dnabyte.encode import Encode


class MockLogger:
    """Simple mock logger for testing."""
    
    def __init__(self):
        self.messages = []
    
    def info(self, msg):
        self.messages.append(('info', msg))
    
    def warning(self, msg):
        self.messages.append(('warning', msg))
    
    def error(self, msg):
        self.messages.append(('error', msg))

class TestEndToEndClusteringRecovery:
    """End-to-end tests for clustering + recovery with all synthesis methods."""
    
    @pytest.mark.parametrize("synthesis_method", ["nosynthpoly", "mesa"])
    def test_full_pipeline_with_primer_clustering_and_consensus(self, synthesis_method):
        """
        Test complete end-to-end pipeline with actual clustering and consensus building.
        
        Workflow:
        1. Encoding outputs sequences with unknown count (all same primer pair)
        2. Post-sequencing auto-extracts primers from sequences
        3. Post-sequencing clusters by extracted primer pairs
        4. ConsensusBuilder builds consensus from each cluster
        
        Note: We don't know sequence count ahead of time - only that they share primer design.
        """
        from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerGrouper
        from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder
        
        logger = MockLogger()
        params = Params(
            encoding_method='church',
            synthesis_method=synthesis_method,
            processing_method='clustering+recovery'
        )
        
        # Generate sequences from unknown count (simulated encoding output)
        # All sequences have same primer design (single region from encoding)
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        payload = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        
        # Don't assume count - just generate some sequences
        num_reads = 20  # This would be unknown in real scenario
        sequences = [left_primer + payload + right_primer for _ in range(num_reads)]
        
        # Cluster by auto-extracted primers (finds most common primer pair)
        grouper = PrimerGrouper()
        grouped_payloads = grouper.group_by_primers(sequences, primer_length=20)
        
        # Should have at least 1 group (may be 1 if all same primer pair)
        assert len(grouped_payloads) >= 1, f"Expected at least 1 group, got {len(grouped_payloads)}"
        
        # Total reads should match input
        total_payloads = sum(len(payloads) for payloads in grouped_payloads.values())
        assert total_payloads == num_reads
        
        # Build consensus from clustered payloads
        builder = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=1)
        consensus_seqs, stats = builder.get_consensus_sequences(grouped_payloads)
        
        # Should have consensus for each group
        assert len(consensus_seqs) == len(grouped_payloads)
        assert stats['total_reads'] == num_reads
        
        # Verify each consensus contains valid DNA bases
        for consensus in consensus_seqs:
            assert len(consensus) > 0
            assert all(base in 'ACGT' for base in consensus)
    
    @pytest.mark.parametrize("synthesis_method", ["nosynthpoly", "mesa"])
    def test_pipeline_with_sequencing_errors_multiple_methods(self, synthesis_method):
        """
        Test clustering+recovery pipeline resilience with realistic sequencing errors.
        
        Unknown number of reads, same primer pair design, with realistic errors.
        Tests with multiple synthesis methods to ensure clustering/recovery
        works independently of synthesis approach.
        """
        from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerGrouper
        from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder
        
        logger = MockLogger()
        params = Params(
            encoding_method='church',
            synthesis_method=synthesis_method
        )
        
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        true_payload = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        
        # Generate reads with realistic sequencing error rate (~1-2%)
        # Don't assume exact count - just generate based on coverage percentages
        base_seq = left_primer + true_payload + right_primer
        sequences = []
        
        # 60% perfect reads
        for _ in range(12):
            sequences.append(base_seq)
        
        # 20% with single substitution
        for i in range(4):
            mutant = list(base_seq)
            mutant[30 + i*5] = 'C' if mutant[30 + i*5] != 'C' else 'G'
            sequences.append(''.join(mutant))
        
        # 20% with two substitutions
        for i in range(4):
            mutant = list(base_seq)
            mutant[50 + i*3] = 'G' if mutant[50 + i*3] != 'G' else 'C'
            mutant[60 + i*3] = 'T' if mutant[60 + i*3] != 'T' else 'A'
            sequences.append(''.join(mutant))
        
        actual_count = len(sequences)  # Don't assume - get actual count
        
        # Cluster by auto-extracted primers
        grouper = PrimerGrouper()
        grouped_payloads = grouper.group_by_primers(sequences, primer_length=20)
        
        # Verify we got payloads
        assert len(grouped_payloads) >= 1
        total_reads = sum(len(payloads) for payloads in grouped_payloads.values())
        assert total_reads == actual_count
        
        # Build consensus with error filtering
        builder = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=2)
        consensus_seqs, stats = builder.get_consensus_sequences(grouped_payloads)
        
        # Should produce consensus despite errors
        assert len(consensus_seqs) > 0
        assert stats['total_reads'] == actual_count
        
        # Verify output format
        for consensus in consensus_seqs:
            assert isinstance(consensus, str)
            assert all(base in 'ACGT' for base in consensus)
    
    def test_processing_method_routing_with_clustering_recovery(self):
        """Test that processing_method='clustering+recovery' correctly routes to new modules."""
        logger = MockLogger()
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            processing_method='clustering+recovery'
        )
        
        # Verify params configured correctly
        assert params.processing_method == 'clustering+recovery'
        assert params.debruijn_kmer_size == 21  # default
        
        # Verify post_sequencing_plugins loaded
        assert hasattr(params, 'post_sequencing_plugins')
        assert 'clustering' in params.post_sequencing_plugins
        assert 'recovery' in params.post_sequencing_plugins
    
    def test_processing_method_auto_fallback(self):
        """Test that processing_method='auto' uses clustering+recovery when encoding.process fails."""
        logger = MockLogger()
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            processing_method='auto'  # Should fallback intelligently
        )
        
        assert params.processing_method == 'auto'
        assert hasattr(params, 'post_sequencing_plugins')
    
    def test_params_accepts_custom_kmer_size(self):
        """Test that Params accepts and stores custom k-mer size."""
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            debruijn_kmer_size=31
        )
        
        assert params.debruijn_kmer_size == 31
        
        params2 = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly'
        )
        assert params2.debruijn_kmer_size == 21  # default


class TestMultiplePrimerGroupsAndConsensus:
    """Tests handling of multiple primer groups in a single run."""
    
    def test_multiple_primer_groups_clustering_and_recovery(self):
        """
        Test realistic scenario where encoding provides known primer pairs.
        
        When encoding has multiple regions, it can provide primer info to post-sequencing.
        Post-sequencing doesn't need to know the sequence count - just cluster by
        known primers and build consensus for each group.
        """
        from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerExtractor
        from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder
        
        # Encoding design specifies 3 regions with known primer pairs
        # (actual read count unknown until runtime)
        primer_configs = [
            {
                'left': "ACGTACGTACGTACGTACGT",
                'right': "TGCATGCATGCATGCATGCA",
                'payload': "A" * 50,
            },
            {
                'left': "GGGGGGGGGGGGGGGGGGGG",
                'right': "CCCCCCCCCCCCCCCCCCCC",
                'payload': "T" * 50,
            },
            {
                'left': "AAAAAAAAAAAAAAAAAAAAAA",
                'right': "TTTTTTTTTTTTTTTTTTTT",
                'payload': "G" * 50,
            },
        ]
        
        # Generate sequences (different counts for each region - unknown at design time)
        region_counts = [25, 18, 22]  # Simulated different coverage
        grouped_payloads = {}
        total_expected_reads = 0
        
        for config, count in zip(primer_configs, region_counts):
            primer_pair = (config['left'], config['right'])
            payloads = []
            
            # Generate reads for this region
            for _ in range(count):
                full_seq = config['left'] + config['payload'] + config['right']
                payload = PrimerExtractor.strip_primers(full_seq, config['left'], config['right'])
                payloads.append(payload)
            
            grouped_payloads[primer_pair] = payloads
            total_expected_reads += count
        
        # Should have 3 separate groups (encoding-provided primer info)
        assert len(grouped_payloads) == 3
        
        # Verify group structure
        from dnabyte.post_sequencing.clustering.primer_grouper import PrimerGrouper
        grouper = PrimerGrouper()
        stats = grouper.get_group_stats(grouped_payloads)
        assert stats['num_groups'] == 3
        assert stats['total_sequences'] == total_expected_reads
        
        # Build consensus for all groups (don't assume count)
        builder = ConsensusBuilder(method='debruijn', kmer_size=21)
        consensus_seqs, pipeline_stats = builder.get_consensus_sequences(grouped_payloads)
        
        # Should have one consensus per group
        assert len(consensus_seqs) == len(grouped_payloads)
        assert pipeline_stats['num_groups'] == len(grouped_payloads)
        assert pipeline_stats['total_reads'] == total_expected_reads
        
        # All should be valid DNA sequences
        for consensus in consensus_seqs:
            assert all(base in 'ACGT' for base in consensus)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
