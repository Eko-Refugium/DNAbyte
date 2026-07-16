"""
End-to-end test demonstrating the plugin architecture where params only specify method names.

This test validates that:
1. Params loads plugins automatically from method names
2. Process class uses plugins generically without hardcoding
3. All defaults come from module.attributes() functions
4. Complete pipeline works with just method names in params
"""

import pytest
from dnabyte.params import Params
from dnabyte.process import Process
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.data_classes.nucleobasecode import NucleobaseCode


class TestParamsPluginArchitecture:
    """Validate plugin architecture with params-based method specification."""
    
    def test_params_specifies_only_method_names(self):
        """
        Test that params only needs method names, not method-specific parameters.
        
        All method-specific defaults come from module.attributes() functions.
        """
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Params should only have method NAMES, not implementations
        assert params.clustering_method == 'primer_grouper'
        assert params.recovery_method == 'debruijn'
        
        # Method-specific defaults should be in params (loaded from module.attributes)
        assert hasattr(params, 'primer_length'), "primer_length should come from clustering.attributes()"
        assert params.primer_length == 20, "Default primer_length from attributes"
        
        assert hasattr(params, 'debruijn_kmer_size'), "debruijn_kmer_size should come from recovery.attributes()"
        assert params.debruijn_kmer_size == 21, "Default kmer_size from attributes"
        
        assert hasattr(params, 'debruijn_min_coverage'), "debruijn_min_coverage should come from recovery.attributes()"
        assert params.debruijn_min_coverage == 1, "Default min_coverage from attributes"
    
    def test_plugins_loaded_from_method_names(self):
        """
        Test that plugin dictionaries are populated based on method names in params.
        
        Demonstrates that params.post_sequencing_plugins contains the discovered modules.
        """
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Plugins should be discovered and stored
        assert hasattr(params, 'post_sequencing_plugins')
        assert 'clustering' in params.post_sequencing_plugins
        assert 'recovery' in params.post_sequencing_plugins
        
        # Each should have the method name as key
        clustering_plugins = params.post_sequencing_plugins['clustering']
        recovery_plugins = params.post_sequencing_plugins['recovery']
        
        assert 'primer_grouper' in clustering_plugins
        assert 'debruijn' in recovery_plugins
        
        # Values should be the actual classes
        assert clustering_plugins['primer_grouper'] is not None
        assert recovery_plugins['debruijn'] is not None
    
    def test_process_uses_plugins_generically(self):
        """
        Test that Process class works with any plugin combination specified in params.
        
        Process doesn't hardcode method names - it gets them from params.
        """
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Create Process instance
        processor = Process(params)
        
        # Verify Process loaded the methods from params
        assert processor.clustering_method == 'primer_grouper'
        assert processor.recovery_method == 'debruijn'
        
        # Verify Process has the plugin dictionaries
        assert 'primer_grouper' in processor.clustering_module
        assert 'debruijn' in processor.recovery_module
    
    def test_full_pipeline_with_params_plugin_architecture(self):
        """
        Simplified test demonstrating the plugin architecture works end-to-end.
        
        Uses mock data instead of full encoding pipeline to focus on plugin architecture:
        1. Process class loads plugins from params
        2. Plugins are called generically through cluster() and recover()
        3. All stats are collected and returned in info dict
        """
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Create Process instance with params-specified methods
        processor = Process(params)
        
        # Verify it loaded the right plugins
        assert processor.clustering_method == 'primer_grouper'
        assert processor.recovery_method == 'debruijn'
        
        # Create mock sequenced data with primers
        left_primer = "ACGTACGTACGTACGTACGT"
        right_primer = "TGCATGCATGCATGCATGCA"
        payload = "A" * 50
        
        # Generate reads
        sequences = []
        for _ in range(20):
            sequences.append(left_primer + payload + right_primer)
        
        # Add some reads with errors (realistic sequencing)
        for i in range(3):
            seq = list(left_primer + payload + right_primer)
            seq[40 + i*10] = 'T' if seq[40 + i*10] != 'T' else 'A'
            sequences.append(''.join(seq))
        
        mock_data = InSilicoDNA(sequences)
        
        # Process using plugins from params
        data_processed, info = processor.process(mock_data)
        
        # Verify processing succeeded with params-based plugins
        assert data_processed is not None
        assert info['status'] == 'success'
        assert info['clustering_method'] == 'primer_grouper'
        assert info['recovery_method'] == 'debruijn'
        
        # Verify stats from both plugins
        assert 'num_groups' in info
        assert 'num_consensus' in info or 'consensus_sequences' in info
        
        # Verify output is valid
        assert data_processed is not None
        assert isinstance(data_processed, NucleobaseCode)
        assert len(data_processed.data) > 0
    
    def test_multiple_method_combinations(self):
        """
        Test that different method combinations all work through generic interface.
        
        Currently supports:
        - clustering: primer_grouper
        - recovery: debruijn
        
        But architecture is extensible for new methods.
        """
        test_cases = [
            ('church', 'nosynthpoly', 'primer_grouper', 'debruijn'),
            ('church', 'mesa', 'primer_grouper', 'debruijn'),
            ('goldman', 'nosynthpoly', 'primer_grouper', 'debruijn'),
            ('goldman', 'mesa', 'primer_grouper', 'debruijn'),
            ('gcplus', 'nosynthpoly', 'primer_grouper', 'debruijn'),
            ('gcplus', 'mesa', 'primer_grouper', 'debruijn'),
        ]
        
        for encoding, synthesis, clustering, recovery in test_cases:
            params = Params(
                encoding_method=encoding,
                synthesis_method=synthesis,
                clustering_method=clustering,
                recovery_method=recovery
            )
            
            # Verify params loaded correctly
            assert params.clustering_method == clustering
            assert params.recovery_method == recovery
            
            # Verify plugins are available
            assert clustering in params.post_sequencing_plugins['clustering']
            assert recovery in params.post_sequencing_plugins['recovery']
            
            # Verify Process can be instantiated
            processor = Process(params)
            assert processor is not None
    
    def test_params_attributes_from_modules(self):
        """
        Test that method-specific parameters come from module.attributes() functions.
        
        This validates the clean separation where:
        - params only specifies method NAMES
        - defaults come from module.attributes()
        - Process doesn't hardcode any defaults
        """
        params = Params(
            encoding_method='church',
            synthesis_method='nosynthpoly',
            clustering_method='primer_grouper',
            recovery_method='debruijn'
        )
        
        # Get defaults from modules directly
        from dnabyte.post_sequencing.clustering.primer_grouper import attributes as clustering_attrs
        from dnabyte.post_sequencing.recovery.debruijn import attributes as recovery_attrs
        
        clustering_defaults = clustering_attrs(params)
        recovery_defaults = recovery_attrs(params)
        
        # Verify params has these attributes
        for key, value in clustering_defaults.items():
            assert hasattr(params, key)
            assert getattr(params, key) == value
        
        for key, value in recovery_defaults.items():
            assert hasattr(params, key)
            assert getattr(params, key) == value
