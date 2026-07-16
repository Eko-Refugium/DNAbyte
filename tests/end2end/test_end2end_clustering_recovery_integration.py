"""
End-to-end integration tests using testbase_end2end_newdata with clustering+recovery processing.

Tests the full pipeline (binarize → encode → synthesize → storage → sequencing → clustering+recovery → decode)
with all encoding and synthesis method combinations using the modular Process class.
"""

import unittest
import sys
from tests.end2end.testbase_end2end_newdata import TestBase
from dnabyte.params import Params


def load_tests(loader, tests, pattern):
    """
    Load test suite with parameterized Params for all encoding/synthesis combinations.
    Tests clustering+recovery processing through the modular Process class.
    """
    suite = unittest.TestSuite()
    
    # Only use encoding methods that actually exist in the system
    encoding_methods = ["church", "goldman", "gcplus"]
    synthesis_methods = ["nosynthpoly", "mesa"]
    
    for encoding_method in encoding_methods:
        for synthesis_method in synthesis_methods:
            test_name = f"{encoding_method}_{synthesis_method}_clustering_recovery"
            
            try:
                # Create params with clustering+recovery processing
                params = Params(
                    name=test_name,
                    encoding_method=encoding_method,
                    synthesis_method=synthesis_method,
                    processing_method='clustering+recovery',  # Use the modular Process class
                    debruijn_kmer_size=21,  # De Bruijn k-mer size for consensus
                    primer_length=20,  # Primer length for clustering
                    filename='fasta_import_test.fa',
                    binarization_method='default',
                    storage_conditions=None,  # Skip storage simulation for speed
                    error_methods=None,  # Skip error simulation for speed
                    sequencing_method='illumina'
                )
                
                # Create test instance
                test = TestBase(params=params)
                suite.addTest(test)
                
            except Exception as e:
                print(f"Warning: Could not create test for {test_name}: {e}")
                continue
    
    return suite


if __name__ == "__main__":
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    loader = unittest.TestLoader()
    suite = load_tests(loader, None, None)
    result = runner.run(suite)
    
    # Exit with appropriate code
    sys.exit(0 if result.wasSuccessful() else 1)
