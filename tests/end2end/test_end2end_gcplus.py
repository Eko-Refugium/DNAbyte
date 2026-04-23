import unittest
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets for GC+ encoding.
# GC+ uses an inner systematic RS code + codebook parity for per-oligo
# edit-error correction.  No primers, no file-level RS outer code in
# these tests — only the inner GC+ code is exercised.
params_list = [
    # Test 1: Basic — no errors, just synthesis copies
    Params(
        name='end2end_gcplus_basic',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='gcplus',
        binarization_method='default',
        sequence_length=200,
        gcplus_k=168,
        gcplus_l=8,
        gcplus_c1=2,

        # error channels
        storage_conditions=None,
        synthesis_method='nosynthpoly',
    ),

    # Test 2: Sequencing errors only
    Params(
        name='end2end_gcplus_seq_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='gcplus',
        binarization_method='default',
        sequence_length=200,
        gcplus_k=168,
        gcplus_l=8,
        gcplus_c1=2,

        # synthesis — create copies without errors
        synthesis_method='nosynthpoly',
        mean=10,
        std_dev=0,

        # storage — none
        storage_conditions=None,

        # sequencing errors
        sequencing_method='kmere',
        kmer_k=1,
        kmer_p_ins=0.001,
        kmer_p_del=0.001,
        kmer_p_sub=0.002,
        kmer_seed=42,
    ),

    # Test 3: Full error pipeline — synthesis + storage + sequencing
    # Higher c1 gives more correction power for the inner code.
    Params(
        name='end2end_gcplus_full_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='gcplus',
        binarization_method='default',
        sequence_length=200,
        gcplus_k=168,
        gcplus_l=8,
        gcplus_c1=4,

        # synthesis — create copies without errors
        synthesis_method='nosynthpoly',
        mean=30,
        std_dev=1,

        # storage
        storage_conditions='biogene',
        years=100,

        # sequencing errors
        sequencing_method='kmere',
        kmer_k=1,
        kmer_p_ins=0.004,
        kmer_p_del=0.004,
        kmer_p_sub=0.004,
        kmer_seed=42,
    ),
]


# Create a parameterized test case
def parameterized_test_generator(params):
    class ParameterizedTestBase(TestBase):
        def __init__(self, methodName='test_logic'):
            super().__init__(methodName, params=params)
    return ParameterizedTestBase


# Dynamically create test cases for each parameter set
def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for params in params_list:
        test_case = parameterized_test_generator(params)
        suite.addTest(test_case('test_logic'))
    return suite


if __name__ == '__main__':
    unittest.main()
