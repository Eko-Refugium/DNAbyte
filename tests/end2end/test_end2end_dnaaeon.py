import unittest
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets for DNA-Aeon encoding.
# DNA-Aeon uses the NOREC4DNA RU10 Raptor fountain code for DNA storage.
# The fountain code provides erasure resilience — some oligos can be lost
# and the data can still be recovered.
params_list = [
    # Test 1: Basic — no errors, just synthesis copies
    Params(
        name='end2end_dna_aeon_basic',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='dna_aeon',
        binarization_method='default',
        sequence_length=200,
        dna_aeon_chunk_size=10,
        dna_aeon_overhead=0.40,
        dna_aeon_insert_header=False,
        dna_aeon_error_correction='crc',
        dna_aeon_use_dna_rules=True,
        dna_aeon_drop_upper_bound=0.5,

        # error channels
        storage_conditions=None,
        synthesis_method='nosynthpoly',
    ),

    # Test 2: Sequencing errors only
    Params(
        name='end2end_dna_aeon_seq_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='dna_aeon',
        binarization_method='default',
        sequence_length=200,
        dna_aeon_chunk_size=10,
        dna_aeon_overhead=0.40,
        dna_aeon_insert_header=False,
        dna_aeon_error_correction='crc',
        dna_aeon_use_dna_rules=True,
        dna_aeon_drop_upper_bound=0.5,

        # synthesis — create copies without errors
        synthesis_method='nosynthpoly',
        mean=100,
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
    Params(
        name='end2end_dna_aeon_full_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='dna_aeon',
        binarization_method='default',
        sequence_length=200,
        dna_aeon_chunk_size=10,
        dna_aeon_overhead=0.80,
        dna_aeon_insert_header=False,
        dna_aeon_error_correction='crc',
        dna_aeon_use_dna_rules=True,
        dna_aeon_drop_upper_bound=0.5,

        # synthesis — create copies without errors
        synthesis_method='nosynthpoly',
        mean=100,
        std_dev=1,

        # storage
        storage_conditions='biogene',
        years=100,

        # sequencing errors
        sequencing_method='kmere',
        kmer_k=1,
        kmer_p_ins=0.002,
        kmer_p_del=0.002,
        kmer_p_sub=0.002,
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
