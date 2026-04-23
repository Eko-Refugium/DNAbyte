import unittest
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets for Church encoding
params_list = [
    # Test 1: Basic - no errors, just synthesis copies
    Params(
        name='end2end_church_basic',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='church',
        binarization_method='default',
        sequence_length=200,
        max_homopolymer=2,
        rs_num=0,
        add_redundancy=True,
        add_primer=True,
        primer_length=20,

        # error channels
        storage_conditions=None,
        synthesis_method='nosynthpoly',
    ),

    # Test 2: Sequencing errors only (substitutions, insertions, deletions)
    Params(
        name='end2end_church_seq_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='church',
        binarization_method='default',
        sequence_length=200,
        max_homopolymer=4,
        rs_num=0,
        add_redundancy=True,
        add_primer=True,
        primer_length=20,

        # synthesis - create copies without errors
        synthesis_method='nosynthpoly',
        mean=10,
        std_dev=0,

        # storage - none
        storage_conditions=None,

        # sequencing errors
        sequencing_method='kmere',
        kmer_k=1,
        kmer_p_ins=0.001,
        kmer_p_del=0.001,
        kmer_p_sub=0.002,
        kmer_seed=42,
    ),

    # Test 3: Full error pipeline - synthesis + storage + sequencing + RS correction
    Params(
        name='end2end_church_full_errors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='church',
        binarization_method='default',
        sequence_length=200,
        max_homopolymer=4,
        rs_num=1,
        add_redundancy=True,
        add_primer=True,
        primer_length=20,

        # synthesis - create copies without errors
        synthesis_method='nosynthpoly',
        mean=30,
        std_dev=1,

        # storage
        storage_conditions='biogene',
        years=1,

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
