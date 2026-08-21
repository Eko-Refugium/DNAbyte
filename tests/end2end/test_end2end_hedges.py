import unittest
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets for Yin-Yang encoding
params_list = [
    # Test 1: Basic - no errors, just synthesis copies
    Params(
        name='end2end_yyc_basic',
        filename='Bohemian_Rhapsody_Lyrics.txt',

        # encoding parameters
        encoding_method='hedges',
        binarization_method='default',
        rs_num=0,
        mean=1,
        add_redundancy=False,
        add_primer=False,
        primer_length=50,
        sequence_length=300,
        max_homopolymer=3,

        # error channels
        storage_conditions=None,
        synthesis_method=None,
        clustering_method='kmere_cluster',
        kmer_size_cluster=250,
        recovery_method='debruijn'
        
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
