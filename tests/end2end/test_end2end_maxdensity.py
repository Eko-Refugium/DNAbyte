# Prevent unittest from auto-discovering this module
__test__ = False

import unittest
from scipy.constants import Avogadro
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets
params_list = [
    Params(
        name='end2end_maxdensity_noEC_noErrors',
        file_paths=['Bohemian_Rhapsody_Lyrics.txt'],

        # binarization method
        binarization_method='compressed',

        # encoding parameters
        encoding_method='max_density',

        # error correction
        inner_error_correction=None,
        outer_error_correction=None,

        # error channels 
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None
    ),
    Params(
        name='end2end_maxdensity_allEC_noErrors',
        file_paths=['Bohemian_Rhapsody_Lyrics.txt'],

        # binarization method
        binarization_method='compressed',

        # encoding parameters
        encoding_method='max_density',

        codeword_length=501,
        dna_barcode_length=34,

        # error correction
        inner_error_correction='ltcode',
        ltcode_header=34,
        percent_of_symbols=2,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        # error channels
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None,

    ),
    Params(
        name='end2end_maxdensity_allErrors',
        file_paths=['Bohemian_Rhapsody_Lyrics.txt'],

        # binarization method
        binarization_method='compressed',

        # encoding parameters
        encoding_method='max_density',
        codeword_length=501,
        dna_barcode_length=34,
        codeword_maxlength_positions=18,  # in bp here not in oligos

        # error correction
        inner_error_correction='ltcode',
        percent_of_symbols=2,
        ltcode_header=34,

        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        # synthesis error channel
        synthesis_method='mesa',
        synthesis_mesa_id=68,
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,

        # storage error channel
        storage_conditions='biogene',
        years=10,

        # sequencing error channel
        sequencing_method='mesa',
        sequencing_mesa_id=68,

    )]

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