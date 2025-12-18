import unittest
from scipy.constants import Avogadro
from dnabyte.params import Params

from tests.testbase_end2end_newdata import TestBase

# Define different parameter sets
params_list = [
    # Params(
    #     name='end2end_maxdensity_noEC_noErrors',
    #     filename='textfile_40b.txt',

    #     # encoding parameters
    #     encoding_method='max_density',

    #     binarization_method='default',

    #     sequencing_method='mesa',

    #     mesa_sequencing_id=41,

    #     # error correction


    # ),
    # Params(
    #     name='end2end_maxdensity_allEC_noErrors',
    #     filename='Bohemian_Rhapsody_Lyrics.txt',

    #     # encoding parameters
    #     encoding_method='max_density',
        
    #     binarization_method='default',
    #     codeword_length=501,
    #     dna_barcode_length=34,

    #     # error correction
    #     inner_error_correction='ltcode',
    #     ltcode_header=34,
    #     percent_of_symbols=2,
    #     outer_error_correction='reedsolomon',
    #     reed_solo_percentage=0.8,

    #     # error channels
    #     storage_conditions=None,
    #     synthesis_method=None,
    #     sequencing_method=None,

    # ),
    Params(
        name='end2end_maxdensity_allErrors',
        filename='Bohemian_Rhapsody_Lyrics.txt',

        # encoding parameters
        encoding_method='max_density',
        codeword_length=1000,
        dna_barcode_length=10,
        codeword_maxlength_positions=4,  # in bp here not in oligos
        
        binarization_method='default',
        # error correction
        inner_error_correction=None,
        percent_of_symbols=3,
        ltcode_header=15,
        index_carry_length=15,

        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        # # synthesis error channel
        # synthesis_method='mesa',
        # mesa_synthesis_id=3,
        # mean=20,
        # std_dev=2,

        # # storage error channel
        # storage_conditions='biogene',
        # years=1000000,

        # # sequencing error channel
        # sequencing_method='mesa',
        # mesa_sequencing_id=41

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