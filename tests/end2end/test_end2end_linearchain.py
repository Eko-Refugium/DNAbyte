import unittest
from scipy.constants import Avogadro

from dnabyte.params import Params
from dnabyte.library import Library
from tests.testbase_end2end_newdata import TestBase

params_list = [
        Params(
        name='end2end_linearchain_noEC_noErrors',
        filename='Bohemian_Rhapsody_Lyrics.txt',

        # binarization parameters
        binarization_method='binarize_default',

        # encoding parameters
        encoding_method='linear_chain',
        library_name='20bp_Lib.csv',

        # error correction
        inner_error_correction=None,
        outer_error_correction=None,

        # error channels 
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None,
        theory='yes'
        )]
#     ),
#     Params(
#         name='end2end_linear_chain',
#         filename='Bohemian_Rhapsody_Lyrics.txt',

#         # encoding parameters
#         encoding_method='linearchain',
#         library_name='20bp_Lib.csv',
#         codeword_length=100,
#         dna_barcode_length=5,  
#         codeword_maxlength_positions=5,

#         # error correction
#         inner_error_correction='ltcode',
#         ltcode_header=4,
#         index_carry_length=3,
#         percent_of_symbols=2,
#         outer_error_correction='reedsolomon',
#         reed_solo_percentage=0.8,

#         # error channels
#         storage_conditions=None,
#         synthesis_method=None,
#         sequencing_method=None
#     ),
#         Params(
#         name='end2end_linear_chain',
#         filename='Bohemian_Rhapsody_Lyrics.txt',

#         # encoding parameters
#         encoding_method='linearchain',
#         library_name='20bp_Lib.csv',
#         codeword_length=100,
#         dna_barcode_length=5,  
#         codeword_maxlength_positions=5,

#         # error correction
#         inner_error_correction='ltcode',
#         ltcode_header=4,
#         index_carry_length=3,
#         percent_of_symbols=2,
#         outer_error_correction='reedsolomon',
#         reed_solo_percentage=0.8,

#         # synthesis error channel
#         mean=20,
#         vol=1000000 / Avogadro,
#         std_dev=1,
#         hybridisation_steps=10000,

#         # storage error channel
#         storage_conditions='biogene',
#         years=100,

#         # sequencing error channel
#         sequencing_method=38
#     )
# ]

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