import unittest
from scipy.constants import Avogadro

from dnabyte.params import Params
from dnabyte.library import Library
from tests.testbase_end2end_newdata import TestBase

params_list = [
    Params(
        name='end2end_polybinom_noEC_noErrors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='poly_binom',
        synthesis_method = 'assembly',
        library_name='lib_positional_22(173)m_35(2)g_24(39)p.csv',

        mean=10,
        std_dev=0,

        # error correction

        binarization_method='default',
        # dna_barcode_length=2,
        # error channels 
        theory='no',
    
    ),
    Params(
        name='end2end_polybinom_noEC_noErrors',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='poly_binom',
        library_name='polymeraselibinparts.txt',
        synthesis_method = 'assembly',

        # error correction
        inner_error_correction='ltcode',
        ltcode_header=4,
        index_carry_length=3,
        percent_of_symbols=2,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        # error channels 
        
        binarization_method='default',
        storage_conditions=None,
        sequencing_method=None,
    ),
    Params(
        name='end2end_positional_binom',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='poly_binom',
        library_name='polymeraselibinparts.txt',

        # error correction
        inner_error_correction=None,

        ltcode_header=4,
        index_carry_length=3,

        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.7,
        binarization_method='default',

        # error channels

        mean=1,
        vol=1000000 / Avogadro,
        std_dev=0,
        hybridisation_steps=10,
        
        dna_barcode_length=5,  
        codeword_maxlength_positions=5,
        years=0,
        storage_conditions='biogene',
        codeword_length=100,
        percent_of_symbols=2,
        
        theory='no',
        synthesis_method=None,
        sequencing_method=None,
        sigma_amount=3,
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