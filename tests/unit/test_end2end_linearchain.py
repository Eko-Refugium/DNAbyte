import unittest
from scipy.constants import Avogadro

from dnabyte.params import Params
from dnabyte.library import Library
from tests.testbase_end2end_newdata import TestBase

params_list = [
        Params(
        name='end2end_linearchain_noEC_noErrors',
        filename='textfile_40b.txt',

        # binarization parameters
        binarization_method='default',

        # encoding parameters
        encoding_method='linear_chain',
        library_name='20bp_Lib.csv',
        synthesis_method = 'assembly',
        codeword_length=30,
        mean=1,
        sdt_dev=0,

        # error correction

        # error channels 
        # storage_conditions=['permafrost', 'biogene', 'newstorage', 'random','roomtemperature'],
        # years=100,
        theory='no'
       
    ),
    Params(
        name='end2end_linear_chain',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='linear_chain',
        library_name='20bp_Lib.csv',
        codeword_length=30,
        dna_barcode_length=2,  
        codeword_maxlength_positions=3,

        # error correction
        inner_error_correction='ltcode',
        ltcode_header=1,
        index_carry_length=1,
        percent_of_symbols=3,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        
        binarization_method='default',
        mean=20,
        std_dev=0,

        # error channels
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None,
        theory='no'
    ),
    Params(
        name='end2end_linear_chain',
        filename='textfile_40b.txt',

        # encoding parameters
        encoding_method='linear_chain',
        library_name='20bp_Lib.csv',
        codeword_length=30,
        dna_barcode_length=2,  
        codeword_maxlength_positions=3,

        # error correction
        inner_error_correction='ltcode',
        ltcode_header=1,
        index_carry_length=1,
        percent_of_symbols=3,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,

        
        binarization_method='default',
        mean=20,
        std_dev=0,

        # error channels
        storage_conditions='biogene',
        years=100,
        sequencing_method='mesa',
        mesa_sequencing_id = 41,
        theory='no'
    )
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