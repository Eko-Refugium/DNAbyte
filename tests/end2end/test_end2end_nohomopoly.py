import unittest
from scipy.constants import Avogadro

from dnabyte.params import Params
from dnabyte.library import Library
from tests.testbase_end2end_newdata import TestBase

params_list = [
    Params(
        name='end2end_nohomopoly_noEC_noErrors',
        file_paths=['Bohemian_Rhapsody_Lyrics.txt'],

        # binarization method
        binarization_method='compressed',

        # encoding parameters
        encoding_method='no_homopolymer',

        # error correction
        inner_error_correction=None,
        outer_error_correction=None,

        # error channels 
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None,
    ),
    Params(
        name='end2end_nohomopoly_allEC_noErrors',
        file_paths=['Bohemian_Rhapsody_Lyrics.txt'],

        # binarization method
        binarization_method='compressed',

        # encoding parameters
        encoding_method='no_homopolymer',
        codeword_length=501,
        dna_barcode_length=34,

        # error correction
        inner_error_correction='ltcode',
        ltcode_header=34,
        percent_of_symbols=4,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.9,

        # error channels
        storage_conditions=None,
        synthesis_method=None,
        sequencing_method=None
    )]#,
    # Params(
    #     name='end2end_synthesis_nohomopoly_allErrors',
    #     filename='Bohemian_Rhapsody_Lyrics.txt',
        
    #     # encoding parameters
    #     encoding_method='nohomopoly',
    #     codeword_length=501,
    #     dna_barcode_length=34, 
    #     codeword_maxlength_positions=18,

    #     # error correction
    #     inner_error_correction='ltcode',
    #     percent_of_symbols=2,
    #     ltcode_header=34,
    #     index_carry_length=35,

    #     outer_error_correction='reedsolomon',
    #     reed_solo_percentage=0.8,

    #     # synthesis error channel
    #     synthesis_method=68,
    #     mean=20,
    #     vol=1000000 / Avogadro,
    #     std_dev=1,
    #     hybridisation_steps=10000,

    #     # storage error channel
    #     storage_conditions='biogene',
    #     years=100,

    #     # sequencing error channel
    #     sequencing_method=38

    # ),
    # Params(
    #     name='end2end_synthesis_nohomopoly',
    #     filename='Bohemian_Rhapsody_Lyrics.txt',
    #     assembly_structure='synthesis',
    #     encoding_scheme='no_homopolymeroddeven_encoding',
    #     library_name='',
    #     mean=100,
    #     vol=1000000 / Avogadro,
    #     std_dev=1,
    #     hybridisation_steps=10000,
    #     inner_error_correction='ltcode',
    #     outer_error_correction='reedsolomon',
    #     dna_barcode_length=34,  # in bp here not in oligos
    #     codeword_maxlength_positions=18,  # in bp here not in oligos
    #     years=100,
    #     storage_conditions='Biogene',
    #     codeword_length=501,  # in bp here not in oligos
    #     percent_of_symbols=2,
    #     index_carry_length=34,  # in bp here not in oligos
    #     synthesis_method=68,
    #     sequencing_method=40,
    #     reed_solo_percentage=0.8,
    #     sigma_amount=None,
    #     theory='no'
    # )


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