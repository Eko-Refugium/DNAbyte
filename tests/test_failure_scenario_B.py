import unittest
from scipy.constants import Avogadro

from dnabyte.params import Params
from dnabyte.library import Library
from tests.testbase_end2end import TestBase

params_list = [
    Params(
        name='FAILURE_SCENARIO_A_MAXDENSITY',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='synthesis',
        encoding_scheme='max_density_encoding',
        library_name='',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=40,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no',
        seed=42
    ),

Params(
        name='FAILURE_SCENARIO_A_NOHOMOPOLYMER',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='synthesis',
        encoding_scheme='no_homopolymeroddeven_encoding',
        library_name='',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=40,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no',
        seed=42
    ),
    
Params(
        name='FAILURE_SCENARIO_A_LINEARCHAIN',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='linear_assembly',
        encoding_scheme='linear_encoding',
        library_name='20bp_Lib.csv',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=40,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no',
        seed=42
    ),
Params(
        name='FAILURE_SCENARIO_A_LINEARBINOMIAL',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='linear_assembly',
        encoding_scheme='binomial_encoding',
        library_name='20bp_Lib.csv',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=40,
        reed_solo_percentage=0.8,
        sigma_amount=3,
        theory='no',
        seed=42
    ),
Params(
        name='FAILURE_SCENARIO_A_POSITIONALLINEAR',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='positional_assembly',
        encoding_scheme='linear_encoding',
        library_name='polymeraselibinparts_200(8)_150(40)_50.csv',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=None,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no',
        seed=42
    ),
Params(
        name='FAILURE_SCENARIO_A_POSITIONALBINOMIAL',
        filename='Bohemian_Rhapsody_Lyrics.txt',
        assembly_structure='positional_assembly',
        encoding_scheme='binomial_encoding',
        library_name='polymeraselibinparts_200(8)_150(40)_50.csv',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=100,
        storage_conditions='permafrost',
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=None,
        sequencing_method=40,
        reed_solo_percentage=0.8,
        sigma_amount=3,
        theory='no',
        seed=42
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