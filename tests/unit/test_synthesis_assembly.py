"""
Comprehensive unit tests for assembly synthesis module.
"""

import unittest
from unittest.mock import Mock, MagicMock, patch
from dnabyte.synthesis.assembly.oligo import Oligo, complement
from dnabyte.synthesis.assembly.oligopool import OligoPool
from dnabyte import BinaryCode, NucleobaseCode, InSilicoDNA
from dnabyte.params import Params
from dnabyte.synthesize import SimulateSynthesis
import itertools
import random
import math

class TestAssembly(unittest.TestCase):

    def setUp(self):
        self.nucleobase_code = NucleobaseCode.random(type='linear_assembly', library='20bp_Lib.csv', m=10, n=5)

    def test_assembly(self):
        params = Params(encoding_method='linear_chain', 
                        codeword_length=100, 
                        dna_barcode_length=20, 
                        assembly_structure='assembly', 
                        synthesis_method='assembly', 
                        mean=100, 
                        std_dev=2, 
                        hybridisation_steps=10000, 
                        library_name='20bp_Lib.csv')
        assembly = SimulateSynthesis(params=params)
        assembled_data, info = assembly.simulate(self.nucleobase_code)
        self.assertIsInstance(assembled_data, InSilicoDNA)

    def test_assembly_linear_binom(self):

        self.nucleobase_code = NucleobaseCode.random(type='linear_binom', library='20bp_Lib.csv', m=15, n=5)

        params = Params(encoding_method='linear_binom', 
                        sigma=5,
                        assembly_structure='assembly', 
                        synthesis_method='assembly', 
                        percent_of_symbols=5,
                        mean=100, 
                        std_dev=2, 
                        hybridisation_steps=10000, 
                        library_name='20bp_Lib.csv')
        assembly = SimulateSynthesis(params=params)
        assembled_data, info = assembly.simulate(self.nucleobase_code)
        print(assembled_data.data)

        self.assertIsInstance(assembled_data, InSilicoDNA)

    def test_assembly_poly_chain(self):
        """Test assembly with poly_chain encoding structure"""
        self.nucleobase_code = NucleobaseCode.random(type='poly_chain', library='polymeraselibinparts.txt', m=10, n=5)

        params = Params(encoding_method='poly_chain', 
                        sigma=5,
                        assembly_structure='assembly', 
                        synthesis_method='assembly', 
                        percent_of_symbols=5,
                        mean=100, 
                        std_dev=2, 
                        hybridisation_steps=10000, 
                        library_name='polymeraselibinparts.txt')
        assembly = SimulateSynthesis(params=params)
        assembled_data, info = assembly.simulate(self.nucleobase_code)
        print(f"Poly chain assembled: {len(assembled_data.data)} sequences")

        self.assertIsInstance(assembled_data, InSilicoDNA)

    def test_assembly_poly_binom(self):
        """Test assembly with poly_binom encoding structure"""
        self.nucleobase_code = NucleobaseCode.random(type='poly_binom', library='polymeraselibinparts.txt', m=15, n=5)

        params = Params(encoding_method='poly_binom', 
                        sigma=5,
                        assembly_structure='assembly', 
                        synthesis_method='assembly', 
                        percent_of_symbols=5,
                        mean=100, 
                        std_dev=2, 
                        hybridisation_steps=10000, 
                        library_name='polymeraselibinparts.txt')
        assembly = SimulateSynthesis(params=params)
        assembled_data, info = assembly.simulate(self.nucleobase_code)
        print(f"Poly binom assembled: {len(assembled_data.data)} sequences")

        self.assertIsInstance(assembled_data, InSilicoDNA)

    


if __name__ == '__main__':
    unittest.main()