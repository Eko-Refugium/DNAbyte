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
        data = [[['ACGTTACTAGT', 'TGCATGCAAGTC'], 'CGTACGTAGCTA'],
                [['GGTACCGTAAGC', 'CCATGGCATGCA'], 'TTAGGCATCGTA'],
                [['TTAACCGGTAGC', 'AAGGTTCCAAGT'], 'GGCCAATTGGCA']]

        self.nucleobase_code = NucleobaseCode(data)

        # Define motif pair dictionary for complement function
        self.motif_pair = {
            'a': 'a*', 'b': 'b*', 'c': 'c*', 'd': 'd*', 'e': 'e*', 'f': 'f*', 'g': 'g*', 'h': 'h*',
            'a*': 'a', 'b*': 'b', 'c*': 'c', 'd*': 'd', 'e*': 'e', 'f*': 'f', 'g*': 'g', 'h*': 'h',
            'A': 'A*', 'B': 'B*', 'C': 'C*', 'D': 'D*', 'E': 'E*', 'F': 'F*', 'G': 'G*', 'H': 'H*',
            'A*': 'A', 'B*': 'B', 'C*': 'C', 'D*': 'D', 'E*': 'E', 'F*': 'F', 'G*': 'G', 'H*': 'H'
        }

        self.nested_list = self.create_random_messages(32, 5)
        self.nucleobase_code.data = self.nested_list


    def create_random_messages(self, positions, k):
        """
        Create k variants of random messages at a given position. 
        """
        if positions > 32:
            raise ValueError("The number of positions must be lower or equal to 32.")
        
        if k > 256:
            raise ValueError("The number of messages must be below 256.")

        data = []

        set1 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'a*', 'b*', 'c*', 'd*', 'e*', 'f*', 'g*', 'h*']
        set2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'A*', 'B*', 'C*', 'D*', 'E*', 'F*', 'G*', 'H*']
        cartesian_product = list(itertools.product(set1, set2))
        
        for pos in range(1, positions + 1):
            rdm = random.sample(cartesian_product, k)

            block = []
            for i in range(k):

                    if pos%2 == 0:
                            index1 = (int(pos/2) + 7) % 16
                            index2 = (int(pos/2) + 8) % 16
                            block.append([(set1[index1], complement(rdm[i][1], self.motif_pair)), rdm[i], (complement(rdm[i][0], self.motif_pair), set2[index2])])

                    else:
                            index = int(math.ceil(pos/2)) - 1
                            block.append([(complement(rdm[i][0], self.motif_pair), set2[index]), rdm[i], (set1[index], complement(rdm[i][1], self.motif_pair))])

            data.append(block)

        return data

    def test_assembly(self):
        # instantiate the SimulateSynthesis class
        # TODO: fix this test with new params structure
        params = Params(encoding_method='linear_chain', 
                        codeword_length=200, 
                        dna_barcode_length=20, 
                        assembly_structure='assembly', 
                        synthesis_method='assembly', 
                        mean=100, 
                        std_dev=2, 
                        hybridisation_steps=10000, 
                        library_name='20bp_Lib.csv')
        assembly = SimulateSynthesis(params=params)
        assembled_data = assembly.simulate(self.nucleobase_code)
        self.assertIsInstance(assembled_data, InSilicoDNA)


if __name__ == '__main__':
    unittest.main()