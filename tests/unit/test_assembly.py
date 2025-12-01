import unittest
from dnabyte.oligo import Oligo, complement
from dnabyte.oligopool import OligoPool
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.insilicodna import InSilicoDNA

from dnabyte.synthesize import SimulateSynthesis
import itertools
import random
import math

class TestAssembly(unittest.TestCase):

    def setUp(self):
        # create Data data object
        self.binary_code = BinaryCode('101010101')
        self.nucleobase_code = NucleobaseCode(self.binary_code)

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
                            block.append([(set1[index1], complement(rdm[i][1])), rdm[i], (complement(rdm[i][0]), set2[index2])])

                    else:
                            index = int(math.ceil(pos/2)) - 1
                            block.append([(complement(rdm[i][0]), set2[index]), rdm[i], (set1[index],complement(rdm[i][1]))])

            data.append(block)

        return data

    # test = [
    #         [('b*','A*'), ('b','B'), ('a','B*'), # pool 1
    #         ('c*','A*'), ('c','C'), ('a','C*'), 
    #         ('d*','A*'), ('d','D'), ('a','D*')],

    #         [('a*','C*'), ('c','C'), ('c*','B*'), # pool 2
    #         ('a*','D*'), ('d','D'), ('d*','B*'), 
    #         ('a*','E*'), ('e','E'), ('e*','B*')],

    #         [('c', 'B'), ('c*', 'C'), ('b', 'C*'), # pool 3
    #         ('d', 'B'), ('d*', 'D'), ('b', 'D*'),
    #         ('e', 'B'), ('e*', 'E'), ('b', 'E*')]
    # ]

    # test2 = [[[[('b*','A*'), ('b','B'), ('a','B*')], ('a*', 'F')], ('g', 'F*')], ('g*', 'C*'), ('d', 'C')]


    def test_assembly(self):
        # instantiate the SimulateSynthesis class
        assembly = SimulateSynthesis(mean=100, std_dev=2, hybridisation_steps=10000)
        assembled_data = assembly.simulate(self.nucleobase_code, mean=100, std_dev=2, hybridisation_steps=20000)
        self.assertIsInstance(assembled_data, InSilicoDNA)


if __name__ == '__main__':
    unittest.main()