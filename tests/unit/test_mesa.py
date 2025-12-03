from dnabyte.synthesis.mesa.synthesize import MESA
from dnabyte.params import Params
from dnabyte.data_classes.insilicodna import InSilicoDNA
import unittest

class TestMESA(unittest.TestCase):
    
    def setUp(self):
        # Initialize MESA synthesis simulator with default parameters
        params = Params(assembly_structure='synthesis', synthesis_method='mesa', mesa_synthesis_id=3, mean=100, std_dev=5)
        self.mesa_simulator = MESA(params=params)

    def test_simulate_synthesis(self):
        # Test data: list of DNA sequences
        test_sequences = ['ACGTAGCTAGCTAGCTAGCT', 'TGCATGCATGCATGCATGCA']

        # Simulate synthesis
        synthesized_sequences, info = self.mesa_simulator.simulate(test_sequences)

        # Check that the output is a list
        self.assertIsInstance(synthesized_sequences, list)

        # Check that the number of synthesized sequences is greater than or equal to input
        self.assertGreaterEqual(len(synthesized_sequences), len(test_sequences))

        # Check that info dictionary contains expected keys
        self.assertIn('average_copy_number', info)
        self.assertIn('number_of_synthesis_errors', info)
        self.assertIn('error_dict', info)