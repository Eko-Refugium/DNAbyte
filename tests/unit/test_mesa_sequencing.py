from dnabyte.sequencing.mesa.sequence import MESA
from dnabyte.params import Params
from dnabyte.data_classes.insilicodna import InSilicoDNA
import unittest

class TestMESA(unittest.TestCase):
    
    def setUp(self):
        # Initialize MESA sequencing simulator with default parameters
        params = Params(sequencing_method='mesa', mesa_sequencing_id=41)
        self.mesa_simulator = MESA(params=params)

    def test_simulate_sequencing(self):
        # Test data: list of DNA sequences
        test_sequences = ['ACGTAGCTAGCTAGCTAGCT', 'TGCATGCATGCATGCATGCA',
                          'GCTAGCTAGCTAGCTAGCTA', 'CATGCATGCATGCATGCATG']

        data = InSilicoDNA(test_sequences)

        # Simulate sequencing
        sequenced_sequences, info = self.mesa_simulator.simulate(data)

        print("Sequenced Sequences:", sequenced_sequences)

        # Check that the output is a list
        self.assertIsInstance(sequenced_sequences, list)

        # Check that the number of sequenced sequences is greater than or equal to input
        self.assertGreaterEqual(len(sequenced_sequences), len(test_sequences))

        # Check that info dictionary contains expected keys
        self.assertIn('average_copy_number', info)
        self.assertIn('number_of_sequencing_errors', info)
        self.assertIn('error_dict', info)