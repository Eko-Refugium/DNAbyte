import unittest

from dnabyte.store import SimulateStorage
from dnabyte.params import Params
from dnabyte import InSilicoDNA

class TestBiogene(unittest.TestCase):
    """ Test cases for the Biogene storage simulation. """

    def setUp(self):
        test_sequences = InSilicoDNA.generate_random_sequences(m=100, n=10000)
        self.data = InSilicoDNA(test_sequences)

    def test_simulate_storage_biogene(self):
        """Test biogene storage simulation"""
        
        # Initialize biogene storage simulator
        params = Params(storage_conditions='biogene', years=100)
        simulator = SimulateStorage(params=params)
        
        # Simulate storage
        stored_sequences, info = simulator.simulate(self.data)
        
        # Check that the output is InSilicoDNA
        self.assertIsInstance(stored_sequences, InSilicoDNA)
        
        # Check that the number of synthesized sequences is greater than or equal to input
        self.assertGreaterEqual(len(self.data.data), len(stored_sequences.data))
        
        # Check that info dictionary contains expected key
        self.assertIn('number_of_strand_breaks', info)
        
        # Verify strand breaks is non-negative
        self.assertGreaterEqual(info['number_of_strand_breaks'], 0)