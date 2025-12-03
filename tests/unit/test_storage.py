import unittest

from dnabyte.store import SimulateStorage
from dnabyte.params import Params
from dnabyte import InSilicoDNA

class TestBiogene(unittest.TestCase):
    """ Test cases for the Biogene storage simulation. """

    def setUp(self):
        self.data = InSilicoDNA.generate_random_sequences(m=100, n=10000)

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

    def test_simulate_storage_permafrost(self):
        """Test permafrost storage simulation"""
        
        # Initialize permafrost storage simulator with fewer years
        # (decay rate is 5.5E-6/nt/yr, higher than biogene)
        params = Params(storage_conditions='permafrost', years=10)
        simulator = SimulateStorage(params=params)
        
        # Simulate storage
        stored_sequences, info = simulator.simulate(self.data)
        
        # Check that the output is InSilicoDNA
        self.assertIsInstance(stored_sequences, InSilicoDNA)
        
        # Check that stored sequences are less than or equal to input
        self.assertGreaterEqual(len(self.data.data), len(stored_sequences.data))
        
        # Check that info dictionary contains expected key
        self.assertIn('number_of_strand_breaks', info)
        
        # Verify strand breaks is non-negative
        self.assertGreaterEqual(info['number_of_strand_breaks'], 0)

    def test_simulate_storage_newstorage(self):
        """Test newstorage (cryogenic) storage simulation"""
        
        # Initialize newstorage simulator
        params = Params(storage_conditions='newstorage', years=100)
        simulator = SimulateStorage(params=params)
        
        # Simulate storage
        stored_sequences, info = simulator.simulate(self.data)
        
        # Check that the output is InSilicoDNA
        self.assertIsInstance(stored_sequences, InSilicoDNA)
        
        # Check that stored sequences are less than or equal to input
        self.assertGreaterEqual(len(self.data.data), len(stored_sequences.data))
        
        # Check that info dictionary contains expected key
        self.assertIn('number_of_strand_breaks', info)
        
        # Verify strand breaks is non-negative
        self.assertGreaterEqual(info['number_of_strand_breaks'], 0)

    def test_simulate_storage_roomtemperature(self):
        """Test room temperature storage simulation"""
        
        # Initialize room temperature storage simulator
        # (half-life = 521 years, lambda = 0.00133/nt/yr)
        params = Params(storage_conditions='roomtemperature', years=0.1)
        simulator = SimulateStorage(params=params)
        
        # Simulate storage
        stored_sequences, info = simulator.simulate(self.data)
        
        # Check that the output is InSilicoDNA
        self.assertIsInstance(stored_sequences, InSilicoDNA)
        
        # Check that stored sequences are less than or equal to input
        self.assertGreaterEqual(len(self.data.data), len(stored_sequences.data))
        
        # Check that info dictionary contains expected key
        self.assertIn('number_of_strand_breaks', info)
        
        # Verify strand breaks is non-negative
        self.assertGreaterEqual(info['number_of_strand_breaks'], 0)

    # def test_simulate_storage_random(self):
    #     """Test random storage simulation (for testing purposes)"""
        
    #     # Initialize random storage simulator
    #     params = Params(storage_conditions='random', years=100)
    #     simulator = SimulateStorage(params=params)
        
    #     # Simulate storage
    #     stored_sequences, info = simulator.simulate(self.data)
        
    #     # Check that the output is a list (random storage may return different format)
    #     self.assertIsInstance(stored_sequences, (list, InSilicoDNA))
        
    #     # Check that info dictionary contains expected key
    #     self.assertIn('number_of_strand_breaks', info)
        
    #     # Verify strand breaks is non-negative
    #     self.assertGreaterEqual(info['number_of_strand_breaks'], 0)


if __name__ == '__main__':
    unittest.main()


