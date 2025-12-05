import unittest

from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from dnabyte import InSilicoDNA


class TestIIDSequencing(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """Generate test sequences once for all tests"""
        cls.test_sequences = InSilicoDNA.generate_random_sequences(m=1000, n=200)
    
    def test_simulate_sequencing_default_error_rate(self):
        """Test IID sequencing with default error rate (0.01)"""
        params = Params(sequencing_method='iid')
        simulator = SimulateSequencing(params=params)
        
        # Simulate sequencing
        sequenced_sequences, info = simulator.simulate(self.test_sequences)

        # Check that the output is InSilicoDNA
        self.assertIsInstance(sequenced_sequences, InSilicoDNA)
        
        # Verify output has same length as input
        self.assertEqual(
            len(sequenced_sequences.data), 
            len(self.test_sequences.data),
            "Output length mismatch"
        )
        
        # Count sequences with differences
        sequences_with_errors = sum(
            1 for orig, seq in zip(self.test_sequences.data, sequenced_sequences.data)
            if orig != seq
        )
        
        error_rate = sequences_with_errors / len(self.test_sequences.data) * 100

        # With 1% error rate and 1000 sequences of 200bp, we expect some errors
        self.assertGreater(
            error_rate, 
            0, 
            "Expected some sequencing errors with 1% error rate"
        )
        
        # Check that info dictionary exists
        self.assertIsInstance(info, dict)
    
    def test_simulate_sequencing_custom_error_rates(self):
        """Test IID sequencing with various custom error rates"""
        test_error_rates = [0.001, 0.005, 0.01, 0.05, 0.1]
        
        for error_rate in test_error_rates:
            with self.subTest(iid_error_rate=error_rate):
                params = Params(
                    sequencing_method='iid',
                    iid_error_rate=error_rate
                )
                simulator = SimulateSequencing(params=params)
                
                # Simulate sequencing
                sequenced_sequences, info = simulator.simulate(self.test_sequences)
                
                # Check that the output is InSilicoDNA
                self.assertIsInstance(sequenced_sequences, InSilicoDNA)
                
                # Verify output has same length as input
                self.assertEqual(
                    len(sequenced_sequences.data), 
                    len(self.test_sequences.data),
                    f"Output length mismatch for error rate {error_rate}"
                )
                
                # Check that params has the correct error rate
                self.assertEqual(params.iid_error_rate, error_rate)
    
    def test_zero_error_rate(self):
        """Test IID sequencing with zero error rate"""
        params = Params(
            sequencing_method='iid',
            iid_error_rate=0.0
        )
        simulator = SimulateSequencing(params=params)
        
        # Simulate sequencing
        sequenced_sequences, info = simulator.simulate(self.test_sequences)
        
        # With 0% error rate, sequences should be identical
        for orig, seq in zip(self.test_sequences.data, sequenced_sequences.data):
            self.assertEqual(orig, seq, "Sequences should be identical with 0% error rate")


if __name__ == '__main__':
    unittest.main()
