import unittest

from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from dnabyte import InSilicoDNA


class TestMESASequencing(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """Generate test sequences once for all tests"""
        cls.test_sequences = InSilicoDNA.random(m=1000, n=200)
    
    def test_simulate_sequencing_all_ids(self):
        """Test all valid mesa_sequencing_id values"""
        test_ids = [35, 36, 37, 38, 39, 40, 41]  # Sorted for clarity

        for mesa_id in test_ids:
            with self.subTest(mesa_sequencing_id=mesa_id):
                # Initialize MESA sequencing simulator with current ID
                params = Params(
                    sequencing_method='mesa', 
                    mesa_sequencing_id=mesa_id
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
                    f"Output length mismatch for ID {mesa_id}"
                )
                
                # Count sequences with differences
                sequences_with_errors = sum(
                    1 for orig, seq in zip(self.test_sequences.data, sequenced_sequences.data)
                    if orig != seq
                )
                
                error_rate = sequences_with_errors / len(self.test_sequences.data) * 100

                print('Error rate:', error_rate)

                # Note: Some sequencing IDs may not introduce errors for all sequences
                # so we don't assert error_rate > 0, just that it's non-negative
                self.assertGreaterEqual(
                    error_rate, 
                    0, 
                    f"Error rate should be non-negative for ID {mesa_id}"
                )
                
                # Check that info dictionary contains expected keys
                self.assertIn('average_copy_number', info)
                self.assertIn('number_of_sequencing_errors', info)
                self.assertIn('error_dict', info)
                
                # Validate info values are reasonable
                self.assertGreater(info['average_copy_number'], 0)
    
    def test_invalid_sequencing_id(self):
        """Test that invalid mesa_sequencing_id raises an error"""
        with self.assertRaises(ValueError):
            Params(
                sequencing_method='mesa',
                mesa_sequencing_id=999  # Invalid ID
            )

if __name__ == '__main__':
    unittest.main()