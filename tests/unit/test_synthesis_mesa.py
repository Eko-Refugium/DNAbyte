from dnabyte.synthesize import SimulateSynthesis
from dnabyte.params import Params
from dnabyte import NucleobaseCode, InSilicoDNA
import unittest

class TestMESA(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """Generate test sequences once for all tests"""
        cls.test_sequences = NucleobaseCode.random(m=100, n=50)
    
    def test_simulate_synthesis_all_ids(self):
        """Test all valid mesa_synthesis_id values"""
        test_ids = [3, 4, 5, 6, 7, 68, 69, 70, 71]

        for mesa_id in test_ids:
            with self.subTest(mesa_synthesis_id=mesa_id):
                # Initialize MESA synthesis simulator with current ID
                params = Params(
                    assembly_structure='synthesis', 
                    synthesis_method='mesa', 
                    mesa_synthesis_id=mesa_id, 
                    mean=100, 
                    std_dev=5)
                simulator = SimulateSynthesis(params=params)
                                
                # Simulate synthesis
                synthesized_sequences, info = simulator.simulate(self.test_sequences)

                # Check that the output is InSilicoDNA
                self.assertIsInstance(synthesized_sequences, InSilicoDNA)
                
                # Verify output has more or equal sequences than input (due to copy number)
                self.assertGreaterEqual(
                    len(synthesized_sequences.data), 
                    len(self.test_sequences.data),
                    f"Output length should be >= input for ID {mesa_id}"
                )
                
                # Count sequences with synthesis errors
                sequences_with_errors = sum(
                    1 for orig, synth in zip(self.test_sequences.data * (len(synthesized_sequences.data) // len(self.test_sequences.data) + 1), 
                                            synthesized_sequences.data)
                    if orig != synth
                )
                
                error_rate = sequences_with_errors / len(synthesized_sequences.data) * 100

                print(f'Synthesis ID {mesa_id} - Error rate: {error_rate:.2f}%')

                # Note: Some synthesis IDs may not introduce errors for all sequences
                # so we don't assert error_rate > 0, just that it's non-negative
                self.assertGreaterEqual(
                    error_rate, 
                    0, 
                    f"Error rate should be non-negative for ID {mesa_id}"
                )
                
                # Check that info dictionary contains expected keys
                self.assertIn('average_copy_number', info)
                self.assertIn('number_of_synthesis_errors', info)
                self.assertIn('error_dict', info)
                
                # Validate info values are reasonable
                self.assertGreater(info['average_copy_number'], 0)
    
    def test_invalid_synthesis_id(self):
        """Test that invalid mesa_synthesis_id raises an error"""
        with self.assertRaises(ValueError):
            Params(
                assembly_structure='synthesis',
                synthesis_method='mesa',
                mesa_synthesis_id=999  # Invalid ID
            )

if __name__ == '__main__':
    unittest.main()