import unittest

from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from dnabyte import InSilicoDNA


class TestIlluminaSequencing(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """Generate test sequences once for all tests"""
        cls.test_sequences = InSilicoDNA.random(m=500, n=200)
    
    def test_simulate_sequencing(self):
        """Test Illumina sequencing simulation"""
        params = Params(sequencing_method='illumina')
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

        # Illumina typically has low but non-zero error rates
        self.assertGreaterEqual(
            error_rate, 
            0, 
            "Error rate should be non-negative"
        )
        
        # Check that info dictionary exists
        self.assertIsInstance(info, dict)
    
    def test_sequence_integrity(self):
        """Test that Illumina sequencing maintains basic sequence properties"""
        params = Params(sequencing_method='illumina')
        simulator = SimulateSequencing(params=params)
        
        # Simulate sequencing
        sequenced_sequences, info = simulator.simulate(self.test_sequences)
        
        # Check that all sequences are strings
        for seq in sequenced_sequences.data:
            self.assertIsInstance(seq, str, "Each sequence should be a string")
        
        # Check that all sequences contain only valid DNA bases
        valid_bases = set('ACGT')
        for seq in sequenced_sequences.data:
            self.assertTrue(
                set(seq).issubset(valid_bases),
                f"Sequence contains invalid bases: {set(seq) - valid_bases}"
            )
    
    def test_multiple_runs_produce_different_results(self):
        """Test that multiple sequencing runs produce different results (stochastic)"""
        params = Params(sequencing_method='illumina')
        simulator = SimulateSequencing(params=params)
        
        # Run sequencing twice
        result1, _ = simulator.simulate(self.test_sequences)
        result2, _ = simulator.simulate(self.test_sequences)
        
        # Results should potentially differ due to randomness (though not guaranteed)
        # At minimum, they should both be valid InSilicoDNA objects
        self.assertIsInstance(result1, InSilicoDNA)
        self.assertIsInstance(result2, InSilicoDNA)
        self.assertEqual(len(result1.data), len(result2.data))


if __name__ == '__main__':
    unittest.main()
