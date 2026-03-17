import unittest

from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from dnabyte import InSilicoDNA


class TestKMERSequencing(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Generate test sequences once for all tests."""
        cls.test_sequences = InSilicoDNA.random(m=200, n=120)

    def test_simulate_sequencing_default(self):
        """Test KMER sequencing simulation with default parameters."""
        params = Params(sequencing_method='kmere', kmer_seed=123)
        simulator = SimulateSequencing(params=params)

        sequenced_sequences, info = simulator.simulate(self.test_sequences)

        self.assertIsInstance(sequenced_sequences, InSilicoDNA)
        self.assertEqual(len(sequenced_sequences.data), len(self.test_sequences.data))

        self.assertIsInstance(info, dict)
        self.assertIn('number_of_sequencing_errors', info)

    def test_zero_error_rates_keep_sequences_identical(self):
        """If insertion/deletion/substitution rates are 0, output should match input."""
        params = Params(
            sequencing_method='kmere',
            kmer_p_ins=0.0,
            kmer_p_del=0.0,
            kmer_p_sub=0.0,
            kmer_seed=123,
        )
        simulator = SimulateSequencing(params=params)

        sequenced_sequences, _ = simulator.simulate(self.test_sequences)

        self.assertEqual(sequenced_sequences.data, self.test_sequences.data)

    def test_optional_custom_tables_allow_none(self):
        """Optional custom transition/substitution tables should allow None."""
        params = Params(
            sequencing_method='kmere',
            kmer_transition_probs=None,
            kmer_substitution_probs=None,
        )

        self.assertTrue(hasattr(params, 'kmer_transition_probs'))
        self.assertTrue(hasattr(params, 'kmer_substitution_probs'))
        self.assertIsNone(params.kmer_transition_probs)
        self.assertIsNone(params.kmer_substitution_probs)

    def test_only_one_custom_table_raises(self):
        """Providing only one custom table should raise ValueError."""
        with self.assertRaises(ValueError):
            Params(
                sequencing_method='kmere',
                kmer_transition_probs={},
            )

        with self.assertRaises(ValueError):
            Params(
                sequencing_method='kmere',
                kmer_substitution_probs={},
            )


if __name__ == '__main__':
    unittest.main()
