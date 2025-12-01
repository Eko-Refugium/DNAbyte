import unittest
from dnabyte.data_classes.insilicodna import InSilicoDNA

class TestInSilicoDNA(unittest.TestCase):
    """Test cases for the InSilicoDNA class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.valid_sequences = [
            "ATCGATCG",
            "GCTAGCTA", 
            "TTAATTAA",
            "CCGGCCGG"
        ]
        
        self.single_sequence = ["ATCG"]
        self.mixed_length_sequences = [
            "AT",
            "GCTAGC", 
            "TTAATTAACCGG"
        ]
    
    def test_valid_initialization(self):
        """Test successful initialization with valid DNA sequences."""
        dna = InSilicoDNA(self.valid_sequences)
        
        self.assertEqual(dna.data, self.valid_sequences)
        self.assertEqual(dna.num_sequences, 4)
        self.assertEqual(dna.total_length, 32)  # 8 * 4 sequences
        self.assertEqual(dna.average_length, 8.0)
        self.assertEqual(dna.file_paths, [])
        self.assertIsNone(dna.size)
    
    def test_initialization_single_sequence(self):
        """Test initialization with single sequence."""
        dna = InSilicoDNA(self.single_sequence)
        
        self.assertEqual(dna.num_sequences, 1)
        self.assertEqual(dna.total_length, 4)
        self.assertEqual(dna.average_length, 4.0)
    
    def test_initialization_mixed_lengths(self):
        """Test initialization with sequences of different lengths."""
        dna = InSilicoDNA(self.mixed_length_sequences)
        
        self.assertEqual(dna.num_sequences, 3)
        self.assertEqual(dna.total_length, 20)  # 2 + 6 + 12
        self.assertEqual(dna.average_length, 20/3)
    
    def test_empty_data_raises_error(self):
        """Test that empty data raises ValueError."""
        with self.assertRaises(ValueError) as context:
            InSilicoDNA([])
        
        self.assertIn("DNA data cannot be empty", str(context.exception))
    
    def test_non_list_data_raises_error(self):
        """Test that non-list data raises TypeError."""
        with self.assertRaises(TypeError) as context:
            InSilicoDNA("ATCG")
        
        self.assertIn("DNA data must be a list", str(context.exception))
    
    def test_non_string_sequence_raises_error(self):
        """Test that non-string sequences raise ValueError."""
        invalid_data = ["ATCG", 123, "GCTA"]
        
        with self.assertRaises(ValueError) as context:
            InSilicoDNA(invalid_data)
        
        self.assertIn("DNA sequence at index 1 must be a string", str(context.exception))
    
    def test_empty_sequence_raises_error(self):
        """Test that empty sequences raise ValueError."""
        invalid_data = ["ATCG", "", "GCTA"]
        
        with self.assertRaises(ValueError) as context:
            InSilicoDNA(invalid_data)
        
        self.assertIn("DNA sequence at index 1 cannot be empty", str(context.exception))
    
    def test_invalid_nucleotides_raise_error(self):
        """Test that sequences with invalid nucleotides raise ValueError."""
        test_cases = [
            (["ATCG", "GXTC"], "{'X'}"),  # Contains X
            (["ATCG", "GCNTA"], "{'N'}"),  # Contains N
            (["ATCG", "GC123"], "{'1', '2', '3'}"),  # Contains numbers
            (["ATCG", "gc-ta"], "{'-', 'c', 'g'}"),  # Contains lowercase and dash
        ]
        
        for invalid_data, expected_chars in test_cases:
            with self.subTest(sequences=invalid_data):
                with self.assertRaises(ValueError) as context:
                    InSilicoDNA(invalid_data)
                
                error_msg = str(context.exception)
                self.assertIn("contains invalid nucleotides", error_msg)
    
    def test_get_sequence(self):
        """Test get_sequence method."""
        dna = InSilicoDNA(self.valid_sequences)
        
        self.assertEqual(dna.get_sequence(0), "ATCGATCG")
        self.assertEqual(dna.get_sequence(1), "GCTAGCTA")
        self.assertEqual(dna.get_sequence(3), "CCGGCCGG")
    
    def test_get_sequence_out_of_range(self):
        """Test get_sequence with out of range indices."""
        dna = InSilicoDNA(self.valid_sequences)
        
        with self.assertRaises(IndexError):
            dna.get_sequence(4)
        
        with self.assertRaises(IndexError):
            dna.get_sequence(-1)
    
    def test_add_sequence_valid(self):
        """Test adding valid sequences."""
        dna = InSilicoDNA(self.valid_sequences.copy())
        original_count = dna.num_sequences
        
        dna.add_sequence("AAATTTGGG")
        
        self.assertEqual(dna.num_sequences, original_count + 1)
        self.assertEqual(dna.data[-1], "AAATTTGGG")
        self.assertGreater(dna.total_length, 32)
    
    def test_add_sequence_invalid(self):
        """Test adding invalid sequences raises errors."""
        dna = InSilicoDNA(self.valid_sequences.copy())
        
        # Test non-string
        with self.assertRaises(ValueError):
            dna.add_sequence(123)
        
        # Test empty string
        with self.assertRaises(ValueError):
            dna.add_sequence("")
        
        # Test invalid nucleotides
        with self.assertRaises(ValueError):
            dna.add_sequence("ATCGX")
    
    def test_remove_sequence(self):
        """Test removing sequences."""
        dna = InSilicoDNA(self.valid_sequences.copy())
        original_count = dna.num_sequences
        original_sequence = dna.get_sequence(1)
        
        dna.remove_sequence(1)
        
        self.assertEqual(dna.num_sequences, original_count - 1)
        self.assertNotEqual(dna.get_sequence(1), original_sequence)
    
    def test_remove_sequence_out_of_range(self):
        """Test removing sequence with invalid index."""
        dna = InSilicoDNA(self.valid_sequences)
        
        with self.assertRaises(IndexError):
            dna.remove_sequence(10)
    
    def test_get_sequence_lengths(self):
        """Test get_sequence_lengths method."""
        dna = InSilicoDNA(self.mixed_length_sequences)
        lengths = dna.get_sequence_lengths()
        
        expected_lengths = [2, 6, 12]
        self.assertEqual(lengths, expected_lengths)
    
    def test_get_nucleotide_counts(self):
        """Test get_nucleotide_counts method."""
        # Use sequences with known composition
        test_sequences = ["AAAA", "CCCC", "GGGG", "TTTT"]
        dna = InSilicoDNA(test_sequences)
        
        counts = dna.get_nucleotide_counts()
        
        expected_counts = {'A': 4, 'C': 4, 'G': 4, 'T': 4}
        self.assertEqual(counts, expected_counts)
    
    def test_get_nucleotide_counts_mixed(self):
        """Test nucleotide counts with mixed sequences."""
        test_sequences = ["ATCG", "AATT"]  # A:3, T:3, C:1, G:1
        dna = InSilicoDNA(test_sequences)
        
        counts = dna.get_nucleotide_counts()
        
        expected_counts = {'A': 3, 'T': 3, 'C': 1, 'G': 1}
        self.assertEqual(counts, expected_counts)
    
    def test_validate(self):
        """Test validate method."""
        dna = InSilicoDNA(self.valid_sequences)
        
        self.assertTrue(dna.validate())
        
        # Test validation after manual data corruption
        dna.data[0] = "ATCGX"  # Invalid nucleotide
        with self.assertRaises(ValueError):
            dna.validate()
    
    def test_len(self):
        """Test __len__ method."""
        dna = InSilicoDNA(self.valid_sequences)
        self.assertEqual(len(dna), 4)
        
        dna_single = InSilicoDNA(self.single_sequence)
        self.assertEqual(len(dna_single), 1)
    
    def test_getitem(self):
        """Test __getitem__ method."""
        dna = InSilicoDNA(self.valid_sequences)
        
        self.assertEqual(dna[0], "ATCGATCG")
        self.assertEqual(dna[1], "GCTAGCTA")

        print('DNA:')
        print(dna[-1])

        self.assertEqual(dna[-1], "CCGGCCGG")
    
    def test_iter(self):
        """Test __iter__ method."""
        dna = InSilicoDNA(self.valid_sequences)
        sequences = list(dna)
        
        self.assertEqual(sequences, self.valid_sequences)
        
        # Test with for loop
        result = []
        for seq in dna:
            result.append(seq)
        self.assertEqual(result, self.valid_sequences)
    
    def test_str_representation(self):
        """Test string representation."""
        dna = InSilicoDNA(self.valid_sequences)
        str_repr = str(dna)
        
        self.assertIn("Type: InSilicoDNA", str_repr)
        self.assertIn("Number of sequences: 4", str_repr)
        self.assertIn("Total length: 32 nucleotides", str_repr)
        self.assertIn("Average length: 8.0 nucleotides", str_repr)
        self.assertIn("Nucleotide distribution:", str_repr)
    
    def test_metrics_update_after_modifications(self):
        """Test that metrics are updated after adding/removing sequences."""
        dna = InSilicoDNA(["ATCG"])
        
        # Initial state
        self.assertEqual(dna.num_sequences, 1)
        self.assertEqual(dna.total_length, 4)
        self.assertEqual(dna.average_length, 4.0)
        
        # After adding sequence
        dna.add_sequence("GGCCAA")
        self.assertEqual(dna.num_sequences, 2)
        self.assertEqual(dna.total_length, 10)
        self.assertEqual(dna.average_length, 5.0)
        
        # After removing sequence
        dna.remove_sequence(0)
        self.assertEqual(dna.num_sequences, 1)
        self.assertEqual(dna.total_length, 6)
        self.assertEqual(dna.average_length, 6.0)
    
    def test_edge_cases(self):
        """Test various edge cases."""
        # Very short sequences
        short_sequences = ["A", "T", "G", "C"]
        dna_short = InSilicoDNA(short_sequences)
        self.assertEqual(dna_short.average_length, 1.0)
        
        # Very long sequence
        long_sequence = ["A" * 1000]
        dna_long = InSilicoDNA(long_sequence)
        self.assertEqual(dna_long.total_length, 1000)

if __name__ == '__main__':
    unittest.main()