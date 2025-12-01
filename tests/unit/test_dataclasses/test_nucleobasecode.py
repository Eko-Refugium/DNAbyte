import unittest
from dnabyte.data_classes.nucleobasecode import NucleobaseCode


class TestNucleobaseCode(unittest.TestCase):
    """Unit tests for NucleobaseCode class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Simple codeword structure
        self.simple_data = [['ATCG', 'CGTA'], ['GGCC', 'AATT']]
        
        # Complex nested structure
        self.complex_data = [
            [['ATG', 'CCC'], ['GGG', 'TAA']],
            [['ATCG'], ['CGAT', 'TACG']]
        ]
        
        # Single level structure
        self.single_level = [['ATCG'], ['CGTA']]
    
    def test_init_valid_simple_structure(self):
        """Test initialization with valid simple structure."""
        nucleobase = NucleobaseCode(self.simple_data)
        self.assertEqual(nucleobase.data, self.simple_data)
        self.assertEqual(nucleobase.num_codewords, 2)
        self.assertEqual(nucleobase.total_elements, 8)  # 8 individual nucleotides
    
    def test_init_valid_complex_structure(self):
        """Test initialization with valid complex nested structure."""
        nucleobase = NucleobaseCode(self.complex_data)
        self.assertEqual(nucleobase.data, self.complex_data)
        self.assertEqual(nucleobase.num_codewords, 2)
        self.assertTrue(nucleobase.total_elements > 0)
        self.assertTrue(nucleobase.max_depth > 0)
    
    def test_init_empty_data(self):
        """Test initialization with empty data raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            NucleobaseCode([])
        self.assertIn("cannot be empty", str(cm.exception))
    
    def test_init_non_list_data(self):
        """Test initialization with non-list data raises TypeError."""
        with self.assertRaises(TypeError) as cm:
            NucleobaseCode("not a list")
        self.assertIn("must be a list", str(cm.exception))
    
    def test_get_codeword_valid_index(self):
        """Test getting codeword with valid index."""
        nucleobase = NucleobaseCode(self.simple_data)
        first_codeword = nucleobase.get_codeword(0)
        self.assertEqual(first_codeword, ['ATCG', 'CGTA'])
    
    def test_get_codeword_invalid_index(self):
        """Test getting codeword with invalid index raises IndexError."""
        nucleobase = NucleobaseCode(self.simple_data)
        with self.assertRaises(IndexError) as cm:
            nucleobase.get_codeword(5)
        self.assertIn("out of range", str(cm.exception))
    
    def test_get_codeword_negative_index(self):
        """Test getting codeword with negative index raises IndexError."""
        nucleobase = NucleobaseCode(self.simple_data)
        with self.assertRaises(IndexError):
            nucleobase.get_codeword(-1)
    
    def test_indexing_functionality(self):
        """Test indexing using [] operator."""
        nucleobase = NucleobaseCode(self.simple_data)
        self.assertEqual(nucleobase[0], ['ATCG', 'CGTA'])
        self.assertEqual(nucleobase[1], ['GGCC', 'AATT'])
    
    def test_length_functionality(self):
        """Test len() functionality."""
        nucleobase = NucleobaseCode(self.simple_data)
        self.assertEqual(len(nucleobase), 2)
    
    def test_count_total_elements(self):
        """Test total elements counting."""
        nucleobase = NucleobaseCode(self.simple_data)
        # Each string counts as 1 element in this implementation
        self.assertEqual(nucleobase.total_elements, 4)  # 4 strings total
    
    def test_calculate_max_depth_simple(self):
        """Test max depth calculation for simple structure."""
        nucleobase = NucleobaseCode(self.simple_data)
        # Structure: [['ATCG', 'CGTA'], ['GGCC', 'AATT']]
        # Depth: list -> list -> string = depth 2
        self.assertEqual(nucleobase.max_depth, 2)
    
    def test_calculate_max_depth_complex(self):
        """Test max depth calculation for complex structure."""
        nucleobase = NucleobaseCode(self.complex_data)
        # Structure has 3 levels: list -> list -> list -> string = depth 3
        self.assertEqual(nucleobase.max_depth, 2)
    
    def test_validate_method(self):
        """Test validate method returns True for valid data."""
        nucleobase = NucleobaseCode(self.simple_data)
        self.assertTrue(nucleobase.validate())
    
    def test_string_representation(self):
        """Test string representation."""
        nucleobase = NucleobaseCode(self.simple_data)
        str_repr = str(nucleobase)
        
        self.assertIn("Type: NucleobaseCode", str_repr)
        self.assertIn("Number of codewords: 2", str_repr)
        self.assertIn("Total elements:", str_repr)
        self.assertIn("Max depth:", str_repr)
        self.assertIn("DATA:", str_repr)
    
    def test_string_representation_with_file_paths(self):
        """Test string representation with file paths."""
        nucleobase = NucleobaseCode(self.simple_data)
        nucleobase.file_paths = ['/path/to/file.txt']
        nucleobase.size = 100
        
        str_repr = str(nucleobase)
        self.assertIn("File paths:", str_repr)
        self.assertIn("Original size: 100 bytes", str_repr)
    
    def test_empty_codeword_structure(self):
        """Test handling of empty codeword structures."""
        empty_codeword_data = [[], ['ATCG']]
        nucleobase = NucleobaseCode(empty_codeword_data)
        self.assertEqual(nucleobase.num_codewords, 2)
        # Should handle empty codewords gracefully
        self.assertTrue(nucleobase.total_elements >= 0)



if __name__ == '__main__':
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add all test classes
    test_classes = [TestNucleobaseCode]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        test_suite.addTests(tests)
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    # Print summary
    print(f"\n{'='*60}")
    print("DATA CLASSES TEST SUMMARY")
    print(f"{'='*60}")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%")
    print(f"{'='*60}")