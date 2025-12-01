import unittest
from dnabyte.data_classes.binarycode import BinaryCode

class TestBinaryCode(unittest.TestCase):
    """Test cases for the BinaryCode class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.valid_bitstream = "10110010"
        self.long_bitstream = "1011001011010010110100101101001011010010"
        self.test_file_paths = ['/test/file1.txt', '/test/file2.txt']
    
    def test_valid_initialization(self):
        """Test successful initialization with valid bitstream."""
        binary = BinaryCode(self.valid_bitstream)
        
        self.assertEqual(binary.data, self.valid_bitstream)
        self.assertEqual(binary.length, len(self.valid_bitstream))
        self.assertEqual(binary.file_paths, [])
        self.assertIsNone(binary.size)
    
    def test_initialization_with_metadata(self):
        """Test initialization with file paths and original size."""
        binary = BinaryCode(
            self.valid_bitstream, 
            file_paths=self.test_file_paths, 
            size=1024
        )
        
        self.assertEqual(binary.file_paths, self.test_file_paths)
        self.assertEqual(binary.size, 1024)
        self.assertEqual(binary.data, self.valid_bitstream)
    
    def test_empty_bitstream_raises_error(self):
        """Test that empty bitstream raises ValueError."""
        with self.assertRaises(ValueError) as context:
            BinaryCode("")
        
        self.assertIn("Bitstream cannot be empty", str(context.exception))
    
    def test_non_string_bitstream_raises_error(self):
        """Test that non-string bitstream raises TypeError."""
        with self.assertRaises(TypeError) as context:
            BinaryCode(123)
        
        self.assertIn("Bitstream must be a string", str(context.exception))
    
    def test_invalid_characters_raise_error(self):
        """Test that invalid characters in bitstream raise ValueError."""
        test_cases = [
            ("101102", "{'2'}"),  # Contains '2'
            ("10112A", "{'2', 'A'}"),  # Contains '2' and 'A'
            ("xyz", "{'x', 'y', 'z'}"),  # Contains non-binary characters
        ]
        
        for bitstream, expected_chars in test_cases:
            with self.subTest(bitstream=bitstream):
                with self.assertRaises(ValueError) as context:
                    BinaryCode(bitstream)
                
                self.assertIn("contains invalid characters", str(context.exception))
    
    def test_mixed_valid_invalid_characters(self):
        """Test bitstream with mix of valid and invalid characters."""
        with self.assertRaises(ValueError) as context:
            BinaryCode("101012101")  # Contains '2'
        
        error_msg = str(context.exception)
        self.assertIn("contains invalid characters", error_msg)
        self.assertIn("{'2'}", error_msg)
    
    def test_len(self):
        """Test __len__ method."""
        binary = BinaryCode(self.valid_bitstream)
        self.assertEqual(len(binary), len(self.valid_bitstream))
        
        binary_long = BinaryCode(self.long_bitstream)
        self.assertEqual(len(binary_long), len(self.long_bitstream))
    
    def test_getitem(self):
        """Test __getitem__ method for indexing and slicing."""
        binary = BinaryCode(self.valid_bitstream)
        
        # Test individual indexing
        self.assertEqual(binary[0], '1')
        self.assertEqual(binary[1], '0')
        self.assertEqual(binary[-1], '0')
        
        # Test slicing
        self.assertEqual(binary[0:4], '1011')
        self.assertEqual(binary[2:6], '1100')
        self.assertEqual(binary[:3], '101')
        self.assertEqual(binary[5:], '010')
    
    def test_iter(self):
        """Test __iter__ method."""
        binary = BinaryCode("101")
        bits = list(binary)
        
        self.assertEqual(bits, ['1', '0', '1'])
        
        # Test with for loop
        result = []
        for bit in binary:
            result.append(bit)
        self.assertEqual(result, ['1', '0', '1'])
    
    def test_str_representation(self):
        """Test string representation."""
        binary = BinaryCode(self.valid_bitstream, 
                           file_paths=self.test_file_paths, 
                           size=100)
        str_repr = str(binary)
        
        self.assertIn("Type: BinaryCode", str_repr)
        self.assertIn("Length: 8 bits", str_repr)
        
        # Test without metadata
        binary_simple = BinaryCode(self.valid_bitstream)
        str_repr_simple = str(binary_simple)
        self.assertIn("Type: BinaryCode", str_repr_simple)
    
    def test_repr_representation(self):
        """Test repr representation."""
        binary = BinaryCode(self.valid_bitstream)
        repr_str = repr(binary)
        
        self.assertIn("BinaryCode", repr_str)
        self.assertIn("bitstream='10110010'", repr_str)
        self.assertIn("length=8", repr_str)
    
    def test_repr_long_bitstream(self):
        """Test repr representation with long bitstream."""
        binary = BinaryCode(self.long_bitstream)
        repr_str = repr(binary)
        
        self.assertIn("BinaryCode", repr_str)
        self.assertIn("...", repr_str)  # Should be truncated
        self.assertIn(f"length={len(self.long_bitstream)}", repr_str)
    
    def test_validate_bitstream_edge_cases(self):
        """Test validation with edge cases."""
        # Single character valid
        binary_single = BinaryCode("1")
        self.assertEqual(binary_single.data, "1")
        
        binary_single_zero = BinaryCode("0")
        self.assertEqual(binary_single_zero.data, "0")
        
        # Very long valid bitstream
        long_valid = "01" * 1000
        binary_long = BinaryCode(long_valid)
        self.assertEqual(len(binary_long), 2000)
    
    def test_attributes_consistency(self):
        """Test that all attributes are consistent after initialization."""
        bitstream = "11001010"
        file_paths = ['/path/to/test.txt']
        original_size = 512
        
        binary = BinaryCode(bitstream, file_paths, original_size)
        
        # Check all attributes are set correctly
        self.assertEqual(binary.data, bitstream)
        self.assertEqual(binary.length, len(bitstream))
        self.assertEqual(binary.file_paths, file_paths)
        self.assertEqual(binary.size, original_size)
        
        # Check consistency
        self.assertEqual(len(binary.data), binary.length)
        self.assertEqual(len(binary), binary.length)

if __name__ == '__main__':
    unittest.main()