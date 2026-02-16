import unittest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
import tarfile
import gzip

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.binarization.text.binarize import TextBinarize
from dnabyte.binarization.compressed.binarize import ArchiveBinarize
from dnabyte.binarization.default.binarize import DefaultBinarize


class TestTextBinarize(unittest.TestCase):
    """Unit tests for TextBinarize class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.binarizer = Binarize(Params(binarization_method='text', file_paths='textfile_40b.txt'))
        self.test_dir = tempfile.mkdtemp()
        
        # Create test text file
        self.text_content = "Hello, DNAbyte!\nThis is a test file with UTF-8 content: áéíóú"
        self.text_file = os.path.join(self.test_dir, 'test.txt')
        with open(self.text_file, 'w', encoding='utf-8') as f:
            f.write(self.text_content)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_init_default_encoding(self):
        """Test initialization with default encoding."""
        binarizer = Binarize(Params(binarization_method='text', file_paths='textfile_40b.txt'))
        self.assertEqual(binarizer.text_encoding, 'utf-8')
    
    def test_init_custom_encoding(self):
        """Test initialization with custom encoding."""
        binarizer = Binarize(Params(binarization_method='text', text_encoding='ascii', file_paths='textfile_40b.txt'))
        self.assertEqual(binarizer.text_encoding, 'ascii')
    
    def test_binarize_single_file_success(self):
        """Test successful binarization of a single text file."""
        data = Data([self.text_file])
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, [self.text_file])
        self.assertEqual(result.size, len(self.text_content.encode('utf-8')))
        self.assertTrue(len(result.data) > 0)
        self.assertTrue(all(bit in '01' for bit in result.data))
    
    def test_binarize_invalid_input_type(self):
        """Test binarization with invalid input type."""
        with self.assertRaises(TypeError):
            self.binarizer.binarize("not a data object")
    
    def test_binarize_multiple_files_error(self):
        """Test error when trying to binarize multiple files."""
        # Create second file
        second_file = os.path.join(self.test_dir, 'test2.txt')
        with open(second_file, 'w') as f:
            f.write("Second file")
        
        data = Data([self.text_file, second_file])
        with self.assertRaises(ValueError) as cm:
            self.binarizer.binarize(data)
        self.assertIn("exactly one file", str(cm.exception))
    
    def test_debinarize_success(self):
        """Test successful debinarization."""
        # First binarize
        data = Data([self.text_file])
        binary_code = self.binarizer.binarize(data)
        
        # Then debinarize
        output_file = os.path.join(self.test_dir, 'restored.txt')
        result = self.binarizer.debinarize(binary_code, output_file)
        
        self.assertTrue(result)
        self.assertTrue(os.path.exists(output_file))
        
        # Verify content
        with open(output_file, 'r', encoding='utf-8') as f:
            restored_content = f.read()
        self.assertEqual(restored_content, self.text_content)
    
    def test_debinarize_invalid_input(self):
        """Test debinarization with invalid input."""
        with self.assertRaises(TypeError):
            self.binarizer.debinarize("not a binary code", "output.txt")
    
    def test_debinarize_invalid_binary_length(self):
        """Test debinarization with invalid binary length."""
        # Create BinaryCode with invalid length (not divisible by 8)
        invalid_binary = BinaryCode(data="1010101", file_paths=[self.text_file])
        
        output_file = os.path.join(self.test_dir, 'output.txt')
        result = self.binarizer.debinarize(invalid_binary, output_file)
        
        self.assertFalse(result)
    
    def test_debinarize_default_output_path(self):
        """Test debinarization with default output path."""
        data = Data([self.text_file])
        binary_code = self.binarizer.binarize(data)
        
        result = self.binarizer.debinarize(binary_code)
        
        self.assertTrue(result)
        expected_path = os.path.join(self.test_dir, 'test_restored.txt')
        self.assertTrue(os.path.exists(expected_path))
    
    def test_roundtrip_binarization(self):
        """Test complete roundtrip: file -> binary -> file."""
        # Original file
        data = Data([self.text_file])
        
        # Binarize
        binary_code = self.binarizer.binarize(data)
        
        # Debinarize to new file
        restored_file = os.path.join(self.test_dir, 'roundtrip.txt')
        success = self.binarizer.debinarize(binary_code, restored_file)
        
        self.assertTrue(success)
        
        # Compare original and restored
        with open(self.text_file, 'r', encoding='utf-8') as f:
            original = f.read()
        with open(restored_file, 'r', encoding='utf-8') as f:
            restored = f.read()
        
        self.assertEqual(original, restored)


if __name__ == '__main__':
    unittest.main()