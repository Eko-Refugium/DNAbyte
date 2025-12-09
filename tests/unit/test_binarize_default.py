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


class TestDefaultBinarize(unittest.TestCase):
    """Unit tests for DefaultBinarize class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create mock params
        self.binarizer = Binarize(Params(binarization_method='binarize_default'))
        self.test_dir = tempfile.mkdtemp()
        
        # Create test files of different types
        self.text_file = os.path.join(self.test_dir, 'test.txt')
        with open(self.text_file, 'w') as f:
            f.write("Test text content")
        
        self.binary_file = os.path.join(self.test_dir, 'test.bin')
        with open(self.binary_file, 'wb') as f:
            f.write(b'\x00\x01\x02\x03\xFF\xFE\xFD')
        
        # Create an image-like file (just bytes that could represent an image)
        self.image_file = os.path.join(self.test_dir, 'test.jpg')
        with open(self.image_file, 'wb') as f:
            # Write some bytes that could represent a simple image header
            f.write(b'\xFF\xD8\xFF\xE0')  # JPEG header-like
            f.write(b'fake image data' * 10)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_init_with_params(self):
        """Test initialization with params."""
        self.assertEqual(self.binarizer.name, 'binarize_default')
    
    def test_binarize_text_file_success(self):
        """Test successful binarization of text file."""
        data = Data([self.text_file])
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, [self.text_file])
        self.assertEqual(result.size, os.path.getsize(self.text_file))
        self.assertTrue(len(result.data) > 0)
        self.assertTrue(all(bit in '01' for bit in result.data))
        # Check that binary length matches file size * 8
        self.assertEqual(len(result.data), os.path.getsize(self.text_file) * 8)
    
    def test_binarize_binary_file_success(self):
        """Test successful binarization of binary file."""
        data = Data([self.binary_file])
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, [self.binary_file])
        self.assertEqual(result.size, os.path.getsize(self.binary_file))
        self.assertTrue(len(result.data) > 0)
    
    def test_binarize_image_file_success(self):
        """Test successful binarization of image-like file."""
        data = Data([self.image_file])
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, [self.image_file])
        self.assertEqual(result.size, os.path.getsize(self.image_file))
    
    def test_binarize_invalid_input_type(self):
        """Test binarization with invalid input type."""
        with self.assertRaises(TypeError):
            self.binarizer.binarize("not a data object")
    
    def test_binarize_multiple_files_error(self):
        """Test error when trying to binarize multiple files."""
        data = Data([self.text_file, self.binary_file])
        with self.assertRaises(ValueError) as cm:
            self.binarizer.binarize(data)
        self.assertIn("exactly one file", str(cm.exception))
    
    def test_debinarize_text_file_success(self):
        """Test successful debinarization of text file."""
        # First binarize
        data = Data([self.text_file])
        binary_code = self.binarizer.binarize(data)
        
        # Then debinarize
        output_file = os.path.join(self.test_dir, 'restored.txt')
        result = self.binarizer.debinarize(binary_code, output_file)
        
        self.assertTrue(result)
        self.assertTrue(os.path.exists(output_file))
        
        # Verify content matches exactly
        with open(self.text_file, 'rb') as f:
            original_bytes = f.read()
        with open(output_file, 'rb') as f:
            restored_bytes = f.read()
        
        self.assertEqual(original_bytes, restored_bytes)
    
    def test_debinarize_binary_file_success(self):
        """Test successful debinarization of binary file."""
        # First binarize
        data = Data([self.binary_file])
        binary_code = self.binarizer.binarize(data)
        
        # Then debinarize
        output_file = os.path.join(self.test_dir, 'restored.bin')
        result = self.binarizer.debinarize(binary_code, output_file)
        
        self.assertTrue(result)
        self.assertTrue(os.path.exists(output_file))
        
        # Verify binary content matches exactly
        with open(self.binary_file, 'rb') as f:
            original_bytes = f.read()
        with open(output_file, 'rb') as f:
            restored_bytes = f.read()
        
        self.assertEqual(original_bytes, restored_bytes)
    
    def test_debinarize_invalid_input(self):
        """Test debinarization with invalid input."""
        with self.assertRaises(TypeError):
            self.binarizer.debinarize("not a binary code")
    
    def test_debinarize_invalid_binary_length(self):
        """Test debinarization with invalid binary length."""
        invalid_binary = BinaryCode(data="1010101", file_paths=[self.text_file])
        
        output_file = os.path.join(self.test_dir, 'output.bin')
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
    
    def test_roundtrip_various_file_types(self):
        """Test roundtrip for various file types."""
        test_files = [self.text_file, self.binary_file, self.image_file]
        
        for original_file in test_files:
            with self.subTest(file=original_file):
                # Binarize
                data = Data([original_file])
                binary_code = self.binarizer.binarize(data)
                
                # Debinarize
                base_name = os.path.basename(original_file)
                restored_file = os.path.join(self.test_dir, f'restored_{base_name}')
                success = self.binarizer.debinarize(binary_code, restored_file)
                
                self.assertTrue(success)
                
                # Compare byte-by-byte
                with open(original_file, 'rb') as f:
                    original_bytes = f.read()
                with open(restored_file, 'rb') as f:
                    restored_bytes = f.read()
                
                self.assertEqual(original_bytes, restored_bytes)

if __name__ == '__main__':
    unittest.main()