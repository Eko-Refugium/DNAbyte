import unittest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock


from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.binarization.text.binarize import TextBinarize
from dnabyte.binarization.compressed.binarize import ArchiveBinarize
from dnabyte.binarization.default.binarize import DefaultBinarize



class TestArchiveBinarize(unittest.TestCase):
    """Unit tests for ArchiveBinarize class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.binarizer = Binarize(Params(binarization_method='compressed'))
        self.test_dir = tempfile.mkdtemp()
        
        # Create multiple test files
        self.files = []
        for i in range(3):
            file_path = os.path.join(self.test_dir, f'test{i}.txt')
            with open(file_path, 'w') as f:
                f.write(f"Test content {i}\nLine 2 of file {i}")
            self.files.append(file_path)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_binarize_multiple_files_success(self):
        """Test successful binarization of multiple files."""
        data = Data(self.files)
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, self.files)
        self.assertTrue(len(result.data) > 0)
        self.assertTrue(all(bit in '01' for bit in result.data))
    
    def test_binarize_single_file_success(self):
        """Test successful binarization of single file."""
        data = Data([self.files[0]])
        result = self.binarizer.binarize(data)
        
        self.assertIsInstance(result, BinaryCode)
        self.assertEqual(result.file_paths, [self.files[0]])
    
    def test_binarize_invalid_input_type(self):
        """Test binarization with invalid input type."""
        with self.assertRaises(TypeError):
            self.binarizer.binarize("not a data object")
    
    def test_debinarize_success(self):
        """Test successful debinarization."""
        # First binarize
        data = Data(self.files)
        binary_code = self.binarizer.binarize(data)
        
        # Then debinarize
        output_dir = os.path.join(self.test_dir, 'output')
        result = self.binarizer.debinarize(binary_code, output_dir)
        
        self.assertTrue(result)
        
        # Check that job folder was created
        job_folders = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
        self.assertEqual(len(job_folders), 1)
        
        job_folder = os.path.join(output_dir, job_folders[0])
        
        # Check that files were restored
        restored_files = os.listdir(job_folder)
        expected_files = [os.path.basename(f) for f in self.files]
        self.assertEqual(set(restored_files), set(expected_files))
    
    def test_debinarize_invalid_input(self):
        """Test debinarization with invalid input."""
        with self.assertRaises(TypeError):
            self.binarizer.debinarize("not a binary code")
    
    def test_debinarize_invalid_binary_length(self):
        """Test debinarization with invalid binary length."""
        invalid_binary = BinaryCode(data="1010101", file_paths=self.files)
        
        result = self.binarizer.debinarize(invalid_binary)
        
        self.assertFalse(result)
    
    def test_roundtrip_binarization(self):
        """Test complete roundtrip: files -> archive -> binary -> archive -> files."""
        # Original files
        data = Data(self.files)
        
        # Binarize
        binary_code = self.binarizer.binarize(data)
        
        # Debinarize
        output_dir = os.path.join(self.test_dir, 'roundtrip')
        success = self.binarizer.debinarize(binary_code, output_dir)
        
        self.assertTrue(success)
        
        # Find the job folder
        job_folders = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
        job_folder = os.path.join(output_dir, job_folders[0])
        
        # Compare original and restored files
        for original_file in self.files:
            filename = os.path.basename(original_file)
            restored_file = os.path.join(job_folder, filename)
            
            with open(original_file, 'r') as f:
                original_content = f.read()
            with open(restored_file, 'r') as f:
                restored_content = f.read()
            
            self.assertEqual(original_content, restored_content)

if __name__ == '__main__':
    unittest.main()