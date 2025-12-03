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
from dnabyte.binarization.text.binarize_text import TextBinarize
from dnabyte.binarization.compressed.binarize_compressed import ArchiveBinarize
from dnabyte.binarization.default.binarize_default import DefaultBinarize


class TestTextBinarize(unittest.TestCase):
    """Unit tests for TextBinarize class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.binarizer = Binarize(Params(binarization_method='binarize_text'))
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
        binarizer = Binarize(Params(binarization_method='binarize_text'))
        self.assertEqual(binarizer.text_encoding, 'utf-8')
    
    def test_init_custom_encoding(self):
        """Test initialization with custom encoding."""
        binarizer = Binarize(Params(binarization_method='binarize_text', text_encoding='ascii'))
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


class TestArchiveBinarize(unittest.TestCase):
    """Unit tests for ArchiveBinarize class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.binarizer = Binarize(Params(binarization_method='binarize_compressed'))
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

class TestBinarizationIntegration(unittest.TestCase):
    """Integration tests comparing different binarization methods."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        
        # Create a test text file that all methods can handle
        self.text_file = os.path.join(self.test_dir, 'test.txt')
        self.text_content = "Integration test content\nWith multiple lines\nAnd special chars: áéíóú"
        with open(self.text_file, 'w', encoding='utf-8') as f:
            f.write(self.text_content)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_all_methods_can_binarize_text_file(self):
        """Test that all three methods can binarize the same text file."""
        data = Data([self.text_file])
        
        # Test TextBinarize
        text_binarizer = Binarize(Params(binarization_method='binarize_text'))
        text_result = text_binarizer.binarize(data)
        self.assertIsInstance(text_result, BinaryCode)
        
        # Test ArchiveBinarize  
        archive_binarizer = Binarize(Params(binarization_method='binarize_compressed'))
        archive_result = archive_binarizer.binarize(data)
        self.assertIsInstance(archive_result, BinaryCode)
        
        # Test DefaultBinarize
        mock_params = MagicMock()
        default_binarizer = Binarize(Params(binarization_method='binarize_default'))
        default_result = default_binarizer.binarize(data)
        self.assertIsInstance(default_result, BinaryCode)
        
        # All should produce valid binary codes
        for result in [text_result, archive_result, default_result]:
            self.assertTrue(all(bit in '01' for bit in result.data))
            self.assertTrue(len(result.data) > 0)
    
    # def test_binary_lengths_comparison(self):
    #     """Compare binary lengths produced by different methods."""
    #     data = Data([self.text_file])
        
    #     text_binarizer = Binarize(Params(binarization_method='binarize_text'))
    #     archive_binarizer = Binarize(Params(binarization_method='binarize_compressed'))
    #     default_binarizer = Binarize(Params(binarization_method='binarize_default'))
        
    #     text_result = text_binarizer.binarize(data)
    #     archive_result = archive_binarizer.binarize(data)
    #     default_result = default_binarizer.binarize(data)
        
    #     # Print for comparison
    #     print(f"\nText result: {text_result}")
    #     print(f"Default result: {default_result}")
    #     print(f"Archive result: {archive_result}")
        
    #     # Print lengths
    #     print(f"\nText length: {len(text_result.data)}")
    #     print(f"Default length: {len(default_result.data)}")
    #     print(f"Archive length: {len(archive_result.data)}")
        
    #     # Print first 100 bits for comparison
    #     print(f"\nText binary (first 100 bits): {text_result.data}")
    #     print(f"Default binary (first 100 bits): {default_result.data}")


    #     # TextBinarize and DefaultBinarize should produce the same length
    #     # (both convert UTF-8 bytes directly to binary)
    #     self.assertEqual(len(text_result.data), len(default_result.data))
        
    #     # ArchiveBinarize might produce different length due to compression
    #     # (could be shorter due to compression, or longer due to archive overhead)
    #     self.assertTrue(len(archive_result.data) > 0)


if __name__ == '__main__':
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add all test classes
    test_classes = [
        TestTextBinarize,
        TestArchiveBinarize, 
        TestDefaultBinarize,
        TestBinarizationIntegration
    ]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        test_suite.addTests(tests)
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    # Print summary
    print(f"\n{'='*50}")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%")
    print(f"{'='*50}")