import unittest
import tempfile
import os
from dnabyte.data_classes.base import Data

class TestData(unittest.TestCase):
    """Test cases for the Data base class."""
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create temporary test files
        self.temp_dir = tempfile.mkdtemp()
        self.test_file1 = os.path.join(self.temp_dir, 'test1.txt')
        self.test_file2 = os.path.join(self.temp_dir, 'test2.txt')
        
        # Create small test files (under 1MB limit)
        with open(self.test_file1, 'w') as f:
            f.write('Hello World')
        with open(self.test_file2, 'w') as f:
            f.write('Test file 2')
    
    def tearDown(self):
        """Clean up after each test method."""
        # Remove temporary files
        for file_path in [self.test_file1, self.test_file2]:
            if os.path.exists(file_path):
                os.remove(file_path)
        if os.path.exists(self.temp_dir):
            os.rmdir(self.temp_dir)
    
    def test_valid_initialization(self):
        """Test successful initialization with valid file paths."""
        data = Data([self.test_file1, self.test_file2])
        
        self.assertEqual(data.file_paths, [self.test_file1, self.test_file2])
        self.assertEqual(len(data.file_paths), 2)
        self.assertGreater(data.size, 0)
    
    def test_empty_file_paths_raises_error(self):
        """Test that empty file_paths list raises ValueError."""
        with self.assertRaises(ValueError) as context:
            Data([])
        
        self.assertIn("file_paths cannot be empty", str(context.exception))
    
    def test_non_list_file_paths_raises_error(self):
        """Test that non-list file_paths raises TypeError."""
        with self.assertRaises(TypeError) as context:
            Data("not_a_list")
        
        self.assertIn("file_paths must be a list", str(context.exception))
    
    def test_non_string_file_path_raises_error(self):
        """Test that non-string file path raises TypeError."""
        with self.assertRaises(TypeError) as context:
            Data([self.test_file1, 123])
        
        self.assertIn("All file paths must be strings", str(context.exception))
    
    def test_non_existent_file_raises_error(self):
        """Test that non-existent file path raises ValueError."""
        with self.assertRaises(ValueError) as context:
            Data([self.test_file1, "/non/existent/file.txt"])
        
        self.assertIn("does not point to an existing file", str(context.exception))
    
    def test_file_size_limit_exceeded(self):
        """Test that files exceeding size limit raise ValueError."""
        large_file = os.path.join(self.temp_dir, 'large.txt')
        
        # Create file larger than 1MB
        with open(large_file, 'w') as f:
            f.write('x' * (Data.MAX_FILE_SIZE + 1))
        
        with self.assertRaises(ValueError) as context:
            Data([large_file])
        
        self.assertIn("exceeds maximum size limit", str(context.exception))
        
        # Clean up
        os.remove(large_file)
    
    def test_from_folder_valid(self):
        """Test from_folder class method with valid folder."""
        data = Data.from_folder(self.temp_dir)
        
        self.assertEqual(len(data.file_paths), 2)
        self.assertIn(self.test_file1, data.file_paths)
        self.assertIn(self.test_file2, data.file_paths)
    
    def test_from_folder_non_string_raises_error(self):
        """Test from_folder with non-string path raises TypeError."""
        with self.assertRaises(TypeError) as context:
            Data.from_folder(123)
        
        self.assertIn("folder_path must be a string", str(context.exception))
    
    def test_from_folder_non_existent_raises_error(self):
        """Test from_folder with non-existent folder raises ValueError."""
        with self.assertRaises(ValueError) as context:
            Data.from_folder("/non/existent/folder")
        
        self.assertIn("does not point to an existing directory", str(context.exception))
    
    def test_from_folder_empty_folder_raises_error(self):
        """Test from_folder with empty folder raises ValueError."""
        empty_dir = tempfile.mkdtemp()
        
        try:
            with self.assertRaises(ValueError) as context:
                Data.from_folder(empty_dir)
            
            self.assertIn("No files found in directory", str(context.exception))
        finally:
            os.rmdir(empty_dir)
    
    def test_calculate_total_bytes(self):
        """Test calculate_total_bytes method."""
        data = Data([self.test_file1, self.test_file2])
        
        expected_size = os.path.getsize(self.test_file1) + os.path.getsize(self.test_file2)
        self.assertEqual(data.calculate_total_bytes(), expected_size)
        self.assertEqual(data.size, expected_size)
    
    def test_str_representation(self):
        """Test string representation of Data object."""
        data = Data([self.test_file1])
        str_repr = str(data)
        
        self.assertIn("Type: Data", str_repr)
        self.assertIn("Total size:", str_repr)
        self.assertIn("Number of files: 1", str_repr)
        self.assertIn(self.test_file1, str_repr)

if __name__ == '__main__':
    unittest.main()