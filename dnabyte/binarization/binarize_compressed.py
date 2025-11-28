from datetime import datetime
import os
import tarfile
import tempfile
import datetime
from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.binarize import Binarize

class ArchiveBinarize(Binarize):
    """
    ArchiveBinarize creates a compressed tar.gz archive from multiple files 
    and converts it to binary representation.
    
    This binarizer can handle multiple files and folders by creating a compressed archive.
    It's ideal when you need to process multiple files while maintaining file structure.
    """
    
    def __init__(self, params):
        """
        Initialize CompressedBinarize with compression settings.
        
        :param compression_level: Compression level (1-9, where 9 is highest compression)
        :param temp_dir: Directory for temporary files (uses system temp if None)
        """
        self.params = params
        self.name = 'archive'

    def binarize(self, data):
        """
        Converts files in the Data object to a binary sequence by creating a compressed archive.

        :param data: A Data object containing file paths
        :return: A BinaryCode object with binary sequence
        :raises TypeError: If data is not a Data object
        :raises ValueError: If data cannot be processed
        :raises FileNotFoundError: If any file cannot be found
        """
        # Validate input
        if not isinstance(data, Data):
            raise TypeError(f"Expected Data object, got {type(data).__name__}")
        
        if not data.file_paths:
            raise ValueError("Data object contains no file paths")
        
        # Validate all files exist before processing
        self._validate_files_exist(data.file_paths)
        
        archive_path = None

        try:
            folder_path = os.path.dirname(data.file_paths[0])
            archive_name = folder_path + '/archive.tar.gz'

            with tarfile.open(archive_name, 'w:gz') as archive:
                for file_path in data.file_paths:
                    archive.add(file_path)
            
            size = os.path.getsize(archive_name)

            with open(archive_name, 'rb') as file:
                binary_data = file.read()
                binary_sequence = ''.join(format(byte, '08b') for byte in binary_data)

        except FileNotFoundError as e:
            raise FileNotFoundError(f"File not found: {e}")
        except Exception as e:
            raise ValueError(f"Error during binarization: {e}")
        finally:
            os.remove(archive_name)

            

        return binary_sequence, size
    
    def _validate_files_exist(self, file_paths):
        """
        Validate that all files exist before processing.
        
        :param file_paths: List of file paths to validate
        :raises FileNotFoundError: If any file doesn't exist
        """
        for file_path in file_paths:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
    
    
    
    def debinarize(self, data, output_directory=None):
        """
        Restores the original files from the binary data by recreating the archive 
        and extracting its contents to a folder named with job_identifier.

        :param data: A BinaryCode object containing the binary sequence
        :param output_directory: Base directory to create job folder in (uses current directory if None)
        :return: Boolean indicating success of restoration
        """
        temp_archive_path = None
        try:
            # Validate input
            if not isinstance(data, BinaryCode):
                raise TypeError(f"Expected BinaryCode object, got {type(data).__name__}")
            
            if len(data.data) % 8 != 0:
                raise ValueError("Binary sequence length must be divisible by 8")

            binary_data = data.data
            self.file_paths = data.file_paths
            job_identifier = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')  

            # Set default output directory if not provided
            if output_directory is None:
                output_directory = os.getcwd()
            
            # Create job-specific folder
            job_folder = os.path.join(output_directory, job_identifier)
            os.makedirs(job_folder, exist_ok=True)

            # Convert binary sequence back to bytes
            byte_data = bytes(int(binary_data[i:i+8], 2) for i in range(0, len(binary_data), 8))

            # Create a temporary tar archive file
            temp_archive_path = os.path.join(job_folder, 'temp_archive.tar.gz')
            with open(temp_archive_path, 'wb') as file:
                file.write(byte_data)

            # Extract the tar archive to the job folder, flattening the directory structure
            with tarfile.open(temp_archive_path, 'r:gz') as archive:
                for member in archive.getmembers():
                    if member.isfile():
                        # Extract only the filename (without path) and place it directly in job_folder
                        member.name = os.path.basename(member.name)
                        archive.extract(member, path=job_folder)

            print(f"Files successfully restored to: {job_folder}")
            return True
            
        except Exception as e:
            print(f"Restoration failed: {e}")
            return False
        finally:
            # Clean up temporary archive
            if temp_archive_path and os.path.exists(temp_archive_path):
                try:
                    os.remove(temp_archive_path)
                except:
                    pass  # Don't fail if cleanup fails

    def __str__(self):
        """String representation of CompressedBinarize."""
        return f"CompressedBinarize(compression_level={self.compression_level}, temp_dir='{self.temp_dir}')"
    
    def __repr__(self):
        """Developer representation of CompressedBinarize."""
        return f"CompressedBinarize(compression_level={self.compression_level}, temp_dir={self.temp_dir!r})"