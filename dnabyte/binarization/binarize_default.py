import os
from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.data_classes.binarycode import BinaryCode

def attributes(params):
    if 'filename' not in params.__dict__ or params.filename is None:
        raise ValueError("filename parameter must be specified")
    else:
        filename = params.filename
        
    return {"filename": filename}

class DefaultBinarize(Binarize):
    """
    DefaultBinarize converts a single file of any type to binary representation by reading raw bytes.
    
    This is a universal binarizer that can handle any single file type (pictures, audio, 
    executables, etc.) by directly converting file bytes to binary.
    
    Note: This binarizer can only handle exactly one file at a time, as multiple files
    cannot be reliably separated during restoration without additional metadata.
    """
    
    def __init__(self, params):
        """
        Initialize DefaultBinarize.
        """
        self.params = params
        self.name = 'default'
    
    def binarize(self, data):
        """
        Converts a single file in the Data object to a binary sequence by reading raw bytes.

        :param data: A Data object containing exactly one file path
        :return: A BinaryCode object with binary sequence
        :raises TypeError: If data is not a Data object
        :raises ValueError: If data contains more than one file or cannot be processed
        :raises FileNotFoundError: If the file cannot be found
        """
        # Validate input
        if not isinstance(data, Data):
            raise TypeError(f"Expected Data object, got {type(data).__name__}")
        
        if not data.file_paths:
            raise ValueError("Data object contains no file paths")
        
        # Check that there's exactly one file
        if len(data.file_paths) != 1:
            raise ValueError(f"DefaultBinarize can only handle exactly one file, got {len(data.file_paths)} files. "
                           f"For multiple files, use CompressedBinarize instead.")
        
        file_path = data.file_paths[0]
        
        try:
            # Binarize the single file
            binary_sequence, size = self._binarize_single_file(file_path)
            #print(binary_sequence[:64], 'binary sequence snippet in binarize_default') 
            
            # Create and return BinaryCode object
            return binary_sequence, size
            
            
        except FileNotFoundError as e:
            raise FileNotFoundError(f"File not found: {e}")
        except Exception as e:
            raise ValueError(f"Error during binarization: {e}")
    
    def _binarize_single_file(self, file_path):
        """
        Binarize a single file by reading its raw bytes.
        
        :param file_path: Path to the file to binarize
        :return: Tuple of (binary_string, original_file_size)
        :raises FileNotFoundError: If file cannot be found
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        try:
            # Get original file size
            size = os.path.getsize(file_path)
            
            # Read file as binary data
            with open(file_path, 'rb') as file:
                file_data = file.read()
            
            # Convert each byte to 8-bit binary and join them
            binary_sequence = ''.join(format(byte, '08b') for byte in file_data)
            
            return binary_sequence, size
            
        except PermissionError:
            raise ValueError(f"Permission denied reading file: {file_path}")
        except Exception as e:
            raise ValueError(f"Error reading file '{file_path}': {e}")

    def debinarize(self, data, output_file_path=None):
        """
        Restores the original file from the binary data.

        :param binary_data: A BinaryCode object containing the binary sequence
        :param output_file_path: Path to save the restored file (optional)
        :return: Boolean indicating success of restoration
        """
        try:
            # Validate input
            if not isinstance(data, BinaryCode):
                raise TypeError(f"Expected BinaryCode object, got {type(data).__name__}")
            
            if len(data.data) % 8 != 0:
                raise ValueError("Binary sequence length must be divisible by 8 for file restoration")
            
            # Determine output file path
            if output_file_path is None:
                output_file_path = self._get_default_output_path(data.file_paths)
            
            # Ensure output directory exists
            output_dir = os.path.dirname(output_file_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Convert binary sequence back to bytes
            byte_data = self._binary_to_bytes(data.data)
            
            # Write bytes to file
            with open(output_file_path, 'wb') as file:
                file.write(byte_data)
            
            # Verify the file was created successfully
            if os.path.exists(output_file_path) and os.path.getsize(output_file_path) > 0:
                return True
            else:
                return False
                
        except Exception as e:
            # Log error if needed, but return False for any failure
            print(f"Restoration failed: {e}")
            return False
    
    def _binary_to_bytes(self, binary_sequence):
        """
        Convert binary sequence back to bytes.
        
        :param binary_sequence: Binary string to convert
        :return: Bytes object
        :raises ValueError: If binary cannot be converted to valid bytes
        """
        try:
            byte_array = bytearray()
            for i in range(0, len(binary_sequence), 8):
                byte = binary_sequence[i:i+8]
                if len(byte) == 8:
                    byte_array.append(int(byte, 2))
                else:
                    raise ValueError(f"Incomplete byte at position {i}: {byte}")
            
            return bytes(byte_array)
            
        except ValueError as e:
            raise ValueError(f"Error converting binary to bytes: {e}")
        except Exception as e:
            raise ValueError(f"Unexpected error during binary conversion: {e}")
    
    def _get_default_output_path(self, original_paths):
        """
        Generate default output file path for restoration.
        
        :param original_paths: Original file paths from BinaryCode
        :return: Default output file path
        """
        if original_paths and len(original_paths) > 0:
            # Use original path with _restored suffix
            original_path = original_paths[0]
            base, ext = os.path.splitext(original_path)
            return f"{base}_restored{ext}"
        else:
            # No original path, create default in current directory
            return "./restored_file.bin"
    
    def get_supported_types(self):
        """
        Get information about supported file types.
        
        :return: String describing supported types
        """
        return "Any single file type (universal binary converter)"
    
    def __str__(self):
        """String representation of DefaultBinarize."""
        return "DefaultBinarize(single file universal binary converter)"
    
    def __repr__(self):
        """Developer representation of DefaultBinarize."""
        return "DefaultBinarize()"