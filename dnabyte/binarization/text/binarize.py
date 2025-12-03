import os
from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.binarize import Binarize

def attributes(params):
    if 'file_paths' not in params.__dict__ or params.file_paths is None:
        raise ValueError("file_paths parameter must be specified")
    else:
        file_paths = params.file_paths

    if 'text_encoding' not in params.__dict__ or params.text_encoding is None:
        text_encoding = 'utf-8'
    else:
        text_encoding = params.text_encoding
        
    return {"file_paths": file_paths, "text_encoding": text_encoding}

class TextBinarize(Binarize):
    """
    TextBinarize converts a single text file to binary representation using UTF-8 encoding.
    
    This binarizer specifically handles one .txt or .text file at a time.
    """
    
    # Supported text file extensions
    SUPPORTED_EXTENSIONS = {'.txt', '.text'}
    
    def __init__(self, params):
        """
        Initialize TextBinarize with specified encoding.
        
        :param encoding: Text encoding to use (default: utf-8)
        """
        self.params = params
        self.text_encoding = params.text_encoding if hasattr(params, 'text_encoding') else 'utf-8'
        self.name = 'text'

    def binarize(self, data):
        """
        Converts a single text file in the Data object to a binary sequence using UTF-8 encoding.

        :param data: A Data object containing exactly one file path
        :return: A BinaryCode object with binary sequence
        :raises TypeError: If data is not a Data object
        :raises ValueError: If data contains more than one file, file is not a text file, or cannot be processed
        :raises FileNotFoundError: If the file cannot be found
        :raises UnicodeDecodeError: If file cannot be decoded with specified encoding
        """
        # Validate input
        if not isinstance(data, Data):
            raise TypeError(f"Expected Data object, got {type(data).__name__}")
        
        # Check that there's exactly one file
        if len(data.file_paths) != 1:
            raise ValueError(f"TextBinarize can only handle exactly one file, got {len(data.file_paths)} files. If you want to binarize multiple files, use an archive format (e.g. CompressedBinarize) instead.")
        
        file_path = data.file_paths[0]
        
        # Validate the file is a text file
        self._validate_text_file(file_path)
        
        try:
            # Binarize the single file
            binary_sequence, size = self._binarize_single_file(file_path)
            
            # Create and return BinaryCode object
            return binary_sequence, size
            
            
        except Exception as e:
            raise ValueError(f"Error during binarization: {e}")
    
    def _validate_text_file(self, file_path):
        """
        Validate that the file is a text file with supported extension.
        
        :param file_path: File path to validate
        :raises ValueError: If file is not a supported text file
        """
        # Check file extension
        _, ext = os.path.splitext(file_path.lower())
        if ext not in self.SUPPORTED_EXTENSIONS:
            raise ValueError(
                f"File '{file_path}' has unsupported extension '{ext}'. "
                f"Supported extensions: {', '.join(sorted(self.SUPPORTED_EXTENSIONS))}"
            )
        
    #     # Additional check: try to detect if it's actually a text file
    #     if not self._is_text_file(file_path):
    #         raise ValueError(f"File '{file_path}' does not appear to be a valid text file")
    
    # def _is_text_file(self, file_path):
    #     """
    #     Check if file appears to be a text file by reading a sample.
        
    #     :param file_path: Path to file to check
    #     :return: True if file appears to be text, False otherwise
    #     """
    #     try:
    #         with open(file_path, 'r', encoding=self.encoding) as file:
    #             # Try to read first 1024 bytes as text
    #             sample = file.read(1024)
    #             return True
    #     except (UnicodeDecodeError, UnicodeError):
    #         return False
    #     except Exception:
    #         # If we can't read it, assume it's not a text file
    #         return False
    
    def _binarize_single_file(self, file_path):
        """
        Binarize a single text file.
        
        :param file_path: Path to the file to binarize
        :return: Tuple of (binary_string, original_file_size)
        :raises FileNotFoundError: If file cannot be found
        :raises UnicodeDecodeError: If file cannot be decoded
        """
        try:
            # Get original file size
            size = os.path.getsize(file_path)
            
            # Read and convert to binary
            with open(file_path, 'r', encoding=self.text_encoding) as file:
                text_data = file.read()
                
            # Convert each character to 8-bit binary
            binary_sequence = ''.join(format(ord(char), '08b') for char in text_data)

            return binary_sequence, size

        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        except UnicodeDecodeError as e:
            raise UnicodeDecodeError(
                f"Cannot decode file '{file_path}' with encoding '{self.text_encoding}': {e}"
            )
        except Exception as e:
            raise ValueError(f"Error reading file '{file_path}': {e}")

    def debinarize(self, binary_data, output_file_path=None):
        """
        Restores the original text file from the binary data.

        :param binary_data: A BinaryCode object containing the binary sequence
        :param output_file_path: Path where to save the restored file (optional)
        :return: Boolean indicating success of restoration
        :raises TypeError: If binary_data is not a BinaryCode object
        """
        try:
            # Validate input
            if not isinstance(binary_data, BinaryCode):
                raise TypeError(f"Expected BinaryCode object, got {type(binary_data).__name__}")
            
            if len(binary_data.data) % 8 != 0:
                raise ValueError("Binary sequence length must be divisible by 8 for text restoration")
            
            # Convert binary back to text
            restored_text = self._binary_to_text(binary_data.data)
            
            # Determine output file path
            if output_file_path is None:
                output_file_path = self._get_default_output_path(binary_data.file_paths)
            
            # Ensure output directory exists
            output_dir = os.path.dirname(output_file_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Write restored text to file
            with open(output_file_path, 'w', encoding=self.text_encoding) as file:
                file.write(restored_text)
            
            # Verify the file was created successfully
            if os.path.exists(output_file_path) and os.path.getsize(output_file_path) > 0:
                return True
            else:
                return False
                
        except Exception as e:
            # Log error if needed, but return False for any failure
            print(f"Restoration failed: {e}")
            return False
    
    def _binary_to_text(self, binary_sequence):
        """
        Convert binary sequence back to text.
        
        :param binary_sequence: Binary string to convert
        :return: Restored text string
        :raises ValueError: If binary cannot be converted to valid text
        """
        try:
            chars = []
            for i in range(0, len(binary_sequence), 8):
                byte = binary_sequence[i:i+8]
                if len(byte) == 8:  # Ensure complete byte
                    char_code = int(byte, 2)
                    chars.append(chr(char_code))
                else:
                    raise ValueError(f"Incomplete byte at position {i}: {byte}")
            
            return ''.join(chars)
            
        except ValueError as e:
            if "chr()" in str(e):
                raise ValueError(f"Invalid character code in binary data: {e}")
            raise
        except Exception as e:
            raise ValueError(f"Error converting binary to text: {e}")
    
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
            return "./restored.txt"
    
    def get_supported_extensions(self):
        """
        Get list of supported file extensions.
        
        :return: Set of supported extensions
        """
        return self.SUPPORTED_EXTENSIONS.copy()
    
    def __str__(self):
        """String representation of TextBinarize."""
        return f"TextBinarize(encoding='{self.text_encoding}', extensions={sorted(self.SUPPORTED_EXTENSIONS)})"
    
    def __repr__(self):
        """Developer representation of TextBinarize."""
        return f"TextBinarize(encoding={self.text_encoding!r})"