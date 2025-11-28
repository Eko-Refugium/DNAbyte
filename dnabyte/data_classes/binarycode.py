from .base import Data

class BinaryCode(Data):
    """
    Represents binary data as a sequence of '0' and '1' characters.
    
    BinaryCode stores and manipulates binary data streams created from file 
    binarization processes. It provides validation, indexing, and iteration 
    capabilities for working with binary representations of original data.
    
    Extends the Data base class to maintain compatibility with the DNAbyte 
    pipeline while providing binary-specific functionality.

    Can be created either:
        - from file data using a Binarize method, or
        - directly by providing a binary string.

    :param bitstream: A string of binary data (only '0' and '1' characters).
    :param file_paths: Optional list of original file paths (for inheritance compatibility).
    :param size: Optional original file size in bytes.
    :raises TypeError: If bitstream is not a string.
    :raises ValueError: If bitstream is empty or contains non-binary characters.
    
    Example:
        >>> binary = BinaryCode("10110010", file_paths=["/path/to/file.txt"])
        >>> print(len(binary))  # 8
        >>> print(binary[0])    # '1'
    """

    def __init__(self, data, file_paths=None, size=None):
        """
        Creates a BinaryCode object from a bitstream.

        :param bitstream: A string of binary data.
        :param file_paths: Optional list of original file paths.
        :param size: Optional original file size in bytes.
        """
        
        #print(data, 'data in binarycode init')
        # Validate the bitstream before setting any attributes
        self._validate_bitstream(data)

        # Set binary-specific attributes
        self.data = data
        self.length = len(data)

        # Handle inheritance - set attributes for parent class compatibility
        self.file_paths = file_paths or []
        self.size = size  # Keep track of original file size

    def _validate_bitstream(self, data):
        """
        Validates that the input is a proper bitstream.

        :param bitstream: The data to validate
        :raises TypeError: If bitstream is not a string
        :raises ValueError: If bitstream is empty or contains invalid characters
        """
        # Check if data is a string
        if not isinstance(data, str):
            raise TypeError(
                f"Bitstream must be a string, got {type(data).__name__}"
            )

        # Check if data is empty
        if len(data) == 0:
            raise ValueError("Bitstream cannot be empty")

        # Check if all characters are either '0' or '1'
        invalid_chars = set(char for char in data if char not in "01")
        if invalid_chars:
            raise ValueError(
                f"Bitstream contains invalid characters: {invalid_chars}. "
                f"Only '0' and '1' are allowed."
            )

    def compare_binary_strings(self, bin_str1, bin_str2):
        """
        Helper function for compare(). Compare two binary strings and return the indices of the differing bits.
        """
        differences = []
        min_length = min(len(bin_str1), len(bin_str2))
        
        for i in range(min_length):
            if bin_str1[i] != bin_str2[i]:
                differences.append(i)
        
        return differences


    def compare(self, data_dec, data_bin, logger=None):
        """
        Compare decoded data with raw data and return the result of the comparison.

        Args:
            data_dec (BinaryCode): The decoded data object.
            data_bin (BinaryCode): The raw data object.

        Returns:
            str: A message indicating whether the data was successfully decoded or not.
        """
        res = self.compare_binary_strings(data_dec.data, data_bin.data)

        if res == [] and len(data_dec.data) == len(data_bin.data):
            return 'SUCCESS', res
        else:
            return 'ERROR', res

    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Length: {self.length} bits\n"
        if self.size is not None:
            output += f"Size: {self.size} bytes\n"
        output += f"Data: {self.data[0:50]}...\n"
        return output

    def __repr__(self):
        """
        Developer-friendly representation of the BinaryCode object.

        :return: Unambiguous string representation
        """
        # Truncate data if too long for readable repr
        if len(self.data) > 20:
            data_repr = f"{self.data[:20]}..."
        else:
            data_repr = self.data

        return f"BinaryCode(bitstream='{data_repr}', length={self.length})"

    def __len__(self):
        """Return the number of bits in the bitstream."""
        return len(self.data)

    def __getitem__(self, index):
        """Allow indexing to get individual bits or slices."""
        return self.data[index]

    def __iter__(self):
        """Allow iteration over individual bits."""
        return iter(self.data)
