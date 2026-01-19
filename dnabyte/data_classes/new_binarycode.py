from bitarray import bitarray
from bitarray.util import urandom

from .base import Data

class BinaryCode(Data):
    """
    Represents binary data as a bitarray for efficient storage and manipulation.

    BinaryCode stores and manipulates binary data streams created from file 
    binarization processes. It provides validation, indexing, and iteration 
    capabilities for working with binary representations of original data.

    Extends the Data base class to maintain compatibility with the DNAbyte 
    pipeline while providing binary-specific functionality.
    """
    def __init__(self, data, file_paths=None, size=None):

        # Validate the bitstream before setting any attributes
        self._validate_bitstream(data)

        if isinstance(data, str):
            # Convert bit string to bytes
            self.data = bitarray(data)
            self.length = len(data)  # number of bits
            self.size = (self.length + 7) // 8  # size in bytes

        elif isinstance(data, bitarray):
            self.data = data
            self.length = len(data)
            self.size = (self.length + 7) // 8  # size in bytes

        else:
            raise TypeError(
                f"Bitstream must be a string or bitarray, got {type(data).__name__}"
            )
        
        self.file_paths = file_paths or []
    
    def to_bitstring(self):
        """Convert bytes back to '01010...' string for compatibility."""
        return self.data.to01()
    
    def _validate_bitstream(self, data):
        """
        Validates that the input is a proper bitstream.

        :param bitstream: The data to validate
        :raises TypeError: If bitstream is not a string
        :raises ValueError: If bitstream is empty or contains invalid characters
        """

        if isinstance(data, str):
            # Check if all characters are either '0' or '1'
            invalid_chars = set(char for char in data if char not in "01")
            if invalid_chars:
                raise ValueError(
                    f"Bitstream contains invalid characters: {invalid_chars}. "
                    f"Only '0' and '1' are allowed."
                )

        # Check if data is empty
        if len(data) == 0:
            raise ValueError("Bitstream cannot be empty")

    def compare_binary_strings(self, ba1, ba2):
        """
        Compare two bitarrays and return the indices of the differing bits.
        
        Uses XOR operation for efficient comparison of bitarrays.
        
        :param ba1: First bitarray to compare
        :param ba2: Second bitarray to compare
        :return: List of indices where bits differ
        """
        differences = []
        min_length = min(len(ba1), len(ba2))
        
        # XOR the two bitarrays to find differences (1 = different, 0 = same)
        # Only compare up to the minimum length
        xor_result = ba1[:min_length] ^ ba2[:min_length]
        
        # Find all indices where the bit is 1 (i.e., different)
        for i in range(len(xor_result)):
            if xor_result[i]:
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
        elif len(data_dec.data) <= len(data_bin.data):
            return 'ERROR_short', res
        else:
            return 'ERROR_long', res

    def random(n):
        """
        Generates a random BinaryCode object with a specified length.

        :param n: Length of the binary string to generate.
        :return: BinaryCode object with random binary data.
        """
        bitstream = urandom(n)
        return BinaryCode(bitstream)
    
    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Length: {self.length} bits\n"
        if self.size is not None:
            output += f"Size: {self.size} bytes\n"
        # Convert bitarray to string for display
        preview = self.data[:50].to01() if len(self.data) > 50 else self.data.to01()
        output += f"Data: {preview}{'...' if len(self.data) > 50 else ''}\n"
        return output

    def __repr__(self):
        """
        Developer-friendly representation of the BinaryCode object.

        :return: Unambiguous string representation
        """
        # Truncate data if too long for readable repr
        if len(self.data) > 20:
            data_repr = f"{self.data[:20].to01()}..."
        else:
            data_repr = self.data.to01()

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