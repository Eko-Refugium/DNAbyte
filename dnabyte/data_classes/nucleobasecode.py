import random
from .base import Data
from dnabyte.library import Library

class NucleobaseCode(Data):
    """
    NucleobaseCode is a subclass of Data. Its main attribute is the data object, which consists of a multi-leveled list. At the highest level, data consists of a list of codewords, there is one list for every codeword in the data object.
    At the lowest level there are lists of motifs or nucleotides, that are to be hybridised in a single pool. The list
    structure represents the hierarchical structure of the pooling process.

    It can be instantiated either by:
        - by an encoder producing the encoded data structure or
        - directly providing the encoded data structure as a parameter.

    :param data: The encoded data structure (multi-level list of codewords).
    :param file_paths: Optional list of original file paths.
    :param size: Optional size information from original data.
    :raises TypeError: If data is not a list.
    :raises ValueError: If data is empty or has invalid structure.
    """

    def __init__(self, data):
        """
        Creates a NucleobaseCode object from encoded data structure.

        :param data: Multi-level list structure representing encoded codewords.
        """
        # Validate the encoded data before setting any attributes
        self._validate_data(data)
        self.data = data
        self.num_codewords = len(data) if data else 0

        # Calculate additional metrics
        self.max_depth = self._calculate_max_depth()

    def _validate_data(self, data):
        """
        Validates that the input is a proper encoded data structure.

        :param data: The data to validate
        :raises TypeError: If data is not a list
        :raises ValueError: If data is empty or has invalid structure
        """
        # Check if data is a list
        if not isinstance(data, list):
            raise TypeError(
                f"Encoded data must be a list, got {type(data).__name__}"
            )

        # Check if data is empty
        if len(data) == 0:
            raise ValueError("Encoded data cannot be empty")

        # TODO: Check structure of the codewords
        # TODO: Check whether each codeword contains only valid nucleotides
        # self._validate_codeword_structure(data)


    # TODO: calculate_max_depth delivers incorrect results in some tests
    def _calculate_max_depth(self):
        """
        Calculates the maximum depth of the nested structure.

        :return: Maximum nesting depth
        """
        if not self.data:
            return 0

        max_depth = 0
        for codeword in self.data:
            depth = self._get_depth_recursive(codeword)
            max_depth = max(max_depth, depth)
        return max_depth

    def _get_depth_recursive(self, structure):
        """
        Recursively calculates depth of nested structure.

        :param structure: The structure to measure depth of
        :return: Depth of structure
        """
        if isinstance(structure, list):
            if not structure:
                return 1
            return 1 + max(self._get_depth_recursive(item) for item in structure)
        else:
            return 0  # Nucleotide has no depth

    def get_codeword(self, index):
        """
        Get a specific codeword by index.

        :param index: Index of the codeword
        :return: The codeword at the specified index
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self.data):
            raise IndexError(
                f"Codeword index {index} out of range (0-{len(self.data)-1})"
            )
        return self.data[index]

    def validate(self):
        """
        Re-validates the current encoded data structure.

        :return: True if valid
        :raises: Various exceptions if invalid
        """
        self._validate_data(self.data)
        return True

    def random(type='synthesis', library = None, n=1000, m=250):
        """
        Generates a random NucleobaseCode object with a specified number of codewords.

        :param type: Type of nucleobase code (e.g., 'DNA', 'RNA')
        :param n: Number of codewords to generate
        :return: NucleobaseCode object with random data
        """
        data = []

        if type == 'synthesis':

            for _ in range(n):
                codeword = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(m))
                data.append(codeword)

        elif type == 'nested_assembly':
            if library is None:
                raise ValueError("Library must be provided for assembly type")
            
            else:
                library = Library(structure='linear_assembly', filename='./tests/testlibraries/20bp_Lib.csv')
            
                # Helper function to create nested list of random depth
                def create_nested_structure(current_depth, max_depth):
                    """Recursively create nested list structure with random depth."""
                    # If we've reached max depth or randomly decide to stop, return a leaf element
                    if current_depth >= max_depth or (current_depth > 0 and random.random() < 0.3):
                        return random.choice(library.library)
                    
                    # Create a list with 1-3 children
                    num_children = random.randint(1, 3)
                    return [create_nested_structure(current_depth + 1, max_depth) for _ in range(num_children)]
                
                # Generate n codewords, each with random nested structure (depth 1-5)
                for _ in range(n):
                    max_depth = random.randint(1, 5)
                    codeword = create_nested_structure(0, max_depth)
                    data.append(codeword)


        elif type == 'linear_assembly':
            if library is None:
                raise ValueError("Library must be provided for assembly type")
            
            else:
                library = Library(structure='linear_assembly', filename='./tests/testlibraries/20bp_Lib.csv')
            
                # Generate n codewords, each with m library elements
                for _ in range(n):
                    codeword = random.choices(library.library, k=m)
                    data.append(codeword)

        else:
            raise ValueError(f"Invalid synthesis type: {type}")
        
        return NucleobaseCode(data)


    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Number of codewords: {self.num_codewords}\n"
        output += f"Max depth: {self.max_depth}\n"

        # Add file paths if they exist
        if hasattr(self, "file_paths") and self.file_paths:
            output += f"File paths: {self.file_paths}\n"

        if hasattr(self, "size") and self.size:
            output += f"Original size: {self.size} bytes\n"

        output += f"DATA: {str(self.data)[:100]}...\n"
        return output

    def __len__(self):
        """Return the number of codewords."""
        return len(self.data)

    def __getitem__(self, index):
        """Allow indexing to get codewords."""
        return self.get_codeword(index)
