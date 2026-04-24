import random
from .base import Data
from dnabyte.library import Library

def complementmap(string):
    """Helper function to get reverse complement of DNA string."""
    compliment = ''
    for i in string:
        if i == 'A':
            compliment += 'T'
        elif i == 'T':
            compliment += 'A'
        elif i == 'C':
            compliment += 'G'
        elif i == 'G':
            compliment += 'C'
    return compliment[::-1]

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

    def random(type='synthesis', library=None, n=1000, m=250):
        """
        Generates a random NucleobaseCode object with a specified number of codewords.

        :param type: Type of nucleobase code generation:
            - 'synthesis': Simple flat DNA strings
            - 'linear_assembly': Flat list of library elements
            - 'linear_chain': List of [oligo, connector, oligo, ...] pattern
            - 'linear_binom': List of [[trio1], [trio2], ...] nested structure
            - 'max_density': Simple DNA string per codeword
            - 'no_homopolymer': Simple DNA string per codeword
            - 'poly_chain': List of nested triplets [[connector, [oligo1, connector, oligo2], connector], ...]
            - 'poly_binom': List of list of triplets [[[connector, [o1, c, o2], connector]], ...]
            - 'nested_assembly': Random depth nested structure
        :param library: Library object or filename string (required for assembly types)
        :param n: Number of codewords to generate
        :param m: Length/number of elements per codeword
        :return: NucleobaseCode object with random data
        """
        data = []
        
        # Helper to ensure we have a Library object
        def get_library(lib):
            import os
            # Get workspace root (3 levels up from this file: data_classes -> dnabyte -> DNAbyte -> workspace)
            workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
            
            if lib is None:
                # Try multiple default locations with absolute paths
                possible_paths = [
                    os.path.join(workspace_root, 'app', 'static', 'libraries', '20bp_Lib.csv'),
                    os.path.join(workspace_root, 'tests', 'testlibraries', '20bp_Lib.csv'),
                    os.path.join(workspace_root, 'DNAbyte', 'tests', 'testlibraries', '20bp_Lib.csv')
                ]
                for path in possible_paths:
                    if os.path.exists(path):
                        return Library(structure='linear_assembly', filename=path)
                # Fallback to first path if none exist (will error downstream)
                return Library(structure='linear_assembly', filename=possible_paths[0])
            elif isinstance(lib, str):
                # If it's a string, treat it as a filename
                if '/' in lib or '\\' in lib:
                    # Full path provided
                    return Library(structure='linear_assembly', filename=lib)
                else:
                    # Just filename - try multiple locations with absolute paths
                    possible_paths = [
                        os.path.join(workspace_root, 'app', 'static', 'libraries', lib),
                        os.path.join(workspace_root, 'tests', 'testlibraries', lib),
                        os.path.join(workspace_root, 'DNAbyte', 'tests', 'testlibraries', lib)
                    ]
                    for path in possible_paths:
                        if os.path.exists(path):
                            return Library(structure='linear_assembly', filename=path)
                    # If none found, use first path (will error if file doesn't exist)
                    return Library(structure='linear_assembly', filename=possible_paths[0])
            else:
                # Assume it's already a Library object
                return lib
        
        # Helper for poly libraries (positional_assembly structure)
        def get_poly_library(lib):
            import os
            workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
            
            if lib is None:
                lib = 'polymeraselibinparts.txt'
            
            if isinstance(lib, str):
                if '/' in lib or '\\' in lib:
                    return Library(structure='positional_assembly', filename=lib)
                else:
                    possible_paths = [
                        os.path.join(workspace_root, 'app', 'static', 'libraries', lib),
                        os.path.join(workspace_root, 'tests', 'testlibraries', lib),
                        os.path.join(workspace_root, 'DNAbyte', 'tests', 'testlibraries', lib)
                    ]
                    for path in possible_paths:
                        if os.path.exists(path):
                            return Library(structure='positional_assembly', filename=path)
                    return Library(structure='positional_assembly', filename=possible_paths[0])
            else:
                return lib

        if type == 'synthesis':
            # Simple flat DNA strings
            for _ in range(n):
                codeword = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(m))
                data.append(codeword)

        elif type == 'max_density' or type == 'no_homopolymer':
            # Simple DNA string per codeword (like synthesis)
            for _ in range(n):
                codeword = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(m))
                data.append(codeword)

        elif type == 'linear_assembly':
            library = get_library(library)
            
            # Generate n codewords, each with m library elements
            for _ in range(n):
                codeword = random.choices(library.library, k=m)
                data.append(codeword)

        elif type == 'linear_chain':
            # Pattern: [oligo, connector, oligo, connector, ..., oligo]
            library = get_library(library)
            
            for _ in range(n):
                codeword = []
                num_oligos = random.randint(5, m)
                for i in range(num_oligos):
                    # Add oligo
                    codeword.append(random.choice(library.library))
                    # Add connector (except after last oligo)
                    if i < num_oligos - 1:
                        connector = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(len(library.library[0])))
                        codeword.append(connector)
                data.append(codeword)

        elif type == 'linear_binom':
            # Pattern: [[trio1], [trio2], ...] where each trio is a simple structure
            library = get_library(library)
            
            for _ in range(n):
                codeword = []
                num_trios = random.randint(5, m // 3)
                for _ in range(num_trios):
                    # Create a simple trio structure
                    trio = [
                        random.choice(library.library),
                        random.choice(library.library),
                        random.choice(library.library)
                    ]
                    codeword.append(trio)
                data.append(codeword)

        elif type == 'poly_chain':
            # Pattern: [[connector, [oligo1, connector_inner, oligo2], connector], ...]
            lib = get_poly_library(library)
            
            messages = lib.messages
            generic = lib.generic
            positions = lib.position
            
            for _ in range(n):
                codeword = []
                num_triplets = random.randint(3, max(3, m // 5))
                for j in range(num_triplets):
                    # Select a random message
                    msg = random.choice(messages)
                    msg_half = len(msg) // 2
                    
                    # Create triplet structure matching encode output
                    # Alternating between even (5'->3') and odd (3'->5') encoding
                    if j % 2 == 0:  # Even: info on 5' to 3' strand
                        dnatriplet = [
                            positions[j % len(positions)] + generic[0],
                            [
                                complementmap(msg[:msg_half]) + complementmap(generic[0]),
                                msg,
                                complementmap(generic[1]) + complementmap(msg[msg_half:])
                            ],
                            generic[1] + positions[(j + 1) % len(positions)]
                        ]
                    else:  # Odd: info on 3' to 5' strand
                        dnatriplet = [
                            complementmap(generic[0]) + complementmap(positions[j % len(positions)]),
                            [
                                generic[0] + msg[:msg_half],
                                complementmap(msg),
                                msg[msg_half:] + generic[1]
                            ],
                            complementmap(positions[(j + 1) % len(positions)]) + complementmap(generic[1])
                        ]
                    
                    codeword.append(dnatriplet)
                data.append(codeword)

        elif type == 'poly_binom':
            # Pattern: [[[triplet1, triplet2, ...]], [[group2...]], ...]
            lib = get_poly_library(library)
            
            messages = lib.messages
            generic = lib.generic
            positions = lib.position
            
            for _ in range(n):
                codeword = []
                num_groups = random.randint(2, max(2, m // 10))
                for j in range(num_groups):
                    oligo_trios = []
                    # Each group has multiple oligos (sigma amount)
                    num_trios = random.randint(2, 5)
                    for k in range(num_trios):
                        # Select a random message
                        msg = random.choice(messages)
                        msg_half = len(msg) // 2
                        
                        # Create triplet structure matching encode output
                        # Alternating based on position j
                        if j % 2 == 0:  # Even: info on 5' to 3' strand
                            dnatriple = [
                                positions[j % len(positions)] + generic[0],
                                [
                                    complementmap(msg[:msg_half]) + complementmap(generic[0]),
                                    msg,
                                    complementmap(generic[1]) + complementmap(msg[msg_half:])
                                ],
                                generic[1] + positions[(j + 1) % len(positions)]
                            ]
                        else:  # Odd: info on 3' to 5' strand
                            dnatriple = [
                                complementmap(generic[0]) + complementmap(positions[j % len(positions)]),
                                [
                                    generic[0] + msg[:msg_half],
                                    complementmap(msg),
                                    msg[msg_half:] + generic[1]
                                ],
                                complementmap(positions[(j + 1) % len(positions)]) + complementmap(generic[1])
                            ]
                        
                        oligo_trios.append(dnatriple)
                    codeword.append(oligo_trios)
                data.append(codeword)

        elif type == 'nested_assembly':
            library = get_library(library)
            
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

        else:
            raise ValueError(f"Invalid synthesis type: {type}. Must be one of: synthesis, linear_assembly, linear_chain, linear_binom, max_density, no_homopolymer, poly_chain, poly_binom, nested_assembly")
        
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
