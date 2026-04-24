import random
from bitarray import bitarray

from .base import Data
from dnabyte.library import Library

class NucleobaseCode(Data):
    """
    NucleobaseCode with efficient packing for DNA sequences.
    
    Two packing strategies:
    1. Nucleotide packing: 2-bit per nucleotide (A=00, C=01, G=10, T=11)
    2. Library index packing: Store indices to library elements (e.g., 8 bits for 256 oligos)
    
    Data structure:
        - Nested lists preserved at all levels except leaves
        - Leaves packed as bitarrays (DNA mode) or integers (library index mode)
    
    Can be instantiated via:
        - NucleobaseCode(data) - from DNA strings (2-bit packing)
        - NucleobaseCode.from_library_indices(data, library) - from indices
    
    :param data: The data structure (multi-level list).
    :param packed: If True, data is already bit-packed; if False, pack data.
    :param library: Library object for assembly-based encodings.
    :param pack_mode: 'dna' (2-bit) or 'library_indices' (index-based).
    :raises TypeError: If data is not a list.
    :raises ValueError: If data is empty or has invalid structure.
    """

    # 2-bit packing scheme for DNA bases
    DNA_TO_BITS = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    BITS_TO_DNA = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    
    # nucleotide complement mapping
    COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, data, library=None):
        """
        Creates a NucleobaseCode object from data structure.

        :param data: Multi-level list structure representing codewords.
                    Can contain: DNA strings, library indices (int), or bitarrays (if pre-packed).
        :param library: Library object. If provided, assumes data contains library indices.
        """
        self.library = library
        
        # Infer pack_mode from library presence
        self.pack_mode = 'library_indices' if library is not None else 'dna'
        
        # Calculate bits needed for library indices
        if self.pack_mode == 'library_indices':
            import math
            self.index_bits = math.ceil(math.log2(len(library.library)))
        else:
            self.index_bits = None
        
        # Auto-detect if data is already packed
        is_packed = self._is_data_packed(data)
        
        if not is_packed:
            # Validate before packing
            if self.pack_mode == 'dna':
                self._validate_data(data)
            # Pack data according to mode
            self.data = self._pack_structure(data)
        else:
            # Data already packed
            self.data = data
        
        self.num_codewords = len(self.data) if self.data else 0
        self.max_depth = self._calculate_max_depth()

    def _is_data_packed(self, data):
        """
        Check if data is already packed (contains bitarrays at leaves).
        
        :param data: Data structure to check
        :return: True if data contains bitarrays, False if strings/ints
        """
        if isinstance(data, bitarray):
            return True
        elif isinstance(data, (str, int)):
            return False
        elif isinstance(data, list):
            if len(data) == 0:
                return False
            # Check first leaf element
            return self._is_data_packed(data[0])
        else:
            return False

    @classmethod
    def from_library_indices(cls, data, library):
        """
        Alternative constructor: Create NucleobaseCode from library indices.
        
        For assembly-based encodings where data is represented as sequences
        of library element indices rather than raw DNA sequences.
        
        :param data: Nested list structure with integer indices at leaves.
        :param library: Library object containing the oligo library.
        :return: NucleobaseCode object with index-packed data.
        
        Example:
            >>> library = Library(structure='linear_assembly', filename='lib.csv')
            >>> data = [[0, 5, 12], [3, 8, 15]]  # indices into library
            >>> nbc = NucleobaseCode.from_library_indices(data, library)
        """
        if library is None:
            raise ValueError("Library must be provided for from_library_indices")
        
        # Simply pass library - pack_mode will be inferred automatically
        return cls(data, library=library)

    def _pack_dna_string(self, dna_str):
        """
        Convert DNA string to 2-bit packed bitarray.
        
        :param dna_str: DNA string (A, C, G, T)
        :return: bitarray with 2 bits per nucleotide
        """
        if not dna_str:
            return bitarray()
        
        # Build bit string: each base -> 2 bits
        bit_string = ''.join(self.DNA_TO_BITS[base] for base in dna_str)
        return bitarray(bit_string)

    def _pack_library_index(self, index):
        """
        Convert library index to fixed-width bitarray.
        
        :param index: Integer index into library
        :return: bitarray representing the index
        """
        if self.index_bits is None:
            raise ValueError("Cannot pack library index without library")
        
        # Convert index to binary with fixed width
        bit_string = format(index, f'0{self.index_bits}b')
        return bitarray(bit_string)

    def _unpack_dna_bitarray(self, ba):
        """
        Convert 2-bit packed bitarray back to DNA string.
        
        :param ba: bitarray with 2 bits per nucleotide
        :return: DNA string
        """
        if len(ba) == 0:
            return ''
        
        # Read 2 bits at a time
        dna = []
        for i in range(0, len(ba), 2):
            two_bits = ba[i:i+2].to01()
            dna.append(self.BITS_TO_DNA[two_bits])
        return ''.join(dna)

    def _unpack_library_index(self, ba):
        """
        Convert bitarray back to library index and look up sequence.
        
        :param ba: bitarray representing library index
        :return: DNA string from library
        """
        if self.library is None:
            raise ValueError("Cannot unpack library index without library")
        
        # Convert bitarray to integer index
        index = int(ba.to01(), 2)
        
        # Look up in library
        if index >= len(self.library.library):
            raise ValueError(f"Library index {index} out of range (library size: {len(self.library.library)})")
        
        return self.library.library[index]

    def _pack_structure(self, data):
        """
        Recursively pack nested structure.
        
        :param data: Nested list structure with strings/indices at leaves
        :return: Same structure with bitarrays at leaves
        """
        if self.pack_mode == 'dna':
            # DNA mode: strings -> 2-bit bitarrays
            if isinstance(data, str):
                return self._pack_dna_string(data)
            elif isinstance(data, bitarray):
                return data  # Already packed
            elif isinstance(data, list):
                return [self._pack_structure(item) for item in data]
            else:
                return data
        
        elif self.pack_mode == 'library_indices':
            # Library index mode: integers -> fixed-width bitarrays
            if isinstance(data, int):
                return self._pack_library_index(data)
            elif isinstance(data, bitarray):
                return data  # Already packed
            elif isinstance(data, list):
                return [self._pack_structure(item) for item in data]
            else:
                return data
        
        else:
            raise ValueError(f"Unknown pack_mode: {self.pack_mode}")

    def _unpack_structure(self, data):
        """
        Recursively unpack structure.
        
        :param data: Nested list structure with bitarrays at leaves
        :return: Same structure with DNA strings at leaves
        """
        if isinstance(data, bitarray):
            # Unpack based on mode
            if self.pack_mode == 'dna':
                return self._unpack_dna_bitarray(data)
            elif self.pack_mode == 'library_indices':
                return self._unpack_library_index(data)
            else:
                raise ValueError(f"Unknown pack_mode: {self.pack_mode}")
        elif isinstance(data, str):
            return data  # Already unpacked
        elif isinstance(data, list):
            return [self._unpack_structure(item) for item in data]
        else:
            return data

    def to_string_format(self):
        """
        Unpack entire structure back to string format for compatibility.
        
        :return: Nested list structure with DNA strings at leaves
        """
        return self._unpack_structure(self.data)

    @staticmethod
    def complement(string):
        """Get complement of DNA string using class COMPLEMENT dict."""
        return ''.join(NucleobaseCode.COMPLEMENT[base] for base in string)

    @staticmethod
    def reverse_complement(string):
        """Get reverse complement of DNA string."""
        return ''.join(NucleobaseCode.COMPLEMENT[base] for base in reversed(string))

    @staticmethod
    def from_fasta(filename, as_list=False):
        """
        Import sequences from a FASTA file and create a NucleobaseCode object.
        
        :param filename: Path to FASTA file
        :param as_list: If True, return sequences in list structure (one list per entry);
                       If False, return flat list of sequences (default)
        :return: NucleobaseCode object with sequences from FASTA file
        :raises FileNotFoundError: If file doesn't exist
        :raises ValueError: If file is empty or contains invalid sequences
        
        Example FASTA format:
            >sequence_1
            ATCGATCG
            >sequence_2
            GCTAGCTA
        """
        import os
        
        if not os.path.exists(filename):
            raise FileNotFoundError(f"FASTA file not found: {filename}")
        
        sequences = []
        current_seq = []
        current_header = None
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Header line
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_seq:
                        seq = ''.join(current_seq).upper()
                        if as_list:
                            sequences.append([seq])  # Wrap in list
                        else:
                            sequences.append(seq)
                        current_seq = []
                    
                    current_header = line[1:]  # Remove '>'
                
                # Sequence line
                else:
                    # Remove whitespace and convert to uppercase
                    current_seq.append(line.replace(' ', '').upper())
        
        # Don't forget the last sequence
        if current_seq:
            seq = ''.join(current_seq).upper()
            if as_list:
                sequences.append([seq])
            else:
                sequences.append(seq)
        
        if not sequences:
            raise ValueError(f"No sequences found in FASTA file: {filename}")
        
        # Create and return NucleobaseCode object (will validate and pack)
        return NucleobaseCode(sequences)

    def _validate_data(self, data):
        """
        Validates that the input is a proper encoded data structure.

        :param data: The data to validate
        :raises TypeError: If data is not a list
        :raises ValueError: If data is empty or has invalid structure
        """
        if not isinstance(data, list):
            raise TypeError(
                f"Encoded data must be a list, got {type(data).__name__}"
            )

        if len(data) == 0:
            raise ValueError("Encoded data cannot be empty")

        # Validate DNA characters in strings
        self._validate_structure_recursive(data)

    def _validate_structure_recursive(self, structure):
        """
        Recursively validate that DNA strings contain only valid bases.
        
        :param structure: Structure to validate
        :raises ValueError: If invalid DNA characters found
        """
        if isinstance(structure, str):
            # Validate DNA string
            valid_bases = set('ACGT')
            invalid_chars = set(char for char in structure if char not in valid_bases)
            if invalid_chars:
                raise ValueError(
                    f"DNA sequence contains invalid characters: {invalid_chars}. "
                    f"Only A, C, G, T are allowed."
                )
        elif isinstance(structure, list):
            for item in structure:
                self._validate_structure_recursive(item)

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
            return 0  # Leaf (bitarray or string) has no depth

    def get_codeword(self, index, unpack=False):
        """
        Get a specific codeword by index.

        :param index: Index of the codeword
        :param unpack: If True, return as strings; if False, return as bitarrays
        :return: The codeword at the specified index
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self.data):
            raise IndexError(
                f"Codeword index {index} out of range (0-{len(self.data)-1})"
            )
        
        codeword = self.data[index]
        if unpack:
            return self._unpack_structure(codeword)
        return codeword

    def validate(self):
        """
        Re-validates the current data structure.

        :return: True if valid
        :raises: Various exceptions if invalid
        """
        # Unpack and validate
        unpacked = self.to_string_format()
        self._validate_data(unpacked)
        return True

    @staticmethod
    def random(type='synthesis', library=None, n=1000, m=250):
        """
        Generates a random NucleobaseCode object with a specified number of codewords.

        :param type: Type of nucleobase code generation (see nucleobasecode.py for types)
        :param library: Library object or filename string (required for assembly types)
        :param n: Number of codewords to generate
        :param m: Length/number of elements per codeword
        :return: NucleobaseCode object with random data (bitarray-packed)
        """
        data = []
        
        # Helper to ensure we have a Library object
        def get_library(lib):
            import os
            workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
            
            if lib is None:
                possible_paths = [
                    os.path.join(workspace_root, 'tests', 'testlibraries', '20bp_Lib.csv'),
                ]
                for path in possible_paths:
                    if os.path.exists(path):
                        return Library(structure='linear_assembly', filename=path)
                return Library(structure='linear_assembly', filename=possible_paths[0])
            elif isinstance(lib, str):
                if '/' in lib or '\\' in lib:
                    return Library(structure='linear_assembly', filename=lib)
                else:
                    possible_paths = [
                        os.path.join(workspace_root, 'tests', 'testlibraries', lib),
                    ]
                    for path in possible_paths:
                        if os.path.exists(path):
                            return Library(structure='linear_assembly', filename=path)
                    return Library(structure='linear_assembly', filename=possible_paths[0])
            else:
                return lib
        
        def get_poly_library(lib):
            import os
            workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
            
            if lib is None:
                lib = 'polymeraselibinparts.txt'
            
            if isinstance(lib, str):
                if '/' in lib or '\\' in lib:
                    return Library(structure='positional_assembly', filename=lib)
                else:
                    possible_paths = [
                        os.path.join(workspace_root, 'tests', 'testlibraries', lib),
                    ]
                    for path in possible_paths:
                        if os.path.exists(path):
                            return Library(structure='positional_assembly', filename=path)
                    return Library(structure='positional_assembly', filename=possible_paths[0])
            else:
                return lib

        if type == 'synthesis':
            for _ in range(n):
                codeword = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(m))
                data.append(codeword)

        elif type == 'max_density' or type == 'no_homopolymer':
            for _ in range(n):
                codeword = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(m))
                data.append(codeword)

        elif type == 'linear_assembly':
            library = get_library(library)
            for _ in range(n):
                codeword = random.choices(library.library, k=m)
                data.append(codeword)

        elif type == 'linear_chain':
            library = get_library(library)
            for _ in range(n):
                codeword = []
                num_oligos = random.randint(5, m)
                for i in range(num_oligos):
                    codeword.append(random.choice(library.library))
                    if i < num_oligos - 1:
                        connector = ''.join(random.choice(['A', 'T', 'C', 'G']) 
                                          for _ in range(len(library.library[0])))
                        codeword.append(connector)
                data.append(codeword)

        elif type == 'linear_binom':
            library = get_library(library)
            for _ in range(n):
                codeword = []
                num_trios = random.randint(5, m // 3)
                for _ in range(num_trios):
                    trio = [
                        random.choice(library.library),
                        random.choice(library.library),
                        random.choice(library.library)
                    ]
                    codeword.append(trio)
                data.append(codeword)

        elif type == 'poly_chain':
            lib = get_poly_library(library)
            messages = lib.messages
            generic = lib.generic
            positions = lib.position
            comp = NucleobaseCode.complement  # Local reference for efficiency
            
            for _ in range(n):
                codeword = []
                num_triplets = random.randint(3, max(3, m // 5))
                for j in range(num_triplets):
                    msg = random.choice(messages)
                    msg_half = len(msg) // 2
                    
                    if j % 2 == 0:
                        dnatriplet = [
                            positions[j % len(positions)] + generic[0],
                            [
                                comp(msg[:msg_half]) + comp(generic[0]),
                                msg,
                                comp(generic[1]) + comp(msg[msg_half:])
                            ],
                            generic[1] + positions[(j + 1) % len(positions)]
                        ]
                    else:
                        dnatriplet = [
                            comp(generic[0]) + comp(positions[j % len(positions)]),
                            [
                                generic[0] + msg[:msg_half],
                                comp(msg),
                                msg[msg_half:] + generic[1]
                            ],
                            comp(positions[(j + 1) % len(positions)]) + comp(generic[1])
                        ]
                    codeword.append(dnatriplet)
                data.append(codeword)

        elif type == 'poly_binom':
            lib = get_poly_library(library)
            messages = lib.messages
            generic = lib.generic
            positions = lib.position
            comp = NucleobaseCode.complement  # Local reference for efficiency
            
            for _ in range(n):
                codeword = []
                num_groups = random.randint(2, max(2, m // 10))
                for j in range(num_groups):
                    oligo_trios = []
                    num_trios = random.randint(2, 5)
                    for k in range(num_trios):
                        msg = random.choice(messages)
                        msg_half = len(msg) // 2
                        
                        if j % 2 == 0:
                            dnatriple = [
                                positions[j % len(positions)] + generic[0],
                                [
                                    comp(msg[:msg_half]) + comp(generic[0]),
                                    msg,
                                    comp(generic[1]) + comp(msg[msg_half:])
                                ],
                                generic[1] + positions[(j + 1) % len(positions)]
                            ]
                        else:
                            dnatriple = [
                                comp(generic[0]) + comp(positions[j % len(positions)]),
                                [
                                    generic[0] + msg[:msg_half],
                                    comp(msg),
                                    msg[msg_half:] + generic[1]
                                ],
                                comp(positions[(j + 1) % len(positions)]) + comp(generic[1])
                            ]
                        oligo_trios.append(dnatriple)
                    codeword.append(oligo_trios)
                data.append(codeword)

        elif type == 'nested_assembly':
            library = get_library(library)
            
            def create_nested_structure(current_depth, max_depth):
                if current_depth >= max_depth or (current_depth > 0 and random.random() < 0.3):
                    return random.choice(library.library)
                num_children = random.randint(1, 3)
                return [create_nested_structure(current_depth + 1, max_depth) 
                       for _ in range(num_children)]
            
            for _ in range(n):
                max_depth = random.randint(1, 5)
                codeword = create_nested_structure(0, max_depth)
                data.append(codeword)

        else:
            raise ValueError(
                f"Invalid synthesis type: {type}. Must be one of: synthesis, "
                f"linear_assembly, linear_chain, linear_binom, max_density, "
                f"no_homopolymer, poly_chain, poly_binom, nested_assembly"
            )
        
        # Return packed version (will pack strings to bitarrays)
        return NucleobaseCode(data)

    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Number of codewords: {self.num_codewords}\n"
        output += f"Max depth: {self.max_depth}\n"

        if hasattr(self, "file_paths") and self.file_paths:
            output += f"File paths: {self.file_paths}\n"

        if hasattr(self, "size") and self.size:
            output += f"Original size: {self.size} bytes\n"

        # Show unpacked version for readability
        unpacked_preview = str(self.to_string_format())[:100]
        output += f"DATA (unpacked): {unpacked_preview}...\n"
        return output

    def __repr__(self):
        return f"NucleobaseCode(codewords={self.num_codewords}, depth={self.max_depth})"

    def __len__(self):
        """Return the number of codewords."""
        return len(self.data)

    def __getitem__(self, index):
        """Allow indexing to get codewords (returns packed bitarray structure)."""
        return self.get_codeword(index, unpack=False)
