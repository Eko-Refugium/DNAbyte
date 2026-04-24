import random
from bitarray import bitarray

from .base import Data
from dnabyte.library import Library

class NucleobaseCode(Data):
    """
    NucleobaseCode with efficient single-bitarray storage.
    
    All DNA sequences are stored in a single bitarray (eliminating per-object overhead).
    The nested list structure is preserved via metadata.
    
    Two packing strategies:
    1. Nucleotide packing: 2-bit per nucleotide (A=00, C=01, G=10, T=11)
    2. Library index packing: Store indices to library elements (e.g., 8 bits for 256 oligos)
    
    Storage format:
        - _packed_data: Single bitarray containing all sequences concatenated
        - _structure: List describing nested structure for reconstruction
    
    Can be instantiated via:
        - NucleobaseCode(data) - from DNA strings (2-bit packing)
        - NucleobaseCode.from_library_indices(data, library) - from indices
    
    :param data: The data structure (multi-level list).
    :param library: Library object for assembly-based encodings.
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
                    Can contain: DNA strings, library indices (int), or packed format.
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
        
        # Check if data is already in packed format (tuple of bitarray and structure)
        if isinstance(data, tuple) and len(data) == 2:
            self._packed_data, self._structure = data
            if not isinstance(self._packed_data, bitarray):
                raise TypeError("Packed data must be a bitarray")
        else:
            # Validate and pack
            if self.pack_mode == 'dna':
                self._validate_data(data)

            self._packed_data, self._structure = self._flatten_and_pack(data)

        self.num_codewords = self._count_top_level()
        self.max_depth = self._calculate_max_depth()

    @property
    def data(self):
        """
        Property that returns packed format for backward compatibility.
        """
        return (self._packed_data, self._structure)

    def _count_top_level(self):
        """Count number of top-level elements (direct children of root list)."""
        if not self._structure:
            return 0
        
        # The first '[' is the root wrapper
        # Count '[' tokens at depth 1 (children of root)
        depth = 0
        count = 0
        
        for item in self._structure:
            if item == '[':
                depth += 1
                if depth == 2:  # Children of root are at depth 2
                    count += 1
            elif item == ']':
                depth -= 1
            elif isinstance(item, int):
                # It's a leaf at current depth
                if depth == 1:  # Direct child of root
                    count += 1
        
        return count

    def _flatten_and_pack(self, data):
        """
        Flatten nested structure and pack all sequences into single bitarray.
        
        Returns:
            packed_data: Single bitarray with all sequences
            structure: List describing the nested structure
                      Format: ['[', '[', length, ']', length, ']', ...]
                      '[' = begin nested list
                      ']' = end nested list  
                      integer = sequence length (in nucleotides or indices)
        """
        all_bits = []
        structure = []
        
        def flatten_recursive(item):
            if isinstance(item, list):
                structure.append('[')
                for sub_item in item:
                    flatten_recursive(sub_item)
                structure.append(']')
            elif isinstance(item, str):
                # DNA string
                bits = ''.join(self.DNA_TO_BITS[base] for base in item)
                all_bits.append(bits)
                structure.append(len(item))  # Store length in nucleotides
            elif isinstance(item, int):
                # Library index
                bits = format(item, f'0{self.index_bits}b')
                all_bits.append(bits)
                structure.append(1)  # One index = length 1
            else:
                raise TypeError(f"Unexpected type in data: {type(item)}")
        
        flatten_recursive(data)
        
        # Combine all bits
        combined = ''.join(all_bits)
        packed_data = bitarray(combined) if combined else bitarray()
        
        return packed_data, structure

    def _reconstruct_structure(self):
        """
        Reconstruct nested list structure from flat bitarray and structure metadata.
        
        Returns:
            Nested list with DNA strings at leaves
        """
        if not self._structure:
            return []
        
        bit_position = 0
        structure_idx = 0
        
        def reconstruct_recursive():
            nonlocal bit_position, structure_idx
            
            token = self._structure[structure_idx]
            
            if token == '[':
                structure_idx += 1
                result = []
                while self._structure[structure_idx] != ']':
                    result.append(reconstruct_recursive())
                structure_idx += 1  # Skip ]
                return result
            
            elif isinstance(token, int):
                # It's a sequence length
                length = token
                structure_idx += 1
                
                if self.pack_mode == 'dna':
                    # Extract bits for this sequence (2 bits per nucleotide)
                    num_bits = length * 2
                    sequence_bits = self._packed_data[bit_position:bit_position + num_bits]
                    bit_position += num_bits
                    
                    # Decode to DNA string
                    bits_str = sequence_bits.to01()
                    dna = ''.join(self.BITS_TO_DNA[bits_str[i:i+2]] for i in range(0, len(bits_str), 2))
                    return dna
                
                elif self.pack_mode == 'library_indices':
                    # Extract bits for this index
                    num_bits = self.index_bits
                    index_bits = self._packed_data[bit_position:bit_position + num_bits]
                    bit_position += num_bits
                    
                    # Decode to library index and lookup
                    index = int(index_bits.to01(), 2)
                    return self.library.library[index]
            
            else:
                raise ValueError(f"Unexpected token in structure: {token}")
        
        return reconstruct_recursive()

    @classmethod
    def from_library_indices(cls, data, library):
        """
        Alternative constructor: Create NucleobaseCode from library indices.
        
        For assembly-based encodings where data is represented as sequences
        of library element indices rather than raw DNA sequences.
        
        :param data: Nested list structure with integer indices at leaves.
        :param library: Library object containing the oligo library.
        :return: NucleobaseCode object with index-packed data.
        """
        if library is None:
            raise ValueError("Library must be provided for from_library_indices")
        
        return cls(data, library=library)

    def unpack(self):
        """
        Unpack the data structure to nested lists with DNA strings.
        
        :return: Nested list structure with DNA strings at leaves
        """
        return self._reconstruct_structure()

    def _validate_data(self, data):
        """
        Validates the nucleobase data structure.
        
        :param data: The data structure to validate
        :raises TypeError: If structure is invalid
        :raises ValueError: If DNA sequences contain invalid characters
        """
        valid_bases = {'A', 'C', 'G', 'T'}
        
        def validate_recursive(item, depth=0):
            if isinstance(item, list):
                if len(item) == 0:
                    raise ValueError(f"Empty list at depth {depth}")
                for sub_item in item:
                    validate_recursive(sub_item, depth + 1)
            elif isinstance(item, str):
                if len(item) == 0:
                    raise ValueError(f"Empty DNA string at depth {depth}")
                invalid = set(item) - valid_bases
                if invalid:
                    raise ValueError(f"Invalid nucleotides {invalid} in sequence: {item[:50]}")
            else:
                raise TypeError(f"Invalid type {type(item)} at depth {depth}. Expected list or str.")
        
        if not isinstance(data, list):
            raise TypeError(f"Data must be a list, got {type(data)}")
        
        validate_recursive(data)

    def _calculate_max_depth(self):
        """Calculate maximum nesting depth."""
        depth = 0
        current_depth = 0
        for token in self._structure:
            if token == '[':
                current_depth += 1
                depth = max(depth, current_depth)
            elif token == ']':
                current_depth -= 1
        return depth

    @staticmethod
    def complement(sequence):
        """
        Get the complement of a DNA sequence.
        
        :param sequence: DNA sequence string
        :return: Complement sequence string
        """
        return ''.join(NucleobaseCode.COMPLEMENT.get(base, base) for base in sequence)
    
    @staticmethod
    def reverse_complement(sequence):
        """
        Get the reverse complement of a DNA sequence.
        
        :param sequence: DNA sequence string
        :return: Reverse complement sequence string
        """
        return NucleobaseCode.complement(sequence)[::-1]

    @classmethod
    def from_fasta(cls, filename, as_list=True):
        """
        Import DNA sequences from a FASTA file.
        
        :param filename: Path to FASTA file
        :param as_list: If True, return as list structure (required for NucleobaseCode)
        :return: NucleobaseCode object
        """
        sequences = []
        current_seq = []
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # Header line - save previous sequence if exists
                    if current_seq:
                        sequences.append(''.join(current_seq).upper())
                        current_seq = []
                else:
                    # Sequence line
                    current_seq.append(line)
            
            # Don't forget last sequence
            if current_seq:
                sequences.append(''.join(current_seq).upper())
        
        if not sequences:
            raise ValueError(f"No sequences found in {filename}")
        
        # Wrap in list structure
        data = [sequences] if as_list else sequences
        return cls(data)

    @staticmethod
    def random(structure_shape, min_length=10, max_length=50):
        """
        Generate random NucleobaseCode with specified structure.
        
        :param structure_shape: Shape of nested structure, e.g., [3, [2, 2]] means
                               top level has 3 elements, first is leaf, second has 2 lists of 2 elements
        :param min_length: Minimum sequence length
        :param max_length: Maximum sequence length
        :return: NucleobaseCode with random sequences
        """
        def generate_recursive(shape):
            if isinstance(shape, int):
                # Generate 'shape' number of random sequences
                result = []
                for _ in range(shape):
                    length = random.randint(min_length, max_length)
                    seq = ''.join(random.choice('ACGT') for _ in range(length))
                    result.append(seq)
                return result
            elif isinstance(shape, list):
                return [generate_recursive(item) for item in shape]
            else:
                raise TypeError(f"Invalid shape element: {shape}")
        
        data = generate_recursive(structure_shape)
        return NucleobaseCode(data)

    def __str__(self):
        """String representation showing packed storage info."""
        output = f"NucleobaseCode (packed={True}, mode={self.pack_mode})\n"
        output += f"  Codewords: {self.num_codewords}\n"
        output += f"  Max depth: {self.max_depth}\n"
        output += f"  Packed size: {len(self._packed_data)} bits ({len(self._packed_data)//8} bytes)\n"
        output += f"  Structure tokens: {len(self._structure)}\n"
        
        # Show unpacked preview
        try:
            unpacked = self.unpack()
            output += f"  Data (unpacked preview): {str(unpacked)[:200]}...\n"
        except:
            output += f"  Data: <unable to unpack>\n"
        
        return output

    def __len__(self):
        """Return number of top-level codewords."""
        return self.num_codewords

    def __getitem__(self, index):
        """Get a codeword by index (requires full unpacking)."""
        unpacked = self.unpack()
        return unpacked[index]

    def __iter__(self):
        """Iterate over codewords (requires full unpacking)."""
        unpacked = self.unpack()
        return iter(unpacked)

    def validate(self):
        """
        Validate the packed data by unpacking and checking.
        
        :return: True if valid
        """
        unpacked = self.unpack()
        self._validate_data(unpacked)
        return True
