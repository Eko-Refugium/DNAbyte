import random
from bitarray import bitarray
from .base import Data

class InSilicoDNA(Data):
    """
    Represents synthesized DNA sequences as oligonucleotide strings.
    
    InSilicoDNA contains a collection of DNA sequences efficiently stored using a single
    bitarray with 2-bit encoding (A=00, C=01, G=10, T=11). All sequences are concatenated
    into one bitarray with boundaries tracked via sequence lengths.

    Inherits from Data to maintain compatibility with the DNAbyte pipeline.

    It can be instantiated either:
        - by a synthesis simulator producing the synthesised DNA sequences, or
        - directly providing the list of DNA sequences as a parameter.

    :param data: List of DNA sequences (strings of nucleotides).
    :raises TypeError: If data is not a list.
    :raises ValueError: If data is empty or contains invalid DNA sequences.
    """
    
    # 2-bit packing scheme for DNA bases
    DNA_TO_BITS = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    BITS_TO_DNA = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}

    # Nucleotide complement mapping
    COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    def __init__(self, data):
        """
        Creates a DNA object from DNA sequences.
        
        Stores all sequences in a single bitarray with boundary tracking.

        :param data: List of DNA sequences (strings or tuple of (bitarray, lengths)).
        """
        # Check if data is already packed (tuple of bitarray and lengths)
        if isinstance(data, tuple) and len(data) == 2:
            self._packed_data, self._sequence_lengths = data
            if not isinstance(self._packed_data, bitarray):
                raise TypeError("Packed data must be a bitarray")
            if not isinstance(self._sequence_lengths, list):
                raise TypeError("Sequence lengths must be a list")
        else:
            # Validate and pack the data
            self._validate_dna_data(data)
            self._packed_data, self._sequence_lengths = self._pack_all_sequences(data)
        
        self.file_paths = []
        self.size = None
        
        # Calculate metrics
        self.num_sequences = len(self._sequence_lengths)
        self.total_length = sum(self._sequence_lengths)
        self.average_length = self.total_length / self.num_sequences if self.num_sequences > 0 else 0

    @property
    def data(self):
        """
        Property for backward compatibility.
        Returns the packed data tuple.
        """
        return (self._packed_data, self._sequence_lengths)
    
    def _get_sequence_bounds(self, index):
        """
        Get the bit boundaries for a sequence in the packed data.
        
        :param index: Sequence index
        :return: (start_bit, end_bit) tuple
        """
        if not 0 <= index < len(self._sequence_lengths):
            raise IndexError(f"Sequence index {index} out of range")
        
        # Calculate cumulative start position
        start_nt = sum(self._sequence_lengths[:index])
        end_nt = start_nt + self._sequence_lengths[index]
        
        # Convert to bit positions (2 bits per nucleotide)
        return start_nt * 2, end_nt * 2
    
    def _pack_all_sequences(self, sequences):
        """
        Pack all DNA sequences into a single bitarray.
        
        :param sequences: List of DNA sequence strings
        :return: Tuple of (packed_bitarray, list_of_lengths)
        """
        all_bits = []
        lengths = []
        
        for seq in sequences:
            # Convert sequence to bits
            bits = ''.join(self.DNA_TO_BITS[nucleotide] for nucleotide in seq)
            all_bits.append(bits)
            lengths.append(len(seq))
        
        # Combine all bits into one bitarray
        combined_bits = ''.join(all_bits)
        packed_data = bitarray(combined_bits)
        
        return packed_data, lengths
    
    def _unpack_sequence(self, index):
        """
        Unpack a single sequence by index.
        
        :param index: Sequence index
        :return: DNA sequence string
        """
        start_bit, end_bit = self._get_sequence_bounds(index)
        sequence_bits = self._packed_data[start_bit:end_bit]
        
        # Convert bits to DNA string
        bits_str = sequence_bits.to01()
        return ''.join(self.BITS_TO_DNA[bits_str[i:i+2]] for i in range(0, len(bits_str), 2))

    def _validate_dna_data(self, data):
        """
        Validates that the input is a proper list of DNA sequences.
        
        :param data: The data to validate
        :raises TypeError: If data is not a list
        :raises ValueError: If data is empty or contains invalid sequences
        """
        # Check if data is a list
        if not isinstance(data, list):
            raise TypeError(f"DNA data must be a list, got {type(data).__name__}")
        
        # Check if data is empty
        if len(data) == 0:
            raise ValueError("DNA data cannot be empty")
        
        # Check if all elements are valid DNA sequences
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        
        for i, sequence in enumerate(data):
            if not isinstance(sequence, str):
                raise ValueError(f"DNA sequence at index {i} must be a string, "
                               f"got {type(sequence).__name__}")
            
            if len(sequence) == 0:
                raise ValueError(f"DNA sequence at index {i} cannot be empty")
            
            # Check if all characters are valid nucleotides
            invalid_chars = set(char for char in sequence if char not in valid_nucleotides)
            if invalid_chars:
                raise ValueError(f"DNA sequence at index {i} contains invalid nucleotides: {invalid_chars}. "
                               f"Only {sorted(valid_nucleotides)} are allowed.")

    def get_sequence(self, index):
        """
        Get a specific DNA sequence by index.
        
        :param index: Index of the sequence
        :return: The DNA sequence at the specified index (as string)
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self._sequence_lengths):
            raise IndexError(f"DNA sequence index {index} out of range (0-{len(self._sequence_lengths)-1})")
        return self._unpack_sequence(index)

    def add_sequence(self, sequence):
        """
        Add a new DNA sequence to the collection.
        
        :param sequence: DNA sequence string to add
        :raises ValueError: If sequence is invalid
        """
        # Validate the sequence
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        if not isinstance(sequence, str):
            raise ValueError(f"DNA sequence must be a string, got {type(sequence).__name__}")
        
        if len(sequence) == 0:
            raise ValueError("DNA sequence cannot be empty")
        
        invalid_chars = set(char for char in sequence if char not in valid_nucleotides)
        if invalid_chars:
            raise ValueError(f"DNA sequence contains invalid nucleotides: {invalid_chars}. "
                           f"Only {sorted(valid_nucleotides)} are allowed.")
        
        # Pack and append to the main bitarray
        bits = ''.join(self.DNA_TO_BITS[nucleotide] for nucleotide in sequence)
        self._packed_data.extend(bits)
        self._sequence_lengths.append(len(sequence))
        
        # Update metrics
        self.num_sequences = len(self._sequence_lengths)
        self.total_length = sum(self._sequence_lengths)
        self.average_length = self.total_length / self.num_sequences

    def remove_sequence(self, index):
        """
        Remove a DNA sequence by index.
        
        :param index: Index of the sequence to remove
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self._sequence_lengths):
            raise IndexError(f"DNA sequence index {index} out of range (0-{len(self._sequence_lengths)-1})")
        
        # Get bounds and remove from bitarray
        start_bit, end_bit = self._get_sequence_bounds(index)
        del self._packed_data[start_bit:end_bit]
        
        # Remove from lengths
        self._sequence_lengths.pop(index)
        
        # Update metrics
        self.num_sequences = len(self._sequence_lengths)
        self.total_length = sum(self._sequence_lengths)
        self.average_length = self.total_length / self.num_sequences if self.num_sequences > 0 else 0

    def get_sequence_lengths(self):
        """
        Get lengths of all DNA sequences.
        
        :return: List of sequence lengths
        """
        return self._sequence_lengths.copy()

    def get_nucleotide_counts(self):
        """
        Get counts of each nucleotide across all sequences.
        
        :return: Dictionary with nucleotide counts
        """
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # Unpack all sequences and count
        for i in range(len(self._sequence_lengths)):
            sequence = self._unpack_sequence(i)
            for nucleotide in sequence:
                if nucleotide in counts:
                    counts[nucleotide] += 1
        return counts

    def validate(self):
        """
        Re-validates the current DNA data structure.
        
        :return: True if valid
        :raises: Various exceptions if invalid
        """
        # Unpack all sequences and validate
        sequences = [self._unpack_sequence(i) for i in range(len(self._sequence_lengths))]
        self._validate_dna_data(sequences)
        return True

    @staticmethod
    def complement(sequence):
        """
        Get the complement of a DNA sequence.
        
        :param sequence: DNA sequence string
        :return: Complement sequence
        """
        return ''.join(InSilicoDNA.COMPLEMENT[base] for base in sequence)
    
    @staticmethod
    def reverse_complement(sequence):
        """
        Get the reverse complement of a DNA sequence.
        
        :param sequence: DNA sequence string
        :return: Reverse complement sequence
        """
        return InSilicoDNA.complement(sequence)[::-1]

    @staticmethod
    def random(m, n):
        """
        Generate m random nucleotide sequences of length n.

        Parameters:
        -----------
        m : int
            Number of sequences to generate
        n : int
            Length of each sequence

        Returns:
        --------
        InSilicoDNA
            An InSilicoDNA object containing m random DNA sequences, each of length n

        """
        if m <= 0:
            raise ValueError("Number of sequences (m) must be positive")
        if n <= 0:
            raise ValueError("Sequence length (n) must be positive")

        sequences = []
        for _ in range(m):
            sequence = ''.join(random.choice('ACGT') for _ in range(n))
            sequences.append(sequence)

        return InSilicoDNA(sequences)

    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Number of sequences: {self.num_sequences}\n"
        output += f"Total length: {self.total_length} nucleotides\n"
        output += f"Average length: {self.average_length:.1f} nucleotides\n"
        output += f"Packed bitarray size: {len(self._packed_data)} bits = {len(self._packed_data) // 8} bytes\n"
        
        # Add file paths if they exist
        if hasattr(self, 'file_paths') and self.file_paths:
            output += f"File paths: {self.file_paths}\n"
        
        if hasattr(self, 'size') and self.size:
            output += f"Original size: {self.size} bytes\n"
        
        # Show nucleotide distribution
        counts = self.get_nucleotide_counts()
        total_nucleotides = sum(counts.values())
        if total_nucleotides > 0:
            output += "Nucleotide distribution: "
            output += ", ".join(f"{nt}: {count} ({100*count/total_nucleotides:.1f}%)" 
                              for nt, count in counts.items())
            output += "\n"
        
        # Show first few sequences as preview
        preview_seqs = [self._unpack_sequence(i) for i in range(min(3, len(self._sequence_lengths)))]
        output += f"DATA (first {len(preview_seqs)} sequences): {str(preview_seqs)[:100]}...\n"
        return output

    def __len__(self):
        """Return the number of DNA sequences."""
        return len(self._sequence_lengths)

    def __getitem__(self, index):
        """Allow indexing to get DNA sequences."""
        return self.get_sequence(index)

    def __iter__(self):
        """Allow iteration over DNA sequences."""
        for i in range(len(self._sequence_lengths)):
            yield self._unpack_sequence(i)
