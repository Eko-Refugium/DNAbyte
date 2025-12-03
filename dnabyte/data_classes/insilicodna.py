import random
from .base import Data

class InSilicoDNA(Data):
    """
    Represents synthesized DNA sequences as oligonucleotide strings.
    
    InSilicoDNA contains a collection of DNA sequences (strings of 'A', 'C', 'G', 'T')
    that represent synthesized DNA. This class is used for DNA synthesis simulation
    and applying error channels like storage simulation and PCR amplification.

    Inherits from Data to maintain compatibility with the DNAbyte pipeline.

    It can be instantiated either:
        - by a synthesis simulator producing the synthesised DNA sequences, or
        - directly providing the list of DNA sequences as a parameter.

    :param data: List of DNA sequences (strings of nucleotides).
    :raises TypeError: If data is not a list.
    :raises ValueError: If data is empty or contains invalid DNA sequences.
    """
    
    def __init__(self, data):
        """
        Creates a DNA object from DNA sequences.

        :param data: List of DNA sequences (strings).
        """

        # Initialize directly from DNA sequences
        self._validate_dna_data(data)
        self.data = data
        self.file_paths = []
        self.size = None
        
        # Calculate metrics
        self.num_sequences = len(self.data) if self.data else 0
        self.total_length = sum(len(seq) for seq in self.data) if self.data else 0
        self.average_length = self.total_length / self.num_sequences if self.num_sequences > 0 else 0

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
            #print((isinstance(sequence, list) and len(sequence) == 2), 'sequence in insilicodna validation')
            if not isinstance(sequence, str) and not (isinstance(sequence, list) and len(sequence) == 2):
                raise ValueError(f"DNA sequence at index {i} must be a string, "
                               f"got {type(sequence).__name__}")
            
            if len(sequence) == 0:
                raise ValueError(f"DNA sequence at index {i} cannot be empty")
            
            # Check if all characters are valid nucleotides
            if isinstance(sequence, list) and len(sequence) == 2:
                invalid_chars = set()
                for seq in sequence:
                    invalid_chars.update(char for char in seq if char not in valid_nucleotides)
            else:
                invalid_chars = set(char for char in sequence if char not in valid_nucleotides)
            if invalid_chars:
                raise ValueError(f"DNA sequence at index {i} contains invalid nucleotides: {invalid_chars}. "
                               f"Only {sorted(valid_nucleotides)} are allowed.")

    def _flatten_encoded_data(self, encoded_data):
        """
        Converts nested encoded data structure to flat list of DNA sequences.
        
        :param encoded_data: Nested list structure from EncodedData
        :return: List of DNA sequence strings
        """
        sequences = []
        for codeword in encoded_data:
            # Flatten each codeword and join nucleotides into sequences
            flattened = self._flatten_recursive(codeword)
            sequence = ''.join(flattened)
            if sequence:  # Only add non-empty sequences
                sequences.append(sequence)
        return sequences

    def _flatten_recursive(self, structure):
        """
        Recursively flattens nested structure to get nucleotides.
        
        :param structure: Nested structure to flatten
        :return: List of nucleotides
        """
        result = []
        if isinstance(structure, list):
            for item in structure:
                result.extend(self._flatten_recursive(item))
        else:
            result.append(structure)
        return result

    def get_sequence(self, index):
        """
        Get a specific DNA sequence by index.
        
        :param index: Index of the sequence
        :return: The DNA sequence at the specified index
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self.data):
            raise IndexError(f"DNA sequence index {index} out of range (0-{len(self.data)-1})")
        return self.data[index]

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
        
        # Add the sequence and update metrics
        self.data.append(sequence)
        self.num_sequences = len(self.data)
        self.total_length = sum(len(seq) for seq in self.data)
        self.average_length = self.total_length / self.num_sequences

    def remove_sequence(self, index):
        """
        Remove a DNA sequence by index.
        
        :param index: Index of the sequence to remove
        :raises IndexError: If index is out of range
        """
        if not 0 <= index < len(self.data):
            raise IndexError(f"DNA sequence index {index} out of range (0-{len(self.data)-1})")
        
        self.data.pop(index)
        self.num_sequences = len(self.data)
        self.total_length = sum(len(seq) for seq in self.data) if self.data else 0
        self.average_length = self.total_length / self.num_sequences if self.num_sequences > 0 else 0

    def get_sequence_lengths(self):
        """
        Get lengths of all DNA sequences.
        
        :return: List of sequence lengths
        """
        return [len(seq) for seq in self.data]

    def get_nucleotide_counts(self):
        """
        Get counts of each nucleotide across all sequences.
        
        :return: Dictionary with nucleotide counts
        """
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for sequence in self.data:
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
        self._validate_dna_data(self.data)
        return True


    def generate_random_sequences(m, n):
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
        list
            List of m random DNA sequences, each of length n

        """
        if m <= 0:
            raise ValueError("Number of sequences (m) must be positive")
        if n <= 0:
            raise ValueError("Sequence length (n) must be positive")

        sequences = []

        for _ in range(m):
            sequence = ''.join(random.choice('ACGT') for _ in range(n))
            sequences.append(sequence)

        return sequences


    def __str__(self):
        output = f"Type: {type(self).__name__}\n"
        output += f"Number of sequences: {self.num_sequences}\n"
        output += f"Total length: {self.total_length} nucleotides\n"
        output += f"Average length: {self.average_length:.1f} nucleotides\n"
        
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
        
        output += f"DATA: {str(self.data)[:100]}...\n"
        return output

    def __len__(self):
        """Return the number of DNA sequences."""
        return len(self.data)

    def __getitem__(self, index):
        """Allow indexing to get DNA sequences."""
        return self.get_sequence(index)

    def __iter__(self):
        """Allow iteration over DNA sequences."""
        return iter(self.data)
