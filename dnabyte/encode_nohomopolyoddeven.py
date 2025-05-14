import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import create_counter_list, decimal_to_binary, MakeReedSolomonCodeSynthesis, makeltcodesynth


class EncodeNoHomoPolyOddEven(Encode):
    """
    This class provides an encoding scheme for DNA sequences that avoid homopolymers.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def binary_to_dna_custom(self, binary_sequence):
        """
        Translates a binary string to a DNA string. 

        For every even position, a 0 and 1 are translated to T and G, respectively. 
        For every odd position, a 0 and 1 are translated to C and A, respectively.
        """
        dna_sequence = []
        for i, bit in enumerate(binary_sequence):
            if i % 2 == 0:  # Even position
                if bit == '0':
                    dna_sequence.append('T')
                else:
                    dna_sequence.append('G')
            else:  # Odd position
                if bit == '0':
                    dna_sequence.append('C')
                else:
                    dna_sequence.append('A')
        return ''.join(dna_sequence)

    def dna_to_binary_custom(self, dna_sequence):
        """
        Translates a DNA string to a binary string. Inverse function of binary_to_dna_custom
        """
        binary_sequence = []
        for i, base in enumerate(dna_sequence):
            if i % 2 == 0:  # Even position
                if base == 'T':
                    binary_sequence.append('0')
                elif base == 'G':
                    binary_sequence.append('1')
            else:  # Odd position
                if base == 'C':
                    binary_sequence.append('0')
                elif base == 'A':
                    binary_sequence.append('1')
        return ''.join(binary_sequence)

    def calculate_codeword_parameters(self, params):
        # Step 1: calculate the number of bits required to store the number of added zeros at the end of the last codeword
        zfill_bits = m.ceil(m.log2(params.codeword_length))

        if params.outer_error_correction == 'reedsolomon':
            # Calculate the adjustment needed to make zfill_bits a multiple of 8
            adjustment = (8 - zfill_bits % 8) % 8
            params.dna_barcode_length -= adjustment
            zfill_bits += adjustment
            
        params.zfill_bits = zfill_bits

        # Step 2: Calculate the message length
        if params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - 2 * params.dna_barcode_length - zfill_bits - params.index_carry_length
        else:
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits

        if params.outer_error_correction == 'reedsolomon':
            while message_length % 8 != 0:
                message_length -= 1
                params.dna_barcode_length += 1   

        params.message_length = message_length

        # Step 3: calculate the bits per codeword
        bits_per_codeword = message_length
        
        if params.outer_error_correction == 'reedsolomon':
            
            bits_per_codeword = m.floor(bits_per_codeword * params.reed_solo_percentage)
            bits_per_ec = (message_length - bits_per_codeword)

            # Adjust bits_per_codeword to be a multiple of 8
            adjustment = (8 - bits_per_codeword % 8) % 8
            bits_per_codeword += adjustment
            bits_per_ec -= adjustment
            params.bits_per_ec = bits_per_ec

        params.bits_per_codeword = bits_per_codeword

        # Step 4: Determine the length of the DNA barcode

        # TODO: Check whether this is correct
        # Determine the length of the DNA barcode
        if self.params.outer_error_correction == 'reedsolomon':
            #barcode = self.params.codeword_length - len(binary_codewords[0])
            dna_barcode_length = params.codeword_length - params.message_length

            params.dna_barcode_length = dna_barcode_length

        return params

    def create_binary_codewords(self, data, params):

        # Step 1: cut the binary string to substrings of length bits_per_codeword
        binary_blocks = []
        for i in range(0, len(data.data), params.bits_per_codeword):
            binary_blocks.append(data.data[i:i + params.bits_per_codeword])
        number_of_zero_padding = len(binary_blocks[-1])
        binary_blocks[-1] = binary_blocks[-1].zfill(params.bits_per_codeword)

        # Step 2: Create the binary codewords
        binary_codewords = []

        for block in binary_blocks:
            codeword = str(decimal_to_binary(0)).zfill(params.zfill_bits) + block
            binary_codewords.append(codeword)

        # Fill the last codeword                
        if number_of_zero_padding == params.bits_per_codeword:
            number_of_zero_padding = 0

        binary_codewords[-1] = str(decimal_to_binary(number_of_zero_padding)).zfill(params.zfill_bits * 1) + binary_codewords[-1][params.zfill_bits * 1:]

        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCodeSynthesis(binary_codewords, params.bits_per_ec)
 
         # Step 3: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcodesynth(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.dna_barcode_length, 1)

        return binary_codewords

    def encode_nohomopoly_odd_even(self, data):
        """
        Encodes the provided data using the no-homopolymer odd-even encoding scheme.

        This method encodes the provided data into DNA sequences using the specified parameters.

        Parameters:
        data (RawData): The raw data to be encoded.

        Returns:
        tuple: A tuple containing the encoded DNA string, the barcode length, and additional info.
        """
        try:

            # Calculate the codeword parameters and save them in the params object
            self.params = self.calculate_codeword_parameters(self.params)

            # Create binary codewords
            binary_codewords = self.create_binary_codewords(data, self.params)

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of binary codewords: {len(binary_codewords)}")
                if not all(len(codeword) == len(binary_codewords[0]) for codeword in binary_codewords):
                    raise ValueError("SANITY CHECK: Not all binary codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each binary codeword: {len(binary_codewords[0])}")

            # create the DNA codewords
            dna_codewords = []

            # TODO: paralize this loop
            for i in range(len(binary_codewords)):
                barcode_base2 = create_counter_list(self.params.dna_barcode_length, 2, i)            
                barcode = ''.join(str(x) for x in barcode_base2)
                binary_codewords[i] = barcode + binary_codewords[i]
                dna_codewords.append(self.binary_to_dna_custom(binary_codewords[i]))

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of DNA codewords: {len(dna_codewords)}")
                if not all(len(codeword) == len(dna_codewords[0]) for codeword in dna_codewords):
                    raise ValueError("SANITY CHECK: Not all DNA codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each DNA codeword: {len(dna_codewords[0])}")   

            # create info return dictionary
            info = {
                "number_of_codewords": len(dna_codewords),
                "length_of_each_codeword": [len(codeword) for codeword in binary_codewords],
                "barcode_length": self.params.dna_barcode_length
                }

        except Exception as e:
            self.logger.error(f'An error occurred during encoding: {e}')
            self.logger.error(traceback.format_exc())

            dna_codewords = None
            info = {}

        return dna_codewords, info