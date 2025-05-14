import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import create_counter_list, decimal_to_binary, MakeReedSolomonCodeSynthesis, makeltcodesynth


class EncodeMaxDensity(Encode):
    """
    This class provides maxdensity encoding for DNA sequences.

    This encoding scheme maps '00' to 'A', '01' to 'G', '10' to 'C', and '11' to 'T'. It has the theoretical maximum density 
    of 2 bits per nucleotide and is, thus, mostly used for theoretical purposes.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def calculate_codeword_parameters(self, params):
        # Step 1: calculate the number of bits required to store the number of added zeros at the end of the last codeword
        # number of bits reserved to save the number of added zeros at the end of the last (and any) codeword. At the beginning of
        # each codeword, the number of zeros is stored in zfill_bits bits. This is necessary to remove the padding zeros at the end 
        # of the last codeword during decoding.
        zfill_bits = m.ceil(m.log2(params.codeword_length))

        if params.outer_error_correction == 'reedsolomon':
            # to apply the reed-solomon error correction efficiently, the dna barcode must have a multiple of 8 bits
            adjustment = (4 - zfill_bits % 4) % 4
            params.dna_barcode_length -= adjustment
            zfill_bits += adjustment

        params.zfill_bits = zfill_bits

      # Step 2: Calculate the message length
        
        if params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits - params.index_carry_length - params.ltcode_header
        else:
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits

        # adjust the bits per codeword such that it is a multiple of 4
        if params.outer_error_correction == 'reedsolomon':
            while message_length % 4 != 0:
                message_length -= 1
                params.dna_barcode_length += 1
        


        params.message_length = message_length

        # Step 3: calculate bits per codeword
        bits_per_codeword = message_length * 2

        if params.outer_error_correction == 'reedsolomon':
            bits_per_codeword = m.floor(bits_per_codeword * params.reed_solo_percentage)
            bits_per_ec = (message_length * 2 - bits_per_codeword)

            # Adjust bits_per_codeword to be a multiple of 8
            adjustment = (8 - bits_per_codeword % 8) % 8
            bits_per_codeword += adjustment
            bits_per_ec -= adjustment
            params.bits_per_ec = bits_per_ec

        params.bits_per_codeword = bits_per_codeword

        # Step 4: Determine the length of the DNA barcode
        # if self.params.outer_error_correction == 'reedsolomon':
        #     # TODO: check if this is correct
        #     #dna_barcode_length = self.params.codeword_length - len(binary_codewords[0]) // 2
        #     dna_barcode_length = params.codeword_length - params.message_length // 2

        #     params.dna_barcode_length = dna_barcode_length

        
        return params

    def create_binary_codewords(self, data, params):
        """
        Generates binary codewords from the provided data. It supports both Reed-Solomon and LT code error correction schemes.
        """
        
        # Step 1: cut the binary string to substrings of length bits_per_codeword 
        binary_blocks = []
        for i in range(0, len(data.data), params.bits_per_codeword):
            binary_blocks.append(data.data[i:i + params.bits_per_codeword])

        number_of_zero_padding = len(binary_blocks[-1])

        binary_blocks[-1] = binary_blocks[-1].zfill(params.bits_per_codeword)
        
        # Step 2: create the binary codewords
        binary_codewords = []

        for block in binary_blocks:
            codeword = str(decimal_to_binary(0)).zfill(params.zfill_bits * 2) + block
            binary_codewords.append(codeword)

        # fill the last binary codeword with zeros

        if number_of_zero_padding == params.bits_per_codeword:
            number_of_zero_padding = 0

        binary_codewords[-1] = str(decimal_to_binary(number_of_zero_padding)).zfill(params.zfill_bits * 2) + binary_codewords[-1][params.zfill_bits * 2:]
        
        # Step 3: apply outer error correction
        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCodeSynthesis(binary_codewords, params.bits_per_ec)

        # Step 4: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcodesynth(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, 2)
        
        return binary_codewords
        
    def base4_list_to_bitstring(self, base4_list):
        """
        Converts a list of base-4 numbers to a bitstring.
        """
        bit_string = ''.join(format(x, '02b') for x in base4_list)
        return bit_string
    
    def binary_to_dna(self, binary_sequence):
        """
        Converts a binary string to a DNA string.
        """
        dna_mapping = {
            '00': 'A',
            '01': 'G',
            '10': 'C',
            '11': 'T'
        }
        dna_sequence = ''.join(dna_mapping[binary_sequence[i:i+2]] for i in range(0, len(binary_sequence), 2))
        return dna_sequence
    
    def encode_max_density(self, data):
        """
        This method encodes the provided data into DNA sequences using the specified parameters.
        """
        try:

            # Calculate the codeword parameters and save them in the params object
            self.params = self.calculate_codeword_parameters(self.params)

            # Creates binary codewords from the provided data given the specified parameters
            binary_codewords = self.create_binary_codewords(data, self.params)

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of binary codewords: {len(binary_codewords)}")
                if not all(len(codeword) == len(binary_codewords[0]) for codeword in binary_codewords):
                    raise ValueError("SANITY CHECK: Not all binary codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each binary codeword: {len(binary_codewords[0])}")

            # add the barcodes to the binary blocks
            dna_codewords = []

            for i in range(len(binary_codewords)):
                # create a barcode for each codeword on base 4 (equivalent to DNA bases)
                barcodes_base4 = create_counter_list(self.params.dna_barcode_length, 4, i)
                # create a barcode for each codeword on base 2
                barcodes_base2 = self.base4_list_to_bitstring(barcodes_base4)
                # translate the barcode to DNA
                dna_codewords.append(self.binary_to_dna(barcodes_base2 + binary_codewords[i]))

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of DNA codewords: {len(dna_codewords)}")
                if not all(len(codeword) == len(dna_codewords[0]) for codeword in dna_codewords):
                    raise ValueError("SANITY CHECK: Not all DNA codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each DNA codeword: {len(dna_codewords[0])}")   

            # create info return dictionary
            info = {
                "number_of_codewords": len(binary_codewords),
                "length_of_each_codeword": len(binary_codewords[0]),
                "barcode_length": self.params.dna_barcode_length,
                "data_length": len(data.data)
            }

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())

            dna_codewords = None
            info = {}

        return dna_codewords, info