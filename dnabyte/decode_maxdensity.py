import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import undoreedsolomonsynthesis, undoltcodesynth, split_string

class DecodeMaxDensity(Encode):
    """
    This class provides maxdensity decoding for DNA sequences based on the specified parameters.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def dna_to_binary(self, dna_string):
        """
        Converts a DNA string to a binary string.
        """
        binary_mapping = {
            'A': '00',
            'G': '01',
            'C': '10',
            'T': '11'
        }
        binary_string = ''.join(binary_mapping[base] for base in dna_string)
        return binary_string
    
    def recreate_binary_codewords(self, data, params):

        check_ltcode, check_reedsolomon = True, True

        # Step 1: undo inner error correction
        if params.inner_error_correction == 'ltcode':
            data, check_ltcode = undoltcodesynth(data, params.index_carry_length, params.ltcode_header, 2)
        # Step 2: undo couter error correction
        if params.outer_error_correction == 'reedsolomon':
            data, check_reedsolomon = undoreedsolomonsynthesis(data, params.bits_per_ec)

        # Step 3: remove the zfill bits
        final_data = []

        for i in range(len(data)):
            zfill_length, actualdata = split_string(data[i], (params.zfill_bits) * 2)
            zfill_length = int(zfill_length, 2)
            if zfill_length != 0:
                actualdata = actualdata[(params.bits_per_codeword - zfill_length):]
            final_data.append(actualdata)

        concatenated_string = ''.join(final_data)
        check = check_ltcode & check_reedsolomon

        return concatenated_string, check
    
    def decode_maxdensity(self, data):
        """
        Decodes the provided corrected data using maxdensity decoding.
        """
        try:

            # Convert DNA sequences to binary strings
            binary_strings = []
            for i in range(len(data.data)):
                binary_strings.append(self.dna_to_binary(data.data[i]))
            decoded_binary, valid = self.recreate_binary_codewords(binary_strings, self.params)

            info = {
                "number_of_codewords": len(binary_strings),
                "data_length": len(data.data)
            }

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {e}")
                self.logger.error(traceback.format_exc())

            decoded_binary = None
            valid = False
            info = {}

        return decoded_binary, valid, info

            
        