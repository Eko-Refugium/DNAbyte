import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.error_correction.auxiliary import undoreedsolomonsynthesis, undoltcodesynth
from dnabyte.encoding.auxiliary import split_string

def decode(data, params, logger=None):
    """
    Decodes the provided corrected data using maxdensity decoding.
    """
    try:

        # Convert DNA sequences to binary strings
        binary_strings = []
        for i in range(len(data.data)):
            binary_strings.append(dna_to_binary(data.data[i]))
        decoded_binary, valid = recreate_binary_codewords(binary_strings, params)

        info = {
            "number_of_codewords": len(binary_strings),
            "data_length": len(data.data)
        }

    except Exception as e:
        if logger:
            logger.error(f"Error during decoding: {e}")
            logger.error(traceback.format_exc())

        decoded_binary = None
        valid = False
        info = {}

    return decoded_binary, valid, info

def dna_to_binary(dna_string):
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

def recreate_binary_codewords(data, params):

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



        
    