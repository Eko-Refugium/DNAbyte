import math as m
import traceback

from dnabyte.error_correction.auxiliary import undoreedsolomonsynthesis, undoltcodesynth
from dnabyte.encoding.auxiliary import split_string

def decode(data, params, logger=None):
    """
    Decodes the provided data using the no-homopolymer odd-even decoding scheme.

    This method decodes the provided DNA sequences into binary data using the specified parameters.

    Parameters:
    data (SequencedData): The sequenced data to be decoded.

    Returns:
    tuple: A tuple containing the decoded binary data, a validity checker, and additional info.
    """
    try:
        binary_codewords = []
        
        for i in range(len(data.data)):
            binary_codewords.append(dna_to_binary_custom(data.data[i]))
        
        decoded_binary, valid = decode_binary_codewords(binary_codewords, params)

        # Additional info
        info = {
            "Length of decoded data": len(decoded_binary),
        }

    except Exception as e:
        if logger:
            logger.error(f"Error during decoding: {str(e)}")
            logger.error(traceback.format_exc())
        
        decoded_binary = None
        valid = False
        info = None

    return decoded_binary, valid, info

def dna_to_binary_custom(dna_string):
    """
    Converts a DNA string to a binary string using a custom mapping to avoid homopolymers.

    Parameters:
    dna_string (str): The DNA string to be converted.

    Returns:
    str: The corresponding binary string.
    """
    binary_string = []
    if dna_string[0] == 'T' or dna_string[0] == 'G':
        for i, base in enumerate(dna_string):
            if i % 2 == 0:  # Even position
                if base == 'T':
                    binary_string.append('0')
                elif base == 'G':
                    binary_string.append('1')
            else:  # Odd position
                if base == 'C':
                    binary_string.append('0')
                elif base == 'A':
                    binary_string.append('1')
    else:
        for i, base in enumerate(dna_string):
            if i % 2 == 0:
                if base == 'C':
                    binary_string.append('0')
                elif base == 'A':
                    binary_string.append('1')
            else:
                if base == 'T':
                    binary_string.append('0')
                elif base == 'G':
                    binary_string.append('1')
    return ''.join(binary_string)


def decode_binary_codewords(data, params):

    check_ltcode, check_reedsolomon = True, True
    print(data)
    if hasattr(params, 'inner_error_correction') and params.inner_error_correction == 'ltcode':
        data, check_ltcode = undoltcodesynth(data, params.index_carry_length, params.ltcode_header, 1)
    print(data)
    if hasattr(params, 'outer_error_correction') and params.outer_error_correction == 'reedsolomon':
        data, check_reedsolomon = undoreedsolomonsynthesis(data, params.bits_per_ec)
    print(data)
    final_data = []

    for i in range(len(data)):
        informationoflength, actualdata = split_string(data[i], (params.zfill_bits) * 1)
        informationoflength = int(informationoflength, 2)
        if informationoflength != 0:
            actualdata = actualdata[(params.bits_per_codeword - informationoflength):]
        final_data.append(actualdata)
    combined_string = ''.join(final_data)
    
    check = check_ltcode & check_reedsolomon
    return combined_string, check


    