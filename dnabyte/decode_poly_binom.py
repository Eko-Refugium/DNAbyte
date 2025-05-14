import csv
import scipy.special
import math as m
import itertools as it
import random
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import bitstring_to_bytes, compare_binary_strings
from dnabyte.auxiliary import undoreedsolomon, undoltcode, check_first_n_entries_are_zero, binarize_and_combine_first_n_elements

class DecodePolyBinom(Encode):
    """
    This class provides a decoding scheme for DNA sequences using a binomial encoding based on the positional library.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def decode_binary_codewords(self, data, params):

        check_ltcode, check_reedsolomon = True, True

        if params.inner_error_correction == 'ltcode':
            data, check_ltcode = undoltcode(data, params.index_carry_length, params.ltcode_header, params.bits_per_unit)

        if params.outer_error_correction == 'reedsolomon':
            tempdata = []
            for i, codewords in enumerate(data):
                codewords = [x for x in codewords if len(bin(x)[2:].zfill(params.bits_per_unit)) == params.bits_per_unit]
                tempdata.append(codewords)
            tempdata = [x for x in tempdata if len(x) > 0]
            data = tempdata

            data, check_reedsolomon = undoreedsolomon(data, params.bits_per_ec, params.bits_per_unit)
            

        final_data = []

        for codewords in data:
                
            if check_first_n_entries_are_zero(codewords, params.codeword_maxlength_positions):
                codewords = codewords[params.codeword_maxlength_positions:]
            else:
                bits_to_pop = binarize_and_combine_first_n_elements(codewords, params.codeword_maxlength_positions, params.bits_per_unit)
                codewords = codewords[: -bits_to_pop]
                codewords = codewords[params.codeword_maxlength_positions:]
            if check_first_n_entries_are_zero(codewords, 1):
                codewords = codewords[1:]
                binaries = [format(number, '0{}b'.format(params.bits_per_unit)) for number in codewords]
            else:
                bits_to_pop = binarize_and_combine_first_n_elements(codewords, 1, params.bits_per_unit)
                binaries = [format(number, '0{}b'.format(params.bits_per_unit)) for number in codewords]
                
                binaries[-1] = binaries[-1][(params.bits_per_unit - bits_to_pop):]
                binaries = binaries[1:]
            
            final_data.append(binaries)
            
        flattened_list = [item for sublist in final_data for item in sublist]
        concatenated_string = ''.join(flattened_list)
        check = check_ltcode & check_reedsolomon

        return concatenated_string, check

    def decode_poly_binom(self, data):
        """
        Decodes the provided data using the chain encoding scheme based on a positional library.
        """
        try:

            oligos = sorted(self.params.library.messages)
                    
            binariesinchuncks = []
            for codeword in range(len(data.data)):
                binary = []
                for i in range(len(data.data[codeword])):
                    makingtheoneslisz= [0] * len(oligos)
                    for l in range(len(data.data[codeword][i])):
                        makingtheoneslisz[oligos.index(data.data[codeword][i][l])] = 1
                    
                    counterforsum = 0
                    sumer = 0
                    onescounter = 0
                    for camlulator in range(len(makingtheoneslisz)):
                        if onescounter != self.params.sigma_amount:
                            
                            if makingtheoneslisz[camlulator] == 0:
                                sumer += scipy.special.binom(len(oligos) - camlulator - 1, self.params.sigma_amount - onescounter - 1)
                                
                            else:
                                onescounter += 1
                        counterforsum += 1     
                    binary.append(int(sumer))
                binariesinchuncks.append(binary)

            decoded_binary, valid = self.decode_binary_codewords(binariesinchuncks, self.params)
            
            info = {
                "length of decoded binary": len(decoded_binary),
                "data_length": len(data.data)
            }
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {e}")
                self.logger.error(traceback.format_exc())

            decoded_binary = None
            valid = False
            info = None

        return decoded_binary, valid, info