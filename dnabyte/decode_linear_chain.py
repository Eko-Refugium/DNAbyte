import traceback
import math as m

from dnabyte.encode import Encode
from dnabyte.auxiliary import undoreedsolomon, undoltcode, check_first_n_entries_are_zero, binarize_and_combine_first_n_elements
from dnabyte.oligo import create_positional_libraries

class DecodeLinearChain(Encode):
    """
    This class provides decoding for DNA sequences based on the linear encoding with chain assembly.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def decode_binary_codewords(self, data, params):

        check_ltcode, check_reedsolomon = True, True

        # undo inner error correction
        if params.inner_error_correction == 'ltcode':
            data, check_ltcode = undoltcode(data,params.index_carry_length, params.ltcode_header, params.bits_per_unit)

        # undo outer error correction
        if params.outer_error_correction == 'reedsolomon':
            tempdata = []
            for i, codewords in enumerate(data):
                codewords = [x for x in codewords if len(bin(x)[2:].zfill(params.bits_per_unit)) == params.bits_per_unit]
                tempdata.append(codewords)
            tempdata = [x for x in tempdata if len(x) > 0]
            data = tempdata

            data, check_reedsolomon = undoreedsolomon(data, params.reed_solo_corr, params.bits_per_unit)

        corrected_list = []
        for codewords in data:

                if check_first_n_entries_are_zero(codewords, params.codeword_maxlength_positions):
                    # remove the first n elements
                    codewords = codewords[params.codeword_maxlength_positions:]
                else:
                    bits_to_pop = binarize_and_combine_first_n_elements(codewords, params.codeword_maxlength_positions, params.bits_per_unit)
                    codewords = codewords[:-bits_to_pop]
                    codewords = codewords[params.codeword_maxlength_positions:]

                if check_first_n_entries_are_zero(codewords, 1):
                    # remove the first element and convert the rest to binary
                    codewords = codewords[1:]
                    binaries = [format(number, '0{}b'.format(params.bits_per_unit)) for number in codewords]
                else:
                    bits_to_pop = binarize_and_combine_first_n_elements(codewords, 1, params.bits_per_unit)
                    binaries = [format(number, '0{}b'.format(params.bits_per_unit)) for number in codewords]
                    binaries[-1] = binaries[-1][(params.bits_per_unit - bits_to_pop):]
                    binaries = binaries[1:]
                corrected_list.append(binaries)
            
        # creating return objects    
        concatenated_string = ''.join([item for sublist in corrected_list for item in sublist])
        check = check_ltcode & check_reedsolomon

        return concatenated_string, check

    def decode_linear_chain(self, data):
        """
        Converts a DNA string to a binary string.
        """
        try:
            DNAsmessages = sorted(self.params.library.library)

            # TODO: sanity checks

            # TODO: find a proper name for decimal_decode
            decimal_decode = []
            for codewords in data.data:
                temp = []
                for i in range(len(codewords)):
                    temp.append(DNAsmessages.index(codewords[i]))	
                decimal_decode.append(temp)

            decoded_binary, valid = self.decode_binary_codewords(decimal_decode, self.params)

            # TODO: sanity checks

            info = {
                "number_of_codewords": len(decoded_binary),
                "data_length": len(data.data)
            }

        except Exception as e:
            traceback.print_exc()
            decoded_binary = None
            valid = False
            info = None

        return decoded_binary, valid, info

