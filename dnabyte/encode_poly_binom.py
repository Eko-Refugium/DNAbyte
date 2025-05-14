import scipy.special
import math as m
import itertools as it
import random
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import create_counter_list, complementmap
from dnabyte.auxiliary import decimal_to_binary, MakeReedSolomonCode, makeltcode

class EncodePolyBinom(Encode):
    """
    This class provides an encoding scheme for DNA sequences using a chain approach based on the poly library.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def combine_elements_randomized(self, list_of_lists):
        """
        Takes a list of lists and generates a list of lists where each new list contains one element from each sublist,
        ensuring that each element from the original sublists is used exactly once, in a randomized order.
        """
        # Ensure all sublists have the same length
        if not all(len(sublist) == len(list_of_lists[0]) for sublist in list_of_lists):
            raise ValueError("All sublists must have the same length")

        # Number of elements in each sublist
        n = len(list_of_lists[0])

        # Generate a list of indices and shuffle it
        indices = list(range(n))
        random.shuffle(indices)

        # Generate the combined lists using the shuffled indices
        combined_lists = []
        for i in indices:
            combined_list = [sublist[i] for sublist in list_of_lists]
            combined_lists.append(combined_list)

        return combined_lists

    def calculate_codeword_parameters(self, params):

        # Step 1: calculate bits per unit
        bits_per_unit = m.floor(m.log2(scipy.special.binom(len(params.library.messages), params.sigma_amount)))

        if params.outer_error_correction == 'reedsolomon':
            # Calculate the adjustment needed to make zfill_bits a multiple of 8
            adjustment = (8 + bits_per_unit % 8) % 8
            bits_per_unit -= adjustment

        params.bits_per_unit = bits_per_unit
        message_length = len(params.library.position) - 1

        # Step 2: calculate message length
        if params.inner_error_correction == 'ltcode':
            message_length = message_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1 - params.index_carry_length - params.ltcode_header
        else:
            message_length = message_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1
        
        

        # Step 3: calculate bits per codeword
        params.bits_per_codeword = message_length

        if params.outer_error_correction == 'reedsolomon':
            
            while message_length * params.bits_per_unit//8 > 255:
                params.bits_per_unit -= 8

            params.bits_per_codeword = m.floor(message_length * params.reed_solo_percentage)
            params.bits_per_ec = message_length - params.bits_per_codeword


        return params

    def create_binary_codewords(self, data, params):
        """
        Generates binary codewords from the provided data. It supports both Reed-Solomon and LT code error correction schemes.
        """
        # Step 1: create blocks of binary data
        binary_blocks = []
        for i in range(0, len(data.data), params.bits_per_unit):
            binary_blocks.append(data.data[i:i + params.bits_per_unit])
        last_block_length = len(binary_blocks[-1])
        binary_blocks[-1] = binary_blocks[-1].zfill(params.bits_per_unit)

        # Step 2: create the binary codewords
        binary_codewords = []

        for i in range(0, len(binary_blocks), params.bits_per_codeword):
            counter = 0
            one_codeword = []
            while counter < params.bits_per_codeword:
                if i + counter >= len(binary_blocks):
                    break
                one_codeword.append(binary_blocks[i+counter])
                counter += 1
            
            # add padding to the codeword
            one_codeword.insert(0, str(decimal_to_binary(0)).zfill(params.bits_per_unit))
            counterhowmanybytesadded = 0
            while len(one_codeword) <= params.bits_per_codeword:
                one_codeword.append('0' * params.bits_per_unit)
                counterhowmanybytesadded += 1
            s = str(decimal_to_binary(counterhowmanybytesadded)).zfill(params.bits_per_unit * params.codeword_maxlength_positions)
            bitsofaddedcounter = [s[i:i + params.bits_per_unit] for i in range(0, len(s), params.bits_per_unit)]
            for i in reversed(range(len(bitsofaddedcounter))):
                one_codeword.insert(0, bitsofaddedcounter[i].zfill(params.bits_per_unit))

            binary_codewords.append(one_codeword)
            if len(binary_blocks) // (params.bits_per_codeword)!= 0:
                ValueError('The last codeword is not the right length')

        # fill the last binary codeword with zeros
        if last_block_length == params.bits_per_unit:
            last_block_length = 0
        binary_codewords[-1][len(bitsofaddedcounter)] = str(decimal_to_binary(last_block_length)).zfill(params.bits_per_unit)

        # Step 3: apply outer error correction
        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCode(binary_codewords, params.bits_per_ec, params.bits_per_unit)
                
        # Step 4: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcode(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, params.bits_per_unit)

        return binary_codewords


    def encode_binom_poly(self, data):
        """
        Encodes the data unsing a binomial encoding scheme.

        Parameters:
        data (RawData): The raw data to be encoded.

        Returns:
        tuple: A tuple containing the encoded DNA codewords and additional info.
        """

        try:
            # Calculate the codeword parameters and save them in the params object
            self.params = self.calculate_codeword_parameters(self.params)

            # Creates binary codewords from the provided data given the specified parameters
            binary_codewords = self.create_binary_codewords(data, self.params)

            # Translate the binary codewords to combinations of library oligos
            sigma_amount = int(self.params.sigma_amount)
            DNA_messages = sorted(self.params.library.messages)
            DNA_generic  = self.params.library.generic
            DNA_positions = self.params.library.position
            

            combinatorial_codewords = []

            #for codeword in range(len(binary_codewords)):
            for codeword_index, binary_codeword in enumerate(binary_codewords):
                possible_oligos = []

                for block in binary_codeword:                    
                    list_of_zeros = [0] * len(DNA_messages)
                    decimal_value = int(block,2)

                    checker = 0
                    comparevalue = 0
                    one_bit_dna = []

                    for k in range(sigma_amount):
                        found = False
                        
                        while not found:
                            comparevalue = comparevalue + scipy.special.binom(len(DNA_messages) - checker - 1, sigma_amount - k - 1)
                            
                            if decimal_value < comparevalue:
                                found = True
                                list_of_zeros[checker] = 1
                                one_bit_dna.append(DNA_messages[checker])
                                comparevalue = comparevalue - scipy.special.binom(len(DNA_messages) - checker - 1, sigma_amount - k - 1)

                            checker += 1
                    
                    possible_oligos.append(one_bit_dna)

                combinatorial_codewords.append(possible_oligos)


            # translate the combinatorial codewords to a nested list of DNA oligo_trios
            dna_codewords = []
            
            #for i in range(len(combinatorial_codewords)): #codeword
            for i, combinatorial_codeword in enumerate(combinatorial_codewords):

                dna_codeword = []
                barcodes = create_counter_list(self.params.dna_barcode_length, len(DNA_messages) - 1, i)

                for j in range(len(combinatorial_codeword) + self.params.dna_barcode_length): 
                    oligo_trios = []
                    is_barcode = j < self.params.dna_barcode_length
                    
                    if is_barcode:
                        message_source = []
                        for k in range(sigma_amount):
                            message_source.append(DNA_messages[barcodes[j]])

                    else:
                        message_source = combinatorial_codeword[j - self.params.dna_barcode_length]
                        

                    
                    for k in range(len(message_source)):
                        is_even = j % 2 == 0

                        if is_even:
                            dnatriple = [
                                DNA_positions[j] + DNA_generic[0],
                                [
                                    complementmap(message_source[k][:len(message_source[k]) // 2]) + complementmap(DNA_generic[0]),
                                    message_source[k],
                                    complementmap(DNA_generic[1]) + complementmap(message_source[k][len(message_source[k]) // 2:])
                                ],
                                DNA_generic[1] + DNA_positions[j + 1]
                            ]
                        else:
                            dnatriple = [
                                complementmap(DNA_generic[0]) + complementmap(DNA_positions[j]),
                                [
                                    DNA_generic[0] + message_source[k][:len(message_source[k]) // 2],
                                    complementmap(message_source[k]),
                                    message_source[k][len(message_source[k]) // 2:] + DNA_generic[1]
                                ],
                                complementmap(DNA_positions[j + 1]) + complementmap(DNA_generic[1])
                            ]
                        
                        oligo_trios.append(dnatriple)
                    
                    dna_codeword.append(oligo_trios)
                
                dna_codewords.append(dna_codeword)
            
            info = {
                "number_of_codewords": len(dna_codewords),
                "length_of_each_codeword": len(dna_codewords[0]),
                "barcode_length": self.params.dna_barcode_length
            }

            if self.params.theory == 'yes':
                theory_output = []
                for i in range(len(dna_codewords)):
                    theory_output.append(self.combine_elements_randomized(dna_codewords[i]))
                return theory_output, info

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())
            oligo_trios = None
            info = None

        return dna_codewords, info