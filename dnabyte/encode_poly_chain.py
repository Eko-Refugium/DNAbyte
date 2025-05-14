import numpy as np
import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import create_counter_list, complementmap
from dnabyte.auxiliary import decimal_to_binary, MakeReedSolomonCode, makeltcode

class EncodePolyChain(Encode):
    """
    This class provides an encoding scheme for DNA sequences using a chain approach based on the poly library.
    """
    def __init__(self, params, logger=None): 
        self.params = params
        self.logger = logger

    def calculate_codeword_parameters(self, params):

        # Step 1: calculate the bits per unit
        bits_per_unit = m.floor(m.log2(len(params.library.messages)))

        if params.outer_error_correction == 'reedsolomon':
            # Calculate the adjustment needed to make zfill_bits a multiple of 8
            adjustment = (8 + bits_per_unit % 8) % 8
            
            # params.dna_barcode_length -= adjustment
            bits_per_unit -= adjustment

        params.bits_per_unit = bits_per_unit

        # Step 2: determine the message length
        max_codeword_length = len(params.library.position) - 1

        if params.inner_error_correction == 'ltcode':
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1 - params.index_carry_length -params.ltcode_header
        else:
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions-1

        params.message_length = message_length

        # Step 3: calculate the bits per codeword
        params.bits_per_codeword = message_length

        if params.outer_error_correction == 'reedsolomon':
            
            while message_length * params.bits_per_unit//8 > 255:
                params.message_length -= 1
                params.barcode_length += 1
            
            params.bits_per_codeword = m.floor(message_length * params.reed_solo_percentage)        
            params.bits_per_ec = message_length - params.bits_per_codeword

        return params


    def create_binary_codewords(self, data, params):
        """
        Generates binary codewords from the provided data. It supports both Reed-Solomon and LT code error correction schemes.
        """
        # Step 1: create the binary blocks
        binary_blocks = []
        for i in range(0, len(data), params.bits_per_unit):
            binary_blocks.append(data[i:i + params.bits_per_unit])
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
            
            one_codeword.insert(0, str(decimal_to_binary(0)).zfill(params.bits_per_unit))
            counterhowmanybytesadded = 0
            while len(one_codeword) <= params.bits_per_codeword:
                one_codeword.append('0' * params.bits_per_unit)
                counterhowmanybytesadded += 1
            s = str(decimal_to_binary(counterhowmanybytesadded)).zfill(params.bits_per_unit * params.codeword_maxlength_positions)
            bitsofaddedcounter=[s[i:i+params.bits_per_unit] for i in range(0, len(s), params.bits_per_unit)]
            for i in reversed(range(len(bitsofaddedcounter))):
                one_codeword.insert(0, bitsofaddedcounter[i].zfill(params.bits_per_unit))

            binary_codewords.append(one_codeword)
            if len(binary_blocks) // (params.bits_per_codeword)!= 0:
                ValueError('The last codeword is not the right length')

        # fill the last binary codeword with zeros
        if last_block_length == params.bits_per_unit:
            last_block_length = 0
        binary_codewords[-1][len(bitsofaddedcounter)] = str(decimal_to_binary(last_block_length)).zfill(params.bits_per_unit)

        # Step 4: apply outer error correction
        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCode(binary_codewords, params.bits_per_ec, params.bits_per_unit) 

        # Step 5: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcode(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, params.bits_per_unit)
            
        return binary_codewords

    def encode_chain_poly(self, data):
        """
        Encodes the provided data using the chain encoding scheme based on a positional library.

        This method encodes the provided data into DNA sequences using the specified parameters.

        Parameters:
        data (RawData): The raw data to be encoded.

        Returns:
        tuple: A tuple containing the encoded DNA codewords and additional info.
        """
        try:
            # Calculate the codeword parameters and save them in the params object
            self.params = self.calculate_codeword_parameters(self.params)

            # Creates binary codewords from the provided data given the specified parameters
            binary_codewords = self.create_binary_codewords(data.data, self.params)

            # Translate the decimal codewords to DNA codewords
            DNA_messages = self.params.library.messages
            DNA_generic  = self.params.library.generic
            DNA_positions = self.params.library.position
            DNA_messages.sort()

            dna_codewords = []

            # TODO: parallelize this loop
            #for i in range(len(decimal_codewords)):

            for i, binary_codeword in enumerate(binary_codewords):
                dna_codeword = []
                barcodes = create_counter_list(self.params.dna_barcode_length, len(DNA_messages) - 1, i)

                decimal_codeword = [int(binary_codeword[j], 2) for j in range(len(binary_codeword))]
                
                
                for j in range(len(decimal_codeword) + self.params.dna_barcode_length):
                    is_barcode = j < self.params.dna_barcode_length
                    is_even = j % 2 == 0

                    # Determine the source of the DNA message
                    if is_barcode:
                        message_index = barcodes[j]
                    else:
                        # decimal_codeword = int(binary_codeword[j - self.params.dna_barcode_length], 2)
                        message_index = decimal_codeword[j - self.params.dna_barcode_length]

                    # Encode the DNA triple
                    if is_even: # Information is on the 5' to 3' strand
                        dnatripple = [
                            DNA_positions[j] + DNA_generic[0],
                            [
                                complementmap(DNA_messages[message_index][:len(DNA_messages[message_index]) // 2]) + complementmap(DNA_generic[0]),
                                DNA_messages[message_index],
                                complementmap(DNA_generic[1]) + complementmap(DNA_messages[message_index][len(DNA_messages[message_index]) // 2:])
                            ],
                            DNA_generic[1] + DNA_positions[j + 1]
                        ]
                    
                    else: # Information is on the 3' to 5' strand
                        dnatripple = [
                            complementmap(DNA_generic[0]) + complementmap(DNA_positions[j]),
                            [
                                DNA_generic[0] + DNA_messages[message_index][:len(DNA_messages[message_index]) // 2],
                                complementmap(DNA_messages[message_index]),
                                DNA_messages[message_index][len(DNA_messages[message_index]) // 2:] + DNA_generic[1]
                            ],
                            complementmap(DNA_positions[j + 1]) + complementmap(DNA_generic[1])
                        ]

                    dna_codeword.append(dnatripple)

                dna_codewords.append(dna_codeword)


            # create info return dictionary
            info = {
                "number_of_codewords": len(dna_codewords),
                "length_of_each_codeword": [len(dna_codewords[0]) for codeword in dna_codewords],
                "barcode_length": self.params.dna_barcode_length,
                "data_length": len(data.data)
            }

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error in encode_chain_poly: {e}")
                self.logger.error(traceback.format_exc())
            dna_codewords = None
            info = None

        return dna_codewords, info
