import numpy as np
import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.auxiliary import create_counter_list, complementmap
from dnabyte.auxiliary import decimal_to_binary, MakeReedSolomonCode, makeltcode

class EncodeLinearChain(Encode):
    """
    This class provides an encoding scheme for DNA sequences using a linear chain approach.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def nest_list_structure(self, lst):
        """
        Nests a list structure by recursively grouping elements into pairs.

        Parameters:
        lst (list): The list to be nested.

        Returns:
        list: The nested list structure.
        """
        if len(lst) <= 1:
            return lst
        nested = lst[0]
        for i in range(1, len(lst)):
            nested = [nested, lst[i]]
        return nested


    def contains_motif_or_complement(self, sublist, string, i):
        """
        Checks if a string contains a motif or its complement from a sublist.

        Parameters:
        sublist (list): The sublist of motifs.
        string (str): The string to be checked.
        i (int): The index of the string.

        Returns:
        bool: True if the string contains a motif or its complement, False otherwise.
        """
        for s in sublist:
            if i % 2 == 0:
                if s != sublist[-1]:
                    if s[:len(s)//2] in string or complementmap(s)[:len(s)//2] in string or s[len(s)//2:] in string or complementmap(s)[len(s)//2:] in string:
                        return True
                else:
                    if s[:len(s)//2] in string or complementmap(s)[:len(s)//2] in string:
                        return True
                
            else:
                
                if s!=sublist[-1]:
                    if s[:len(s)//2] in string or complementmap(s)[:len(s)//2] in string or s[len(s)//2:] in string or complementmap(s)[len(s)//2:] in string:
                            return True
                else:
                    if s[len(s)//2:] in string or complementmap(s)[len(s)//2:] in string:
                        return True
        return False

    def process_strings(self, strings):
        """
        Processes a list of strings by grouping them into sublists based on motifs and their complements.

        Parameters:
        strings (list): The list of strings to be processed.

        Returns:
        list: The processed list of strings grouped into sublists.
        """
        result = []
        i = 0
        
        while i + 3 <= len(strings):
            sublist = strings[i:i+3]
            i += 3
            
            while i < len(strings):
                if not self.contains_motif_or_complement(sublist, strings[i], i):
                    sublist.append(strings[i])
                    i += 1
                else:
                    break
            
            result.append(sublist)
        
        # If there are fewer than three strings remaining, leave them in the outer list
        if i < len(self,strings):
            for j in strings[i:]:
                result.append(j)
        
        return result

    def calculate_codeword_parameters(self, params):

        # Step 1: Calculate the number of bits per unit
        bits_per_unit = m.floor(m.log2(len(params.library.library)))

        # Adjust for Reed-Solomon error correction
        if params.outer_error_correction == 'reedsolomon':
            bits_per_unit -= bits_per_unit % 8

        params.bits_per_unit = bits_per_unit

        # Step 2: Calculate the message length
        if params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1 - params.index_carry_length - params.ltcode_header
        else:
            message_length = params.codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1

        params.message_length = message_length

        # TODO: rename these variables such that they are consistent between different encoding schemes
        if params.outer_error_correction == 'reedsolomon':
            
            while message_length * params.bits_per_unit//8 > 255:
                message_length -= 1
        # calculate the number of bits for the info
        reed_solo_info = m.floor(params.message_length * params.reed_solo_percentage)
        # calculate the number of bits for the correction
        reed_solo_corr = params.message_length - reed_solo_info

        params.reed_solo_info = reed_solo_info
        params.reed_solo_corr = reed_solo_corr

        return params

    def create_binary_codewords(self, data, params):

        # Step 1: create binary blocks of length bits_per_unit
        binary_blocks = [data.data[i:i + params.bits_per_unit] for i in range(0, len(data.data), params.bits_per_unit)]
        len_of_last_block = len(binary_blocks[-1])
        binary_blocks[-1] = binary_blocks[-1].zfill(params.bits_per_unit)

        # fill the last binary block with zeros
        if len_of_last_block == params.bits_per_unit:
            len_of_last_block = 0

        # TODO: implement a version that calculates the parameters first and then creates the binary codewords irrespective of 
        # whether there is an outer error correction or not

        if params.outer_error_correction == 'reedsolomon':

            
            # create binary codewords
            binary_codewords = []
            for i in range(0, len(binary_blocks), params.reed_solo_info):

                counter = 0
                one_codeword = []
                while counter < params.reed_solo_info:
                    if i + counter >= len(binary_blocks):
                        break
                    one_codeword.append(binary_blocks[i + counter])
                    counter += 1
                
                one_codeword.insert(0, str(decimal_to_binary(0)).zfill(params.bits_per_unit))

                counterhowmanybytesadded = 0
                while len(one_codeword) <= params.reed_solo_info:
                    one_codeword.append('0' * params.bits_per_unit)
                    counterhowmanybytesadded += 1


                s = str(decimal_to_binary(counterhowmanybytesadded)).zfill(params.bits_per_unit * params.codeword_maxlength_positions)
                bitsofaddedcounter = [s[i:i + params.bits_per_unit] for i in range(0, len(s), params.bits_per_unit)]
                for i in reversed(range(len(bitsofaddedcounter))):
                    one_codeword.insert(0, bitsofaddedcounter[i].zfill(params.bits_per_unit))

                binary_codewords.append(one_codeword)
                if len(binary_blocks) // (params.reed_solo_info) != 0:
                    ValueError('The last codeword is not the right length')
            binary_codewords[-1][len(bitsofaddedcounter)] = str(decimal_to_binary(len_of_last_block)).zfill(params.bits_per_unit)
            binary_codewords = MakeReedSolomonCode(binary_codewords, params.reed_solo_corr, params.bits_per_unit)
        
        else:
            codeword_length = params.message_length

            # create binary codewords
            binary_codewords = []
            for i in range(0, len(binary_blocks), codeword_length):
                counter = 0
                one_codeword = []
                while counter < codeword_length:
                    if i + counter >= len(binary_blocks):
                        break
                    one_codeword.append(binary_blocks[i+counter])
                    counter += 1
                
                one_codeword.insert(0, str(decimal_to_binary(0)).zfill(params.bits_per_unit))
                counterhowmanybytesadded = 0
                while len(one_codeword) <= codeword_length:
                    one_codeword.append('0' * params.bits_per_unit)
                    counterhowmanybytesadded += 1
                s = str(decimal_to_binary(counterhowmanybytesadded)).zfill(params.bits_per_unit*params.codeword_maxlength_positions)
                bitsofaddedcounter=[s[i:i + params.bits_per_unit] for i in range(0, len(s), params.bits_per_unit)]
                for i in reversed(range(len(bitsofaddedcounter))):
                    one_codeword.insert(0, bitsofaddedcounter[i].zfill(params.bits_per_unit))

                binary_codewords.append(one_codeword)
                if len(binary_blocks)//(codeword_length)!= 0:
                    ValueError('The last codeword is not the right length')
            binary_codewords[-1][1] = str(decimal_to_binary(len_of_last_block)).zfill(params.bits_per_unit)

        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcode(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, params.bits_per_unit)            
            
        return binary_codewords


    def encode_linear_chain(self, data):
        """
        Encodes the provided data using the linear chain encoding scheme.

        This method encodes the provided data into DNA sequences using the specified parameters.

        Parameters:
        data (RawData): The raw data to be encoded.

        Returns:
        tuple: A tuple containing the encoded DNA codewords and additional info.
        """
        try:
            DNAs = self.params.library.library
            oligolen = len(DNAs[0])
            DNAs.sort()

            # Calculate the codeword parameters
            self.params = self.calculate_codeword_parameters(self.params)

            # Create binary codewords
            binary_codewords = self.create_binary_codewords(data, self.params)
            # Translate binary codewords to decimal codewords
            decimal_codewords = []
            for codeword in binary_codewords:
                decimal_codeword = [int(binary_codewords, 2) for binary_codewords in codeword]
                decimal_codewords.append(decimal_codeword)
            
            # Translate decimal codewords to DNA codewords
            dna_codewords = []

            # TODO: paralize this loop
            for i in range(len(decimal_codewords)):
                barcodes = create_counter_list(self.params.dna_barcode_length, len(DNAs) - 1, i)
                dna_codeword = []

                # TODO: explain the logic of the following code
                for j in range(len(decimal_codewords[i]) + self.params.dna_barcode_length - 1):
                    if j < self.params.dna_barcode_length:
                        if j + 1 == self.params.dna_barcode_length:
                            dna_codeword.append(DNAs[barcodes[j]])
                            connector = complementmap(DNAs[decimal_codewords[i][j - self.params.dna_barcode_length + 1]][:oligolen // 2]) + complementmap(DNAs[barcodes[j]][oligolen // 2:])
                            dna_codeword.append(connector)
                        else:
                            dna_codeword.append(DNAs[barcodes[j]])
                            connector = complementmap(DNAs[barcodes[j + 1]][:oligolen // 2]) + complementmap(DNAs[barcodes[j]][oligolen // 2:])
                            dna_codeword.append(connector)
                    else:
                        dna_codeword.append(DNAs[decimal_codewords[i][j - self.params.dna_barcode_length]])
                        connector = complementmap(DNAs[decimal_codewords[i][j - self.params.dna_barcode_length + 1]][:oligolen // 2]) + complementmap(DNAs[decimal_codewords[i][j - self.params.dna_barcode_length]][oligolen // 2:])
                        dna_codeword.append(connector)
                dna_codeword.append(DNAs[decimal_codewords[i][-1]])
                dna_codewords.append(dna_codeword)

            # Create additional info
            info = {
                "number_of_codewords": len(dna_codewords),
                "length_of_each_codeword": [len(dna_codewords[0])],
                "barcode_length": self.params.dna_barcode_length,
                "data_length": len(data.data)
            }

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())
            dna_codewords = None
            info = None

            

        return dna_codewords, info



