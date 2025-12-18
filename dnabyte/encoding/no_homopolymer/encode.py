import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.encoding.auxiliary import create_counter_list, decimal_to_binary, check_parameter, check_library
from dnabyte.error_correction.auxiliary import MakeReedSolomonCodeSynthesis, makeltcodesynth
from dnabyte.encoding.no_homopolymer.decode import decode as decode_function
from dnabyte.encoding.no_homopolymer.process import process as process_function


class NoHomoPoly(Encode):
    """
    This class provides an encoding scheme for DNA sequences that avoid homopolymers.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
        self._decode_function = decode_function
        self._process_function = process_function

    def encode(self, data):
        """
        Encodes the provided data using the no-homopolymer odd-even encoding scheme.

        This method encodes the provided data into DNA sequences using the specified parameters.

        Parameters:
        data (BinaryCode): The raw data to be encoded.

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
                dna_codeword = self.binary_to_dna_custom(binary_codewords[i])
                dna_codewords.append(dna_codeword)
            
            # Calculate the actual barcode length in DNA (nucleotides) and in binary (after decoding)
            # The first DNA codeword's barcode length can be determined by encoding just the barcode
            if dna_codewords:
                # Store the original binary barcode length (before DNA encoding) for LT code
                original_binary_barcode_length = self.params.dna_barcode_length
                
                barcode_base2_first = create_counter_list(self.params.dna_barcode_length, 2, 0)
                barcode_first = ''.join(str(x) for x in barcode_base2_first)
                dna_barcode_first = self.binary_to_dna_custom(barcode_first)
                actual_dna_barcode_length = len(dna_barcode_first)
                # Convert DNA barcode back to binary to get the binary barcode length for decoding
                binary_barcode_decoded = self.dna_to_binary_custom(dna_barcode_first)
                actual_binary_barcode_length = len(binary_barcode_decoded)
                # Update params with the actual binary barcode length for decoding
                self.params.dna_barcode_length = actual_binary_barcode_length
                # Store the original for use in LT code encoding (which happens with binary, not DNA)
                self.params.original_binary_barcode_length = original_binary_barcode_length
            else:
                actual_dna_barcode_length = self.params.dna_barcode_length
                actual_binary_barcode_length = self.params.dna_barcode_length
                self.params.original_binary_barcode_length = self.params.dna_barcode_length

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of DNA codewords: {len(dna_codewords)}")
                if not all(len(codeword) == len(dna_codewords[0]) for codeword in dna_codewords):
                    raise ValueError("SANITY CHECK: Not all DNA codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each DNA codeword: {len(dna_codewords[0])}")   

            # create info return dictionary
            codeword_lengths = [len(codeword) for codeword in dna_codewords]
            info = {
                "number_of_codewords": len(dna_codewords),
                "min_codeword_length": min(codeword_lengths) if codeword_lengths else 0,
                "max_codeword_length": max(codeword_lengths) if codeword_lengths else 0,
                "barcode_length": self.params.dna_barcode_length
                }

        except Exception as e:
            self.logger.error(f'An error occurred during encoding: {e}')
            self.logger.error(traceback.format_exc())

            dna_codewords = None
            info = {}

        return dna_codewords, info

    def decode(self, data):
        """
        Decodes the provided data using nohomopolymer decoding.
        Wraps the standalone decode function to use self.params.
        """
        return self._decode_function(data, self.params, self.logger)

    def process(self, data):
        """
        Processes the provided data using nohomopolymer processing.
        Wraps the standalone process function to use self.params.
        """
        return self._process_function(data, self.params, self.logger)

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

    def calculate_codeword_parameters(self, params):
        # Step 1: calculate the number of bits required to store the number of added zeros at the end of the last codeword
        zfill_bits = m.ceil(m.log2(params.codeword_length))

        if hasattr(params, 'outer_error_correction') and params.outer_error_correction == 'reedsolomon':
            # Calculate the adjustment needed to make zfill_bits a multiple of 8
            adjustment = (8 - zfill_bits % 8) % 8
            params.dna_barcode_length -= adjustment
            zfill_bits += adjustment
            
        params.zfill_bits = zfill_bits

        # Step 2: Calculate the message length
        if hasattr(params, 'inner_error_correction') and params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits - params.index_carry_length - params.ltcode_header
        else:
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits

        if hasattr(params, 'outer_error_correction') and params.outer_error_correction == 'reedsolomon':
            while message_length % 8 != 0:
                message_length -= 1
                params.dna_barcode_length += 1   

        params.message_length = message_length

        # Step 3: calculate the bits per codeword
        bits_per_codeword = message_length
        
        if hasattr(params, 'outer_error_correction') and params.outer_error_correction == 'reedsolomon':
            
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
        # if hasattr(self.params, 'outer_error_correction') and self.params.outer_error_correction == 'reedsolomon':
        #     dna_barcode_length = params.codeword_length - params.message_length
        #     params.dna_barcode_length = dna_barcode_length

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

        if hasattr(params, 'outer_error_correction') and params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCodeSynthesis(binary_codewords, params.bits_per_ec)
 
         # Step 3: apply inner error correction
        if hasattr(params, 'inner_error_correction') and params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcodesynth(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, 1)

        return binary_codewords

def attributes(inputparams):

    encoding_method = inputparams.encoding_method
    assembly_structure = 'synthesis'
    
    codeword_length = check_parameter(parameter="codeword_length",
                                      default=500,
                                      min=20,
                                      max=2000,
                                      inputparams=inputparams)

    dna_barcode_length = check_parameter(parameter="dna_barcode_length",
                                         default=m.ceil(codeword_length * 0.15),
                                         min=2,
                                         max=codeword_length - 2,
                                         inputparams=inputparams)

    codeword_maxlength_positions = check_parameter(parameter="codeword_maxlength_positions",
                                                   default=m.ceil(codeword_length * 0.15),
                                                   min=2,
                                                   max=codeword_length - dna_barcode_length,
                                                   inputparams=inputparams)

    checker = codeword_length - codeword_maxlength_positions - dna_barcode_length

    # parameter group: inner_error_correction
    if getattr(inputparams, 'inner_error_correction', None) == 'ltcode':

        index_carry_length = check_parameter(parameter="index_carry_length",
                                             default=m.ceil(codeword_length * 0.15),
                                             min=1,
                                             max=0.9 * codeword_length,
                                             inputparams=inputparams)

        ltcode_header = check_parameter(parameter="ltcode_header",
                                        default=m.ceil(codeword_length * 0.15),
                                        min=1,
                                        max=0.3 * codeword_length,
                                        inputparams=inputparams)

        percent_of_symbols = check_parameter(parameter="percent_of_symbols",
                                             default=2,
                                             min=1,
                                             max=5,
                                             inputparams=inputparams)

        checker = codeword_length - codeword_maxlength_positions - dna_barcode_length - index_carry_length - ltcode_header

    elif getattr(inputparams, 'inner_error_correction', None) is None:
        pass
    else: 
        raise ValueError("Invalid inner_error_correction method")


    # parameter group: outer_error_correction
    if getattr(inputparams, 'outer_error_correction', None) == 'reedsolomon':

        reed_solo_percentage = check_parameter(parameter="reed_solo_percentage",
                                              default=0.8,
                                              min=0.2,
                                              max=0.95,
                                              inputparams=inputparams)

    elif getattr(inputparams, 'outer_error_correction', None) is None:
        pass
    else:
        raise ValueError("Invalid outer_error_correction method")

    # TODO: does the outer_error_correction affect the checker calculation?
    
    if checker < 0:
        raise ValueError("The codeword length is too small for the given parameters.")
    
    # Build return dictionary with only relevant parameters
    attributes = {
        "encoding_method": encoding_method, 
        "assembly_structure": assembly_structure,
        "codeword_length": codeword_length, 
        "dna_barcode_length": dna_barcode_length, 
        "codeword_maxlength_positions": codeword_maxlength_positions,
    }
    
    # Add inner error correction parameters only if ltcode is used
    if hasattr(inputparams, 'inner_error_correction') and inputparams.inner_error_correction == 'ltcode':
        attributes.update({
            "inner_error_correction": 'ltcode',
            "percent_of_symbols": percent_of_symbols,
            "index_carry_length": index_carry_length,
            "ltcode_header": ltcode_header,
        })
    
    # Add outer error correction parameters only if reedsolomon is used
    if hasattr(inputparams, 'outer_error_correction') and inputparams.outer_error_correction == 'reedsolomon':
        attributes.update({
            "outer_error_correction": 'reedsolomon',
            "reed_solo_percentage": reed_solo_percentage,
        })
    
    return attributes

