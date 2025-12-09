import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.encoding.auxiliary import create_counter_list, decimal_to_binary
from dnabyte.error_correction.auxiliary import MakeReedSolomonCodeSynthesis, makeltcodesynth
from dnabyte.encoding.no_homopolymer.decode import decode
from dnabyte.encoding.no_homopolymer.process import process

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
                dna_codewords.append(self.binary_to_dna_custom(binary_codewords[i]))

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of DNA codewords: {len(dna_codewords)}")
                if not all(len(codeword) == len(dna_codewords[0]) for codeword in dna_codewords):
                    raise ValueError("SANITY CHECK: Not all DNA codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each DNA codeword: {len(dna_codewords[0])}")   

            # create info return dictionary
            info = {
                "number_of_codewords": len(dna_codewords),
                "length_of_each_codeword": [len(codeword) for codeword in binary_codewords],
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

        if params.outer_error_correction == 'reedsolomon':
            # Calculate the adjustment needed to make zfill_bits a multiple of 8
            adjustment = (8 - zfill_bits % 8) % 8
            params.dna_barcode_length -= adjustment
            zfill_bits += adjustment
            
        params.zfill_bits = zfill_bits

        # Step 2: Calculate the message length
        if params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - 2 * params.dna_barcode_length - zfill_bits - params.index_carry_length
        else:
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits

        if params.outer_error_correction == 'reedsolomon':
            while message_length % 8 != 0:
                message_length -= 1
                params.dna_barcode_length += 1   

        params.message_length = message_length

        # Step 3: calculate the bits per codeword
        bits_per_codeword = message_length
        
        if params.outer_error_correction == 'reedsolomon':
            
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
        if self.params.outer_error_correction == 'reedsolomon':
            #barcode = self.params.codeword_length - len(binary_codewords[0])
            dna_barcode_length = params.codeword_length - params.message_length

            params.dna_barcode_length = dna_barcode_length

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

        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCodeSynthesis(binary_codewords, params.bits_per_ec)
 
         # Step 3: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcodesynth(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.dna_barcode_length, 1)

        return binary_codewords

def attributes(inputparams):

    encoding_method = 'no_homopolymeroddeven_encoding'
    assembly_structure = 'synthesis'
    
    # parameter: codeword_length
    if not hasattr(inputparams, 'codeword_length') or inputparams.codeword_length is None:
        codeword_length = 500
    # TODO: set proper limits
    elif not (20 < inputparams.codeword_length <= 2000):
        raise ValueError("codeword_length must be between 20 and 2000")
    else:
        codeword_length = inputparams.codeword_length

    # parameter: dna_barcode_length
    if not hasattr(inputparams, 'dna_barcode_length') or inputparams.dna_barcode_length is None:
        dna_barcode_length = m.ceil(codeword_length * 0.15)
    # TODO: set proper limits
    elif not (5 < inputparams.dna_barcode_length <= codeword_length // 2):
        raise ValueError("dna_barcode_length must be between 5 and codeword_length/2")
    else:
        dna_barcode_length = inputparams.dna_barcode_length

    # parameter: codeword_maxlength_positions
    if not hasattr(inputparams, 'codeword_maxlength_positions') or inputparams.codeword_maxlength_positions is None:
        codeword_maxlength_positions = m.ceil(codeword_length * 0.15)
    # TODO: set proper limits
    elif not (10 < inputparams.codeword_maxlength_positions <= codeword_length - dna_barcode_length):
        raise ValueError("codeword_maxlength_positions must be between 10 and codeword_length - dna_barcode_length")
    else:
        codeword_maxlength_positions = inputparams.codeword_maxlength_positions

    checker = codeword_length - codeword_maxlength_positions - dna_barcode_length

    # parameter group: inner_error_correction
    if hasattr(inputparams, 'inner_error_correction') and inputparams.inner_error_correction == 'ltcode':

        # parameter: index_carry_length
        if not hasattr(inputparams, 'index_carry_length') or inputparams.index_carry_length is None:
            index_carry_length = m.ceil(codeword_length * 0.15) # default
        #TODO: set proper range for index_carry_length
        elif not (0.5 * codeword_length < inputparams.index_carry_length < 0.9 * codeword_length):
            raise ValueError(f"index_carry_length must be greater than 0.5 * codeword_length and less than 0.9 * codeword_length, got {inputparams.index_carry_length}")
        else:
            index_carry_length = inputparams.index_carry_length

        # parameter: ltcode_header
        if not hasattr(inputparams, 'ltcode_header') or inputparams.ltcode_header is None:
            ltcode_header = m.ceil(codeword_length * 0.15) # default
        # TODO: set proper range for ltcode_header
        elif not (0.05 * codeword_length < inputparams.ltcode_header < 0.3 * codeword_length):
            raise ValueError(f"ltcode_header must be greater than 0.05 * codeword_length and less than 0.2 * codeword_length, got {inputparams.ltcode_header}")
        else:
            ltcode_header = inputparams.ltcode_header

        # parameter: percent_of_symbols
        if not hasattr(inputparams, 'percent_of_symbols') or inputparams.percent_of_symbols is None:
            percent_of_symbols = 2 # default
        # TODO: set proper range for percent_of_symbols
        elif not (1 < inputparams.percent_of_symbols < 5):
            raise ValueError(f"percent_of_symbols must be greater than 1 and less than 5, got {inputparams.percent_of_symbols}")
        else:
            percent_of_symbols = inputparams.percent_of_symbols

        checker = codeword_length - codeword_maxlength_positions - dna_barcode_length - index_carry_length - ltcode_header

    elif inputparams.inner_error_correction is None or not hasattr(inputparams, 'inner_error_correction'):
        pass
    else: 
        raise ValueError("Invalid inner_error_correction method")


    # parameter group: outer_error_correction
    if hasattr(inputparams, 'outer_error_correction') and inputparams.outer_error_correction == 'reedsolomon':

        # parameter: reed_solo_percentage
        if not hasattr(inputparams, 'reed_solo_percentage') or inputparams.reed_solo_percentage is None:
            reed_solo_percentage = 0.8 # default
        # TODO: set proper range for reed_solo_percentage
        elif not (0.5 < inputparams.reed_solo_percentage < 0.95):
            raise ValueError(f"reed_solo_percentage must be greater than 0.5 and less than 0.95, got {inputparams.reed_solo_percentage}")
        else:
            reed_solo_percentage = inputparams.reed_solo_percentage
    elif inputparams.outer_error_correction is None:
        pass
    else:
        raise ValueError("Invalid outer_error_correction method")

    # TODO: does the outer_error_correction affect the checker calculation?
    
    if checker < 0:
        raise ValueError("The codeword length is too small for the given parameters.")
    
    # Build return dictionary with only relevant parameters
    result = {
        "encoding_method": encoding_method, 
        "assembly_structure": assembly_structure,
        "codeword_length": codeword_length, 
        "dna_barcode_length": dna_barcode_length, 
        "codeword_maxlength_positions": codeword_maxlength_positions,
    }
    
    # Add inner error correction parameters only if ltcode is used
    if hasattr(inputparams, 'inner_error_correction') and inputparams.inner_error_correction == 'ltcode':
        result.update({
            "inner_error_correction": 'ltcode',
            "percent_of_symbols": percent_of_symbols,
            "index_carry_length": index_carry_length,
            "ltcode_header": ltcode_header,
        })
    
    # Add outer error correction parameters only if reedsolomon is used
    if hasattr(inputparams, 'outer_error_correction') and inputparams.outer_error_correction == 'reedsolomon':
        result.update({
            "outer_error_correction": 'reedsolomon',
            "reed_solo_percentage": reed_solo_percentage,
        })
    
    return result

