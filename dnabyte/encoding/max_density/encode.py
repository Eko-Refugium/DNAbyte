import math as m
import traceback

from dnabyte.encode import Encode
from dnabyte.error_correction.auxiliary import MakeReedSolomonCodeSynthesis, makeltcodesynth
from dnabyte.encoding.auxiliary import create_counter_list, check_parameter, check_library
from dnabyte.encoding.max_density.decode import decode as decode_function
from dnabyte.encoding.max_density.process import process as process_function

class MaxDensity(Encode):
    """
    This class provides maxdensity encoding for DNA sequences.

    This encoding scheme maps '00' to 'A', '01' to 'G', '10' to 'C', and '11' to 'T'. It has the theoretical maximum density 
    of 2 bits per nucleotide and is, thus, mostly used for theoretical purposes.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
        self._decode_function = decode_function
        self._process_function = process_function

    def encode(self, data):
        """
        This method encodes the provided data into DNA sequences using the specified parameters.
        """
        try:

            # Calculate the codeword parameters and save them in the params object
            self.params = self.calculate_codeword_parameters(self.params)

            # Creates binary codewords from the provided data given the specified parameters
            binary_codewords = self.create_binary_codewords(data, self.params)

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of binary codewords: {len(binary_codewords)}")
                if not all(len(codeword) == len(binary_codewords[0]) for codeword in binary_codewords):
                    raise ValueError("SANITY CHECK: Not all binary codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each binary codeword: {len(binary_codewords[0])}")

            # add the barcodes to the binary blocks
            dna_codewords = []

            for i in range(len(binary_codewords)):
                # create a barcode for each codeword on base 4 (equivalent to DNA bases)
                barcodes_base4 = create_counter_list(self.params.dna_barcode_length, 4, i)
                # create a barcode for each codeword on base 2
                barcodes_base2 = self.base4_list_to_bitstring(barcodes_base4)
                # translate the barcode to DNA
                dna_codewords.append(self.binary_to_dna(barcodes_base2 + binary_codewords[i]))

            # Sanity checks
            if self.params.debug:
                self.logger.info(f"SANITY CHECK: Number of DNA codewords: {len(dna_codewords)}")
                if not all(len(codeword) == len(dna_codewords[0]) for codeword in dna_codewords):
                    raise ValueError("SANITY CHECK: Not all DNA codewords are of equal length.")
                else:
                    self.logger.info(f"SANITY CHECK: Length of each DNA codeword: {len(dna_codewords[0])}")   

            # create info return dictionary
            info = {
                "number_of_codewords": len(binary_codewords),
                "length_of_each_codeword": len(binary_codewords[0]),
                "barcode_length": self.params.dna_barcode_length,
                "data_length": len(data.data)
            }

        except Exception as e:

            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())

            dna_codewords = None
            info = {}

        return dna_codewords, info

    def decode(self, data):
        """
        Decodes the provided data using maxdensity decoding.
        Wraps the standalone decode function to use self.params.
        """
        return self._decode_function(data, self.params, self.logger)

    def process(self, data):
        """
        Processes the provided data using maxdensity processing.
        Wraps the standalone process function to use self.params.
        """
        return self._process_function(data, self.params, self.logger)

    def calculate_codeword_parameters(self, params):
        # Step 1: calculate the number of bits required to store the number of added zeros at the end of the last codeword
        # number of bits reserved to save the number of added zeros at the end of the last (and any) codeword. At the beginning of
        # each codeword, the number of zeros is stored in zfill_bits bits. This is necessary to remove the padding zeros at the end 
        # of the last codeword during decoding.

        zfill_bits = m.ceil(m.log2(params.codeword_length))

        if params.outer_error_correction == 'reedsolomon':
            # to apply the reed-solomon error correction efficiently, the dna barcode must have a multiple of 8 bits
            adjustment = (4 - zfill_bits % 4) % 4
            params.dna_barcode_length -= adjustment
            zfill_bits += adjustment

        params.zfill_bits = zfill_bits

      # Step 2: Calculate the message length

        if params.inner_error_correction == 'ltcode':
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits - params.index_carry_length - params.ltcode_header
        else:
            message_length = params.codeword_length - params.dna_barcode_length - zfill_bits

        # adjust the bits per codeword such that it is a multiple of 4
        if params.outer_error_correction == 'reedsolomon':
            while message_length % 4 != 0:
                message_length -= 1
                params.dna_barcode_length += 1
        
        params.message_length = message_length

        # Step 3: calculate bits per codeword
        bits_per_codeword = message_length * 2

        if params.outer_error_correction == 'reedsolomon':
            bits_per_codeword = m.floor(bits_per_codeword * params.reed_solo_percentage)
            bits_per_ec = (message_length * 2 - bits_per_codeword)

            # Adjust bits_per_codeword to be a multiple of 8
            adjustment = (8 - bits_per_codeword % 8) % 8
            bits_per_codeword += adjustment
            bits_per_ec -= adjustment
            params.bits_per_ec = bits_per_ec

        params.bits_per_codeword = bits_per_codeword

        # Step 4: Determine the length of the DNA barcode
        # if self.params.outer_error_correction == 'reedsolomon':
        #     # TODO: check if this is correct
        #     #dna_barcode_length = self.params.codeword_length - len(binary_codewords[0]) // 2
        #     dna_barcode_length = params.codeword_length - params.message_length // 2

        #     params.dna_barcode_length = dna_barcode_length

        return params

    def create_binary_codewords(self, data, params):
        """
        Generates binary codewords from the provided data. It supports both Reed-Solomon and LT code error correction schemes.
        """
        # Step 1: cut the binary string to substrings of length bits_per_codeword 
        binary_blocks = []
        for i in range(0, len(data.data), params.bits_per_codeword):
            binary_blocks.append(data.data[i:i + params.bits_per_codeword])

        number_of_zero_padding = len(binary_blocks[-1])

        binary_blocks[-1] = binary_blocks[-1].zfill(params.bits_per_codeword)
        
        # Step 2: create the binary codewords
        binary_codewords = []
        
        for block in binary_blocks:
            codeword = str(self.decimal_to_binary(0)).zfill(params.zfill_bits * 2) + block
            binary_codewords.append(codeword)

        # fill the last binary codeword with zeros

        if number_of_zero_padding == params.bits_per_codeword:
            number_of_zero_padding = 0

        binary_codewords[-1] = str(self.decimal_to_binary(number_of_zero_padding)).zfill(params.zfill_bits * 2) + binary_codewords[-1][params.zfill_bits * 2:]
        
        # Step 3: apply outer error correction
        if params.outer_error_correction == 'reedsolomon':
            binary_codewords = MakeReedSolomonCodeSynthesis(binary_codewords, params.bits_per_ec)

        # Step 4: apply inner error correction
        if params.inner_error_correction == 'ltcode':
            binary_codewords = makeltcodesynth(binary_codewords, params.percent_of_symbols, params.index_carry_length, params.ltcode_header, 2)
        
        return binary_codewords
        
    def base4_list_to_bitstring(self, base4_list):
        """
        Converts a list of base-4 numbers to a bitstring.
        """
        bit_string = ''.join(format(x, '02b') for x in base4_list)
        return bit_string
    
    def decimal_to_binary(self, n): 
        return bin(n).replace("0b", "") 

    def binary_to_dna(self, binary_sequence):
        """
        Converts a binary string to a DNA string.
        """
        dna_mapping = {
            '00': 'A',
            '01': 'G',
            '10': 'C',
            '11': 'T'
        }
        dna_sequence = ''.join(dna_mapping[binary_sequence[i:i+2]] for i in range(0, len(binary_sequence), 2))
        return dna_sequence
    

def attributes(inputparams):

    encoding_method = inputparams.encoding_method
    assembly_structure = 'synthesis'

    codeword_length = check_parameter(parameter="codeword_length",
                                      default=500,
                                      min=20,
                                      max=10000,
                                      inputparams=inputparams)
                            
    dna_barcode_length = check_parameter(parameter="dna_barcode_length",
                                        default=m.ceil(codeword_length * 0.15),
                                        min=5,
                                        max=1000,
                                        inputparams=inputparams)

    codeword_maxlength_positions = check_parameter(parameter="codeword_maxlength_positions",
                                                  default=m.ceil(codeword_length * 0.15),
                                                  min=1,
                                                  max=codeword_length - dna_barcode_length - 1,
                                                  inputparams=inputparams)

    checker = codeword_length - codeword_maxlength_positions - dna_barcode_length

    # parameter group: inner_error_correction
    if hasattr(inputparams, 'inner_error_correction') and inputparams.inner_error_correction == 'ltcode':

        index_carry_length = check_parameter(parameter="index_carry_length",
                                            default=m.ceil(codeword_length * 0.15),
                                            min=0.5 * codeword_length,
                                            max=0.9 * codeword_length,
                                            inputparams=inputparams)
        
        ltcode_header = check_parameter(parameter="ltcode_header",
                                       default=m.ceil(codeword_length * 0.15),
                                       min=0.05 * codeword_length,
                                       max=0.3 * codeword_length,
                                       inputparams=inputparams)

        percent_of_symbols = check_parameter(parameter="percent_of_symbols",
                                            default=2,
                                            min=1,
                                            max=5,
                                            inputparams=inputparams)

        checker = codeword_length - codeword_maxlength_positions - dna_barcode_length - index_carry_length - ltcode_header

    elif inputparams.inner_error_correction is None or not hasattr(inputparams, 'inner_error_correction'):
        pass
    else: 
        raise ValueError("Invalid inner_error_correction method")


    # parameter group: outer_error_correction
    if hasattr(inputparams, 'outer_error_correction') and inputparams.outer_error_correction == 'reedsolomon':

        reed_solo_percentage = check_parameter(parameter="reed_solo_percentage",
                                              default=0.8,
                                              min=0.5,
                                              max=0.95,
                                              inputparams=inputparams)

    elif inputparams.outer_error_correction is None or not hasattr(inputparams, 'outer_error_correction'):
        pass
    else:
        raise ValueError("Invalid outer_error_correction method")

    # TODO: does the outer_error_correction affect the checker calculation?

    # Check correct combination of codeword parameters
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
