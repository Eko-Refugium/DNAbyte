from pyxdameraulevenshtein import damerau_levenshtein_distance
import random
import math as m
from collections import defaultdict, Counter
from typing import List, Dict, Tuple

from dnabyte.data import RawData, EncodedData, CorrectedData, DecodedData, SequencedData
from dnabyte.oligo import translate_nested_list, create_positional_libraries
from dnabyte.auxiliary import complementmap, transpose_lists, transpose_matrix, split_string_into_chunks, dna_to_base3_no_homopolymers, binary_to_dna_no_homopolymers
from dnabyte.encode import Encode
from dnabyte.data import CorrectedData

from dnabyte.process_maxdensity import ProcessMaxDensity
from dnabyte.process_nohomopolyoddeven import ProcessNoHomoPolyOddEven
from dnabyte.process_linear_chain import ProcessLinearChain
from dnabyte.process_linear_binom import ProcessLinearBinom
from dnabyte.process_poly_chain import ProcessPolyChain
from dnabyte.process_poly_binom import ProcessPolyBinom


class Process:

    """
    This class processes the sequenced data to data that can be decoded.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
        self.assembly_structure = params.assembly_structure
        

    def process(self, data):

        if not isinstance(data, SequencedData):
            raise TypeError("raw_data must be an instance of RawData")
        
        if self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'linear_encoding':
            pro = ProcessLinearChain(self.params, logger=self.logger)
            processed_data, info = pro.process_linear_chain(data)

        elif self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'binomial_encoding':
            pro = ProcessLinearBinom(self.params, logger=self.logger)
            processed_data, info = pro.encode_linear_binom(data)

        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'linear_encoding':
            pro = ProcessPolyChain(self.params, logger=self.logger)
            processed_data, info = pro.encode_chain_poly(data)

        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'binomial_encoding':
            pro = ProcessPolyBinom(self.params, logger=self.logger)
            processed_data, info = pro.encode_binom_poly(data)

        elif self.assembly_structure == 'synthesis' and self.encoding_scheme == 'max_density_encoding':
            pro = ProcessMaxDensity(self.params, logger=self.logger)
            processed_data, self.dna_barcode_length_temp, info = pro.encode_max_density(data)

        elif self.assembly_structure == 'synthesis' and self.encoding_scheme == 'no_homopolymeroddeven_encoding':
            pro = ProcessNoHomoPolyOddEven(self.params, logger=self.logger)
            processed_data, self.dna_barcode_length_temp, info = pro.encode_nohomopoly_odd_even(data)

        else:
            raise ValueError("Invalid encoding scheme or library structure")

            obj = EncodedData(processed_data)
            obj.file_paths = raw_data.file_paths

            return obj, info
