from dnabyte.data import RawData, EncodedData, CorrectedData, DecodedData, SequencedData
from dnabyte.library import Library

class Encode:
    """
    This class is responsible for encoding the raw data into sequences of oligos.

    There are three steps involved in encoding:
    1. Splitting raw data into codewords. This step also includes error correction and barcoding.
    2. Encoding the codewords into sequences of oligos. This is called the blueprint.
    3. Subsetting the blueprint into pooling sequences. 
    
    """

    def __init__(self, params, logger=None):
        
        self.assembly_structure = params.assembly_structure
        self.encoding_scheme = params.encoding_scheme
        self.inner_error_correction = params.inner_error_correction
        self.outer_error_correction = params.outer_error_correction
        self.dna_barcode_length = params.dna_barcode_length
        self.codeword_maxlength_positions = params.codeword_maxlength_positions
        self.sigmaamount = params.sigma_amount
        self.codewordlength = params.codeword_length
        self.percentofsymbols = params.percent_of_symbols
        self.indexcarrylength = params.index_carry_length
        self.reedsolopercentage = params.reed_solo_percentage
        self.library_name = params.library_name
        self.library = params.library
        self.theory = params.theory
        self.logger = logger
        self.params = params


    def encode(self, raw_data):
        """
        Encodes the raw data based on the library structure and encoding scheme.

        This method encodes the provided raw data using the specified encoding scheme
        and library structure. The encoding scheme and library structure determine
        which encoding function is used.

        Parameters:
        raw_data (RawData): The raw data to be encoded. Must be an instance of RawData.

        Returns:
        EncodedData: An EncodedData object containing the encoded data.

        Raises:
        TypeError: If raw_data is not an instance of RawData.
        ValueError: If the encoding scheme or library structure is invalid.

        Encoding Schemes:
        - 'linear_assembly' structure with 'linear_encoding' encoding: Uses Encode_chain_simp.
        - 'linear_assembly' structure with 'binomial_encoding' encoding: Uses Encode_binom_simp.
        - 'positional_assembly' structure with 'linear_encoding' encoding: Uses Encode_chain_poly.
        - 'positional_assembly' structure with 'binomial_encoding' encoding: Uses Encode_binom_poly.
        """

        from dnabyte.encode_maxdensity import EncodeMaxDensity
        from dnabyte.encode_nohomopolyoddeven import EncodeNoHomoPolyOddEven
        from dnabyte.encode_linear_chain import EncodeLinearChain
        from dnabyte.encode_linear_binom import EncodeLinearBinom
        from dnabyte.encode_poly_chain import EncodePolyChain
        from dnabyte.encode_poly_binom import EncodePolyBinom


        if not isinstance(raw_data, RawData):
            raise TypeError("raw_data must be an instance of RawData")

        if self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'linear_encoding':
            enc = EncodeLinearChain(self.params, logger=self.logger)
            encoded_data, info = enc.encode_linear_chain(raw_data)

        elif self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'binomial_encoding':
            enc = EncodeLinearBinom(self.params, logger=self.logger)
            encoded_data, info = enc.encode_linear_binom(raw_data)

        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'linear_encoding':
            enc = EncodePolyChain(self.params, logger=self.logger)
            encoded_data, info = enc.encode_chain_poly(raw_data)

        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'binomial_encoding':
            enc = EncodePolyBinom(self.params, logger=self.logger)
            encoded_data, info = enc.encode_binom_poly(raw_data)

        elif self.assembly_structure == 'synthesis' and self.encoding_scheme == 'max_density_encoding':
            enc = EncodeMaxDensity(self.params, logger=self.logger)
            encoded_data, info = enc.encode_max_density(raw_data)

        elif self.assembly_structure == 'synthesis' and self.encoding_scheme == 'no_homopolymeroddeven_encoding':
            enc = EncodeNoHomoPolyOddEven(self.params, logger=self.logger)
            encoded_data, info = enc.encode_nohomopoly_odd_even(raw_data)

        else:
            raise ValueError("Invalid encoding scheme or library structure")

        obj = EncodedData(encoded_data)
        obj.file_paths = raw_data.file_paths

        return obj, info


    def decode(self, corrected_data):

        """
        Decodes the corrected data based on the assembly structure and encoding scheme.

        This method decodes the provided corrected data using the specified assembly structure
        and encoding scheme. The assembly structure and encoding scheme determine
        which decoding function is used.

        Parameters:
        corrected_data (CorrectedData): The corrected data to be decoded. Must be an instance of CorrectedData.

        Returns:
        DecodedData: A DecodedData object containing the decoded data.

        Raises:
        TypeError: If corrected_data is not an instance of CorrectedData.
        ValueError: If the encoding scheme or library structure is invalid.

        Assembly Structures:
        - 'linear_assembly' with 'linear_encoding' encoding: Uses Decode_simple_chain.
        - 'linear_assembly' with 'binomial_encoding' encoding: Uses Decode_simple_binom.
        - 'positional_assembly' with 'linear_encoding' encoding: Uses Decode_poly_chain.
        - 'positional_assembly' with 'binomial_encoding' encoding: Uses Decode_poly_binom.
        """

        from dnabyte.decode_linear_chain import DecodeLinearChain
        from dnabyte.decode_linear_binom import DecodeLinearBinom
        from dnabyte.decode_poly_chain import DecodePolyChain
        from dnabyte.decode_poly_binom import DecodePolyBinom
        from dnabyte.decode_maxdensity import DecodeMaxDensity
        from dnabyte.decode_nohomopolyoddeven import DecodeNoHomoPolyOddEven

        if not isinstance(corrected_data, CorrectedData):
            raise TypeError("corrected_data must be an instance of CorrectedData")

        if self.params.assembly_structure == 'linear_assembly' and self.params.encoding_scheme == 'linear_encoding':
            dec = DecodeLinearChain(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_linear_chain(corrected_data)

        elif self.params.assembly_structure == 'linear_assembly' and self.params.encoding_scheme == 'binomial_encoding':
            dec = DecodeLinearBinom(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_linear_binom(corrected_data)

        elif self.params.assembly_structure == 'positional_assembly' and self.params.encoding_scheme == 'linear_encoding':
            dec = DecodePolyChain(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_poly_chain(corrected_data)

        elif self.params.assembly_structure == 'positional_assembly' and self.params.encoding_scheme == 'binomial_encoding':
            dec = DecodePolyBinom(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_poly_binom(corrected_data)

        elif self.params.assembly_structure == 'synthesis' and self.params.encoding_scheme == 'max_density_encoding':
            dec = DecodeMaxDensity(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_maxdensity(corrected_data)

        elif self.params.assembly_structure == 'synthesis' and self.params.encoding_scheme == 'no_homopolymeroddeven_encoding':
            dec = DecodeNoHomoPolyOddEven(self.params, logger=self.logger)
            decoded_data, valid, info = dec.decode_nohomopolyoddeven(corrected_data)

        else:
            raise ValueError("Invalid encoding scheme or library structure")

        obj = DecodedData(corrected_data)
        obj.data = decoded_data

        return obj, valid, info
    

    def process(self, data):
        
        """
        Processes the sequenced data based on the assembly structure and encoding scheme.

        This method processes the provided sequenced data using the specified assembly structure
        and encoding scheme. The assembly structure and encoding scheme determine
        which processing function is used.

        Parameters:
        data (SequencedData): The sequenced data to be processed. Must be an instance of SequencedData.

        Returns:
        CorrectedData: A CorrectedData object containing the processed data.

        Raises:
        TypeError: If data is not an instance of SequencedData.
        ValueError: If the encoding scheme or library structure is invalid.
        
        Assembly Structures:
        - 'linear_assembly' with 'linear_encoding' encoding: Uses Process.
        - 'linear_assembly' with 'binomial_encoding' encoding: Uses Process.
        - 'positional_assembly' with 'linear_encoding' encoding: Uses Process.
        - 'positional_assembly' with 'binomial_encoding' encoding: Uses Process.
        """

        from dnabyte.process_maxdensity import ProcessMaxDensity
        from dnabyte.process_nohomopolyoddeven import ProcessNoHomoPolyOddEven
        from dnabyte.process_linear_chain import ProcessLinearChain
        from dnabyte.process_linear_binom import ProcessLinearBinom
        from dnabyte.process_poly_chain import ProcessPolyChain
        from dnabyte.process_poly_binom import ProcessPolyBinom

        if data is not None and isinstance(data, SequencedData):
                    
            if self.assembly_structure == 'synthesis' and self.encoding_scheme == 'max_density_encoding':
                pro = ProcessMaxDensity(self.params, logger=self.logger)
                processed_data, info = pro.process_maxdensity(data)

            if self.assembly_structure == 'synthesis' and self.encoding_scheme == 'no_homopolymeroddeven_encoding':
                pro = ProcessNoHomoPolyOddEven(self.params, logger=self.logger)
                processed_data, info = pro.process_nohomopolyoddeven(data)

            if self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'linear_encoding':
                pro = ProcessLinearChain(self.params, logger=self.logger)
                processed_data, info = pro.process_linear_chain(data)

            elif self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'binomial_encoding':
                pro = ProcessLinearBinom(self.params, logger=self.logger)
                processed_data, info = pro.process_linear_binom(data)
            
            elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'linear_encoding':
                pro = ProcessPolyChain(self.params, logger=self.logger)
                processed_data, info = pro.process_poly_chain(data)
                
            elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'binomial_encoding':
                pro = ProcessPolyBinom(self.params, logger=self.logger)
                processed_data, info = pro.process_poly_binom(data)
        
        obj = CorrectedData(data)
        obj.data = processed_data

        return obj, info