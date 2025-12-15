from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.library import Library
import importlib

class Encode:
    """
    This class is responsible for encoding the raw data into sequences of oligos.

    There are three steps involved in encoding:
    1. Splitting raw data into codewords. This step also includes error correction and barcoding.
    2. Encoding the codewords into sequences of oligos. 
    3. Subsetting the blueprint into pooling sequences. 
    
    """

    def __init__(self, params, logger=None):
        
        self.encoding_method = params.encoding_method
        self.logger = logger
        self.params = params
        self.encoding_plugins = params.encoding_plugins
        self.debug = params.debug

    def encode(self, data):
        """
        Encodes the raw data based on the library structure and encoding method.

        This method encodes the provided raw data using the specified encoding method
        and library structure. The encoding method and library structure determine
        which encoding function is used.

        Parameters:
        data (BinaryCode): The raw data to be encoded. Must be an instance of BinaryCode.

        Returns:
        NucleobaseCode: A NucleobaseCode object containing the encoded data.

        Raises:
        TypeError: If data is not an instance of BinaryCode.
        ValueError: If the encoding method or library structure is invalid.
        """
        if isinstance(data, BinaryCode):
            try:
                encode_class = self.encoding_plugins[self.encoding_method]
                plugin = encode_class(self.params, logger=self.logger)
                encoded_data, info = plugin.encode(data)
                obj = NucleobaseCode(encoded_data)
                obj.file_paths = data.file_paths
                return obj, info
            except KeyError:
                raise ValueError(f"Encoding method '{self.encoding_method}' not found in plugins.")
        else:
            raise TypeError("data must be an instance of BinaryCode")


    def decode(self, data):

        """
        Decodes the corrected data based on the assembly structure and encoding method.

        This method decodes the provided corrected data using the specified assembly structure
        and encoding method. The assembly structure and encoding method determine
        which decoding function is used.

        Parameters:
        data (NucleobaseCode): The corrected data to be decoded. Must be an instance of NucleobaseCode.

        Returns:
        BinaryCode: A BinaryCode object containing the decoded data.

        Raises:
        TypeError: If data is not an instance of NucleobaseCode.
        ValueError: If the encoding method or library structure is invalid.
        """
        if isinstance(data, NucleobaseCode):
            try:
                encode_class = self.encoding_plugins[self.encoding_method]
                plugin = encode_class(self.params, logger=self.logger)
                #decoded_data, valid, info = plugin.decode(data=data, params=self.params, logger=self.logger)
                decoded_data, valid, info = plugin.decode(data=data)
                obj = BinaryCode(decoded_data)
                obj.file_paths = data.file_paths
                return obj, valid, info
            except KeyError:
                raise ValueError(f"Decode Module for encoding method '{self.encoding_method}' not found in plugins.")
        else:
            raise TypeError("data must be an instance of NucleobaseCode")

    def process(self, data):
        
        """
        Processes the sequenced data based on the assembly structure and encoding method.

        This method processes the provided sequenced data using the specified assembly structure
        and encoding method. The assembly structure and encoding method determine
        which processing function is used.

        Parameters:
        data (InSilicoDNA): The sequenced data to be processed. Must be an instance of InSilicoDNA.

        Returns:
        BinaryCode: A BinaryCode object containing the processed data.

        Raises:
        TypeError: If data is not an instance of InSilicoDNA.
        ValueError: If the encoding method or library structure is invalid.
        """
        if isinstance(data, InSilicoDNA):
            try:
                print(self.encoding_method, self.encoding_plugins)
                encode_class = self.encoding_plugins[self.encoding_method]
                plugin = encode_class(self.params, logger=self.logger)
                processed, info = plugin.process(data)
                obj = NucleobaseCode(processed)
                obj.file_paths = data.file_paths
                return obj, info
            except KeyError:
                raise ValueError(f"Process Module for encoding method '{self.encoding_method}' not found in plugins.")
        else:
            raise TypeError("data must be an instance of InSilicoDNA")

