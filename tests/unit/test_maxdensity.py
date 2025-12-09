import unittest
import random
import logging

from dnabyte import BinaryCode, NucleobaseCode, InSilicoDNA
from dnabyte.encoding.max_density.encode import MaxDensity
from dnabyte.params import Params

class TestMaxDensityEncodingDecoding(unittest.TestCase):

    def setUp(self):
        # Set up logging configuration
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        # Define parameters for encoding and decoding
        self.params = Params(
            encoding_method='max_density',
            assembly_structure='synthesis',
            inner_error_correction=None,
            outer_error_correction=None,
            dna_barcode_length=10,
            codeword_maxlength_positions=100,
            codeword_length=200
        )

    def test_encode_decode_maxdensity(self):

        # Generate a random bitstream        
        binary_code = BinaryCode.random(500)

        # Encode the data
        coder = MaxDensity(self.params, logger=self.logger)
        encoded_data, encode_info = coder.encode(binary_code)

        # Ensure encoding was successful
        self.assertIsNotNone(encoded_data, "Encoding failed")
        self.assertIsNotNone(encode_info, "Encoding failed")
        
        # TODO: this test should run without applying the processing step
        # we, thus, need to move some steps of the processing into the encoding step
        # in particular the step that removes the barcode and index from the sequences
        encoded_data = InSilicoDNA(encoded_data) 
        procesed_data, info = coder.process(encoded_data)

        # Prepare corrected data for decoding
        nucleobase_code = NucleobaseCode(procesed_data)

        # Decode the data
        decoded_data, checkervalid, decode_info = coder.decode(nucleobase_code)

        # Ensure decoding was successful
        self.assertTrue(checkervalid, "Decoding failed")
        self.assertIsNotNone(decoded_data, "Decoding failed")
        self.assertIsNotNone(decode_info, "Decoding failed")

        # Compare the original data with the decoded data
        self.assertEqual(binary_code.data, decoded_data, "Decoded data does not match the original data")



if __name__ == '__main__':
    unittest.main()