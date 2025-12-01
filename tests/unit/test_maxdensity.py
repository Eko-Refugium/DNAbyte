import unittest
import random
import logging

from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.encoding.max_density.encode import MaxDensity
from dnabyte.encoding.max_density.decode import MaxDensity
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
            assembly_structure='synthesis',
            encoding_scheme='max_density_encoding',
            inner_error_correction='ltcode',
            outer_error_correction='reedsolomon',
            dna_barcode_length=10,
            codeword_maxlength_positions=100,
            codeword_length=200,
            index_carry_length=5,
            percent_of_symbols=0.1,
            reed_solo_percentage=0.2
        )

    def generate_random_bitstream(self, length):
        return ''.join(random.choice('01') for _ in range(length))

    def test_encode_decode_maxdensity(self):
        # Generate a random bitstream
        original_data = self.generate_random_bitstream(1000)
        
        binary_code = BinaryCode(original_data)

        # Encode the data
        encoder = MaxDensity(self.params, logger=self.logger)
        encoded_data, barcode, encode_info = encoder.encode_maxdensity(binary_code)

        # Ensure encoding was successful
        self.assertIsNotNone(encoded_data, "Encoding failed")
        self.assertIsNotNone(barcode, "Encoding failed")
        self.assertIsNotNone(encode_info, "Encoding failed")

        # Prepare corrected data for decoding
        nucleobase_code = NucleobaseCode(encoded_data)

        # Decode the data
        decoder = MaxDensity(self.params, logger=self.logger)
        decoded_data, checkervalid, decode_info = decoder.decode_maxdensity(nucleobase_code)

        # Ensure decoding was successful
        self.assertTrue(checkervalid, "Decoding failed")
        self.assertIsNotNone(decoded_data, "Decoding failed")
        self.assertIsNotNone(decode_info, "Decoding failed")

        # Compare the original data with the decoded data
        self.assertEqual(original_data, decoded_data, "Decoded data does not match the original data")

if __name__ == '__main__':
    unittest.main()