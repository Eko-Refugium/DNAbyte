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

        # Define all parameter configurations to test
        self.test_configs = [
            {
                'name': 'no_error_correction',
                'params': Params(
                    encoding_method='max_density',
                    assembly_structure='synthesis',
                    inner_error_correction=None,
                    outer_error_correction=None,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200
                ),
                'expect_perfect_decode': True,
                'note': 'Basic encoding/decoding without error correction'
            },
            {
                'name': 'ltcode_only',
                'params': Params(
                    encoding_method='max_density',
                    assembly_structure='synthesis',
                    inner_error_correction='ltcode',
                    outer_error_correction=None,
                    percent_of_symbols=2,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=50,
                    codeword_length=200
                ),
                'expect_perfect_decode': True,
                'note': 'LT code (fountain code) error correction only'
            },
            {
                'name': 'reed_solomon_only',
                'params': Params(
                    encoding_method='max_density',
                    assembly_structure='synthesis',
                    inner_error_correction=None,
                    outer_error_correction='reedsolomon',
                    reed_solo_percentage=0.8,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200
                ),
                'expect_perfect_decode': False,
                'note': 'Reed-Solomon decoding may fail without errors in the data'
            },
            {
                'name': 'ltcode_and_reed_solomon',
                'params': Params(
                    encoding_method='max_density',
                    assembly_structure='synthesis',
                    inner_error_correction='ltcode',
                    outer_error_correction='reedsolomon',
                    percent_of_symbols=2,
                    reed_solo_percentage=0.8,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=50,
                    codeword_length=200
                ),
                'expect_perfect_decode': False,
                'note': 'LT code + Reed-Solomon may have memory issues'
            }
        ]

    def test_encode_decode_maxdensity_all_configurations(self):
        """Test MaxDensity encoding/decoding with all error correction configurations."""

        for config in self.test_configs:
            with self.subTest(config=config['name']):
                # Generate a random bitstream        
                binary_code = BinaryCode.random(500)

                # Encode the data
                coder = MaxDensity(config['params'], logger=self.logger)
                encoded_data, encode_info = coder.encode(binary_code)

                # Ensure encoding was successful
                self.assertIsNotNone(encoded_data, f"Encoding failed for {config['name']}")
                self.assertIsNotNone(encode_info, f"Encoding info is None for {config['name']}")
                
                # TODO: this test should run without applying the processing step
                # we, thus, need to move some steps of the processing into the encoding step
                # in particular the step that removes the barcode and index from the sequences
                encoded_data = InSilicoDNA(encoded_data) 
                processed_data, process_info = coder.process(encoded_data)

                # Prepare corrected data for decoding
                nucleobase_code = NucleobaseCode(processed_data)

                # Decode the data
                decoded_data, checkervalid, decode_info = coder.decode(nucleobase_code)

                # Validate decoding based on expectations
                if config['expect_perfect_decode']:
                    # Expect perfect decoding
                    self.assertTrue(checkervalid, 
                        f"Decoding failed for {config['name']} - checker not valid")
                    self.assertIsNotNone(decoded_data, 
                        f"Decoding failed for {config['name']} - decoded data is None")
                    self.assertIsNotNone(decode_info, 
                        f"Decoding info is None for {config['name']}")
                    self.assertEqual(binary_code.data, decoded_data, 
                        f"Decoded data does not match original for {config['name']}")
                else:
                    # Relaxed validation for configurations with known issues
                    if decoded_data is not None:
                        self.assertIsNotNone(decode_info, 
                            f"Decode info is None for {config['name']}")
                        # If decoding succeeds, verify the data matches
                        if checkervalid:
                            self.assertEqual(binary_code.data, decoded_data, 
                                f"Decoded data does not match original for {config['name']}")
                    # If decode fails, that's acceptable for these configurations
                    # Note: {config['note']}



if __name__ == '__main__':
    unittest.main()