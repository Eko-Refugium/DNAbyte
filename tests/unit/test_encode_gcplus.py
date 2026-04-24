import unittest
import logging

from dnabyte import BinaryCode, NucleobaseCode, InSilicoDNA
from dnabyte.encoding.gcplus.encode import GCPlus
from dnabyte.params import Params


class TestGCPlusEncodingDecoding(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        # GC+ uses a systematic inner RS code + codebook parity for per-oligo
        # edit-error correction.  gcplus_k and gcplus_l define the codeword
        # structure; gcplus_c1 controls error correction capacity.
        self.test_configs = [
            {
                'name': 'c1_2_low_correction',
                'params': Params(
                    encoding_method='gcplus',
                    assembly_structure='synthesis',
                    gcplus_k=168,
                    gcplus_l=8,
                    gcplus_c1=2,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'GC+ with c1=2 (minimal correction capacity)'
            },
            {
                'name': 'c1_4_higher_correction',
                'params': Params(
                    encoding_method='gcplus',
                    assembly_structure='synthesis',
                    gcplus_k=168,
                    gcplus_l=8,
                    gcplus_c1=4,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'GC+ with c1=4 (higher correction capacity)'
            },
            {
                'name': 'c1_2_with_ltcode',
                'params': Params(
                    encoding_method='gcplus',
                    assembly_structure='synthesis',
                    inner_error_correction='ltcode',
                    outer_error_correction=None,
                    percent_of_symbols=2,
                    gcplus_k=168,
                    gcplus_l=8,
                    gcplus_c1=2,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=50,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'GC+ with LT code fountain correction'
            },
            {
                'name': 'c1_2_with_reed_solomon',
                'params': Params(
                    encoding_method='gcplus',
                    assembly_structure='synthesis',
                    inner_error_correction=None,
                    outer_error_correction='reedsolomon',
                    reed_solo_percentage=0.8,
                    gcplus_k=168,
                    gcplus_l=8,
                    gcplus_c1=2,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': False,
                'note': 'GC+ with Reed-Solomon may fail without errors in the data'
            },
        ]

    def test_encode_decode_gcplus_all_configurations(self):
        """Test GC+ encoding/decoding with all error correction configurations."""

        for config in self.test_configs:
            with self.subTest(config=config['name']):
                # Generate a random bitstream
                binary_code = BinaryCode.random(500)

                # Encode the data
                coder = GCPlus(config['params'], logger=self.logger)
                encoded_data, encode_info = coder.encode(binary_code)

                self.assertIsNotNone(encoded_data, f"Encoding failed for {config['name']}")
                self.assertIsNotNone(encode_info, f"Encoding info is None for {config['name']}")

                # Process the encoded data
                encoded_data = InSilicoDNA(encoded_data)
                processed_data, process_info = coder.process(encoded_data)

                # Prepare for decoding
                nucleobase_code = NucleobaseCode(processed_data)

                # Decode the data
                decoded_data, checkervalid, decode_info = coder.decode(nucleobase_code)

                # Validate decoding based on expectations
                if config['expect_perfect_decode']:
                    self.assertTrue(checkervalid,
                        f"Decoding failed for {config['name']} - checker not valid")
                    self.assertIsNotNone(decoded_data,
                        f"Decoding failed for {config['name']} - decoded data is None")


if __name__ == '__main__':
    unittest.main()
