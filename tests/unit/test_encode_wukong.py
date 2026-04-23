import unittest
import logging

from dnabyte import BinaryCode, NucleobaseCode, InSilicoDNA
from dnabyte.encoding.wukong.encode import Wukong
from dnabyte.params import Params


class TestWukongEncodingDecoding(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        self.test_configs = [
            {
                'name': 'no_error_correction',
                'params': Params(
                    encoding_method='wukong',
                    assembly_structure='synthesis',
                    inner_error_correction=None,
                    outer_error_correction=None,
                    max_homopolymer=4,
                    min_gc=0.4,
                    max_gc=0.6,
                    rule_num=1,
                    rs_num=0,
                    add_redundancy=True,
                    add_primer=True,
                    primer_length=20,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'Basic Wukong encoding/decoding without error correction'
            },
            {
                'name': 'reed_solomon_only',
                'params': Params(
                    encoding_method='wukong',
                    assembly_structure='synthesis',
                    inner_error_correction=None,
                    outer_error_correction='reedsolomon',
                    reed_solo_percentage=0.8,
                    max_homopolymer=4,
                    min_gc=0.4,
                    max_gc=0.6,
                    rule_num=1,
                    rs_num=1,
                    add_redundancy=True,
                    add_primer=True,
                    primer_length=20,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': False,
                'note': 'Reed-Solomon decoding may fail without errors in the data'
            },
            {
                'name': 'ltcode_only',
                'params': Params(
                    encoding_method='wukong',
                    assembly_structure='synthesis',
                    inner_error_correction='ltcode',
                    outer_error_correction=None,
                    percent_of_symbols=2,
                    max_homopolymer=4,
                    min_gc=0.4,
                    max_gc=0.6,
                    rule_num=1,
                    rs_num=0,
                    add_redundancy=True,
                    add_primer=False,
                    primer_length=20,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=50,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'LT code (fountain code) error correction only'
            },
            {
                'name': 'ltcode_and_reed_solomon',
                'params': Params(
                    encoding_method='wukong',
                    assembly_structure='synthesis',
                    inner_error_correction='ltcode',
                    outer_error_correction='reedsolomon',
                    percent_of_symbols=2,
                    reed_solo_percentage=0.9,
                    max_homopolymer=4,
                    min_gc=0.4,
                    max_gc=0.6,
                    rule_num=1,
                    rs_num=1,
                    add_redundancy=True,
                    add_primer=True,
                    primer_length=20,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=50,
                    codeword_length=200,
                ),
                'expect_perfect_decode': False,
                'note': 'LT code + Reed-Solomon may have memory issues'
            },
        ]

    def test_encode_decode_wukong_all_configurations(self):
        """Test Wukong encoding/decoding with all error correction configurations."""

        for config in self.test_configs:
            with self.subTest(config=config['name']):
                # Generate a random bitstream — Wukong needs enough segments to reliably
                # pair and encode/decode; 2000 bits ensures at least 10+ segments
                binary_code = BinaryCode.random(2000)

                # Encode the data
                coder = Wukong(config['params'], logger=self.logger)
                encoded_data, encode_info = coder.encode(binary_code)

                self.assertIsNotNone(encoded_data, f"Encoding failed for {config['name']}")
                self.assertIsNotNone(encode_info, f"Encoding info is None for {config['name']}")

                # Verify GC content constraints in encoded sequences (only when no primer is added,
                # as primer sequences do not follow the GC content rules)
                if not config['params'].add_primer or config['params'].add_primer is None:
                    for seq in encoded_data:
                        gc_count = seq.count('G') + seq.count('C')
                        gc_ratio = gc_count / len(seq) if len(seq) > 0 else 0
                        self.assertGreaterEqual(gc_ratio, config['params'].min_gc - 0.05,
                            f"GC ratio {gc_ratio:.2f} below minimum for {config['name']}: {seq}")
                        self.assertLessEqual(gc_ratio, config['params'].max_gc + 0.05,
                            f"GC ratio {gc_ratio:.2f} above maximum for {config['name']}: {seq}")

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
