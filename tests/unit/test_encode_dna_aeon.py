import unittest
import logging

from dnabyte import BinaryCode, NucleobaseCode, InSilicoDNA
from dnabyte.encoding.dna_aeon.encode import DNAAeon
from dnabyte.params import Params


class TestDNAAeonEncodingDecoding(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        # DNA-Aeon uses the NOREC4DNA RU10 Raptor fountain code.
        # crc error correction is used to detect corrupted packets.
        # The fountain code provides erasure resilience — overhead controls
        # how many redundant packets are generated above the minimum.
        self.test_configs = [
            {
                'name': 'crc_low_overhead',
                'params': Params(
                    encoding_method='dna_aeon',
                    assembly_structure='synthesis',
                    dna_aeon_chunk_size=10,
                    dna_aeon_overhead=0.40,
                    dna_aeon_insert_header=False,
                    dna_aeon_error_correction='crc',
                    dna_aeon_use_dna_rules=True,
                    dna_aeon_drop_upper_bound=0.5,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'DNA-Aeon with CRC and 40% overhead'
            },
            {
                'name': 'crc_high_overhead',
                'params': Params(
                    encoding_method='dna_aeon',
                    assembly_structure='synthesis',
                    dna_aeon_chunk_size=10,
                    dna_aeon_overhead=0.80,
                    dna_aeon_insert_header=False,
                    dna_aeon_error_correction='crc',
                    dna_aeon_use_dna_rules=True,
                    dna_aeon_drop_upper_bound=0.5,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'DNA-Aeon with CRC and 80% overhead for higher resilience'
            },
            {
                'name': 'crc_larger_chunks',
                'params': Params(
                    encoding_method='dna_aeon',
                    assembly_structure='synthesis',
                    dna_aeon_chunk_size=20,
                    dna_aeon_overhead=0.40,
                    dna_aeon_insert_header=False,
                    dna_aeon_error_correction='crc',
                    dna_aeon_use_dna_rules=True,
                    dna_aeon_drop_upper_bound=0.5,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'DNA-Aeon with larger chunk size (20 bytes)'
            },
            {
                'name': 'nocode_no_error_correction',
                'params': Params(
                    encoding_method='dna_aeon',
                    assembly_structure='synthesis',
                    dna_aeon_chunk_size=10,
                    dna_aeon_overhead=0.40,
                    dna_aeon_insert_header=False,
                    dna_aeon_error_correction='nocode',
                    dna_aeon_use_dna_rules=True,
                    dna_aeon_drop_upper_bound=0.5,
                    dna_barcode_length=10,
                    codeword_maxlength_positions=100,
                    codeword_length=200,
                ),
                'expect_perfect_decode': True,
                'note': 'DNA-Aeon with no packet-level error correction (nocode)'
            },
        ]

    def test_encode_decode_dna_aeon_all_configurations(self):
        """Test DNA-Aeon encoding/decoding with all configurations."""

        for config in self.test_configs:
            with self.subTest(config=config['name']):
                # Generate a random bitstream
                binary_code = BinaryCode.random(500)

                # Encode the data
                coder = DNAAeon(config['params'], logger=self.logger)
                encoded_data, encode_info = coder.encode(binary_code)

                self.assertIsNotNone(encoded_data, f"Encoding failed for {config['name']}")
                self.assertIsNotNone(encode_info, f"Encoding info is None for {config['name']}")

                # Verify encoded data is non-empty
                self.assertGreater(len(encoded_data), 0,
                    f"No sequences encoded for {config['name']}")

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
