import tempfile
import os
import traceback

from dnabyte.encode import Encode
from dnabyte.encoding.goldman.StorageD.codec import GoldmanEncode


class Goldman(Encode):
    """
    This class provides Goldman encoding for DNA sequences.
    Uses the GoldmanEncode from dnabyte.encoding.goldman.StorageD.codec for encoding.

    Goldman encoding uses a Huffman-based ternary encoding with a rotating code
    to avoid homopolymers. It uses 4x overlapping segments for error tolerance.
    Goldman does NOT use primers, RS codes, or the standard split/merge pipeline.

    Reference:
    Goldman, N., et al. (2013). "Towards practical, high-capacity, low-maintenance
    information storage in synthesized DNA." Nature 494(7435): 77-80.
    DOI: 10.1038/nature11875
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def encode(self, data):
        """
        Encodes binary data using Goldman encoding.
        Writes binary data to temporary file, encodes using GoldmanEncode,
        deletes temporary file, and returns DNA sequences.
        """
        temp_file_path = None
        output_dir = None
        try:
            # Create temporary directory for output
            output_dir = tempfile.mkdtemp()

            # Pack bitstream into actual bytes for efficient encoding
            # data.data is a string of '0' and '1' characters
            bitstream = data.data
            # Pad to multiple of 8 if needed
            padded_len = len(bitstream)
            if padded_len % 8 != 0:
                bitstream = bitstream + '0' * (8 - padded_len % 8)

            byte_data = bytes(int(bitstream[i:i+8], 2) for i in range(0, len(bitstream), 8))

            temp_file = tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.bin')
            temp_file_path = temp_file.name
            temp_file.write(byte_data)
            temp_file.close()

            if self.logger:
                self.logger.info(f"Binary data written to temporary file: {temp_file_path}")

            # Use Goldman encoding
            # Goldman only takes input_file_path, output_dir, sequence_length
            encode_worker = GoldmanEncode(
                input_file_path=temp_file_path,
                output_dir=output_dir,
                sequence_length=getattr(self.params, 'sequence_length', 200),
            )

            # Perform encoding (Goldman overrides common_encode entirely)
            result_file = encode_worker.common_encode()

            if self.logger:
                self.logger.info(f"Encoding completed. Result file: {result_file}")

            # Read result DNA sequences and metadata from the FASTA
            # Goldman metadata format: >indexLen:{},addLen:{},fileExtension:{}
            dna_codewords = []
            original_fasta_metadata = ""

            if os.path.exists(result_file):
                with open(result_file, 'r') as f:
                    raw_content = f.read()

                lines = [line.strip() for line in raw_content.split('\n')]

                for line in lines:
                    if line.startswith('>') and 'indexLen:' in line:
                        # This is the Goldman metadata header
                        original_fasta_metadata = line
                    elif line.startswith('>'):
                        # Skip sequence headers like >seq_1
                        continue
                    else:
                        # Only add DNA sequences (non-empty, non-header)
                        if line and not line.startswith('seq_'):
                            dna_codewords.append(line)

            # Store original FASTA metadata on params so decode can use it
            self.params.goldman_fasta_metadata = original_fasta_metadata
            # Goldman metadata doesn't include totalBit, so store it separately
            self.params.goldman_total_bits = len(data.data)

            total_bits = len(data.data)

            # Create info dictionary
            info = {
                "number_of_codewords": len(dna_codewords),
                "result_file": result_file,
                "data_length": len(data.data),
                "barcode_length": 0,  # Goldman has no primers/barcodes
                "metadata": original_fasta_metadata,
                "left_primer": "",
                "right_primer": "",
                "total_bits": total_bits,
                "original_fasta_metadata": original_fasta_metadata
            }

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())

            dna_codewords = None
            info = {}

        finally:
            # Delete temporary file
            if temp_file_path and os.path.exists(temp_file_path):
                try:
                    os.remove(temp_file_path)
                    if self.logger:
                        self.logger.info(f"Temporary file deleted: {temp_file_path}")
                except Exception as e:
                    if self.logger:
                        self.logger.warning(f"Could not delete temporary file: {str(e)}")

        return dna_codewords, info

    def decode(self, data):
        """
        Decodes DNA sequences using Goldman decoding.
        """
        try:
            from dnabyte.encoding.goldman.decode import decode as decode_function
            result = decode_function(data, self.params, self.logger)
            return result
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    def process(self, data):
        """
        Processes DNA sequences using Goldman processing.
        """
        try:
            from dnabyte.encoding.goldman.process import process as process_function
            return process_function(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during processing: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, {}


def attributes(inputparams):
    """
    Validates and returns Goldman encoding attributes based on input parameters.
    """
    encoding_method = getattr(inputparams, 'encoding_method', 'goldman')
    assembly_structure = 'synthesis'

    # Goldman specific parameters with defaults
    sequence_length = getattr(inputparams, 'sequence_length', 200)

    return {
        'encoding_method': encoding_method,
        'assembly_structure': assembly_structure,
        'sequence_length': sequence_length,
    }
