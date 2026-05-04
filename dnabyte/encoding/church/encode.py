import tempfile
import os
import traceback

from dnabyte.encode import Encode
from dnabyte.encoding.church.StorageD.codec import ChurchEncode


class Church(Encode):
    """
    This class provides Church encoding for DNA sequences.
    Uses the ChurchEncode from dnabyte.encoding.church.StorageD.codec for encoding.
    
    Reference:
    Church, G. M., et al. (2012). "Next-generation digital information storage in DNA."
    Science 337(6102): 1628. DOI: 10.1126/science.1226355
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def encode(self, data):
        """
        Encodes binary data using Church encoding.
        Writes binary data to temporary file, encodes using ChurchEncode,
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

            # Use Church encoding
            # Church does not use rule_num, min_gc, max_gc (unlike Wukong)
            encode_worker = ChurchEncode(
                input_file_path=temp_file_path,
                output_dir=output_dir,
                sequence_length=getattr(self.params, 'sequence_length', None) or 200,
                max_homopolymer=getattr(self.params, 'max_homopolymer', None) or 6,
                rs_num=getattr(self.params, 'rs_num', None) if getattr(self.params, 'rs_num', None) is not None else 0,
                add_redundancy=getattr(self.params, 'add_redundancy', None) if getattr(self.params, 'add_redundancy', None) is not None else False,
                add_primer=getattr(self.params, 'add_primer', None) if getattr(self.params, 'add_primer', None) is not None else False,
                primer_length=getattr(self.params, 'primer_length', None) or 20
            )

            # Perform encoding
            result_file = encode_worker.common_encode()

            if self.logger:
                self.logger.info(f"Encoding completed. Result file: {result_file}")

            # Read result DNA sequences and metadata from the original FASTA
            dna_codewords = []
            left_primer = ""
            right_primer = ""
            original_fasta_metadata = ""

            if os.path.exists(result_file):
                with open(result_file, 'r') as f:
                    raw_content = f.read()

                lines = [line.strip() for line in raw_content.split('\n')]

                for line in lines:
                    if line.startswith('>') and 'totalBit:' in line:
                        # This is the original metadata header from ChurchEncode
                        original_fasta_metadata = line
                        if 'leftPrimer:' in line:
                            try:
                                left_primer = line.split('leftPrimer:')[1].split(',')[0]
                            except:
                                left_primer = ""
                        if 'rightPrimer:' in line:
                            try:
                                right_primer = line.split('rightPrimer:')[1].split(',')[0]
                            except:
                                right_primer = ""
                    elif line.startswith('>'):
                        # Skip sequence headers like >seq_1
                        continue
                    else:
                        # Only add DNA sequences (non-empty, non-header)
                        if line and not line.startswith('seq_'):
                            dna_codewords.append(line)

            # Store original FASTA metadata on params so decode can use it
            self.params.church_fasta_metadata = original_fasta_metadata

            # Calculate parameters
            barcode_length = getattr(self.params, 'primer_length', 20)
            total_bits = len(data.data)

            # Get primers from params if not found in result
            if not left_primer:
                left_primer = getattr(self.params, 'left_primer', '')
            if not right_primer:
                right_primer = getattr(self.params, 'right_primer', '')

            metadata = original_fasta_metadata

            # Create info dictionary
            info = {
                "number_of_codewords": len(dna_codewords),
                "result_file": result_file,
                "data_length": len(data.data),
                "barcode_length": barcode_length,
                "metadata": metadata,
                "left_primer": left_primer,
                "right_primer": right_primer,
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
        Decodes DNA sequences using Church decoding.
        """
        try:
            from dnabyte.encoding.church.decode import decode as decode_function
            result = decode_function(data, self.params, self.logger)
            return result
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    def process(self, data):
        """
        Processes DNA sequences using Church processing.
        """
        try:
            from dnabyte.encoding.church.process import process as process_function
            return process_function(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during processing: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, {}


def attributes(inputparams):
    """
    Validates and returns Church encoding attributes based on input parameters.
    """
    encoding_method = getattr(inputparams, 'encoding_method', 'church')
    assembly_structure = 'synthesis'

    # Church specific parameters with defaults
    sequence_length = getattr(inputparams, 'sequence_length', 200)
    max_homopolymer = getattr(inputparams, 'max_homopolymer', 6)
    rs_num = getattr(inputparams, 'rs_num', 0)
    add_redundancy = getattr(inputparams, 'add_redundancy', False)
    add_primer = getattr(inputparams, 'add_primer', True)
    primer_length = getattr(inputparams, 'primer_length', 20)

    return {
        'encoding_method': encoding_method,
        'assembly_structure': assembly_structure,
        'sequence_length': sequence_length,
        'max_homopolymer': max_homopolymer,
        'rs_num': rs_num,
        'add_redundancy': add_redundancy,
        'add_primer': add_primer,
        'primer_length': primer_length,
    }
