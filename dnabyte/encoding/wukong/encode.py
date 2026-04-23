import tempfile
import os
import traceback

from dnabyte.encode import Encode
from dnabyte.encoding.wukong.StorageD.codec import WukongEncode


class Wukong(Encode):
    """
    This class provides Wukong encoding for DNA sequences.
    Uses the WukongEncode from dnabyte.encoding.wukong.StorageD.codec for encoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def encode(self, data):
        """
        Encodes binary data using Wukong encoding.
        Writes binary data to temporary file, encodes using WukongEncode, 
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

            # Use Wukong encoding
            add_primer = getattr(self.params, 'add_primer', True)
            primer_length = getattr(self.params, 'primer_length', 20) if add_primer else 0
            encode_worker = WukongEncode(
                input_file_path=temp_file_path,
                output_dir=output_dir,
                sequence_length=getattr(self.params, 'sequence_length', 200),
                max_homopolymer=getattr(self.params, 'max_homopolymer', 4),
                min_gc=getattr(self.params, 'min_gc', 0.4),
                max_gc=getattr(self.params, 'max_gc', 0.6),
                rule_num=getattr(self.params, 'rule_num', 1),
                rs_num=getattr(self.params, 'rs_num', 0),
                add_redundancy=getattr(self.params, 'add_redundancy', True),
                add_primer=add_primer,
                primer_length=primer_length
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
                        # This is the original metadata header from WukongEncode
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
            self.params.wukong_fasta_metadata = original_fasta_metadata

            # Calculate parameters
            barcode_length = getattr(self.params, 'primer_length', 20)
            rs_num = getattr(self.params, 'rs_num', 0)
            add_redundancy = 1 if getattr(self.params, 'add_redundancy', True) else 0
            total_bits = len(data.data)

            # Get primers from params if not found in result
            if not left_primer:
                left_primer = getattr(self.params, 'left_primer', 'TTAAGGTTGCCGGGTTCACC')
            if not right_primer:
                right_primer = getattr(self.params, 'right_primer', 'CGCGAGGAGTCAGAGGTTAT')

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
        Decodes DNA sequences using Wukong decoding.
        """
        try:
            from dnabyte.encoding.wukong.decode import decode as decode_function
            result = decode_function(data, self.params, self.logger)
            return result
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    def process(self, data):
        """
        Processes DNA sequences using Wukong processing.
        """
        try:
            from dnabyte.encoding.wukong.process import process as process_function
            return process_function(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during processing: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, {}



def attributes(inputparams):
    """
    Validates and returns Wukong encoding attributes based on input parameters.
    """
    encoding_method = getattr(inputparams, 'encoding_method', 'wukong')
    assembly_structure = 'synthesis'

    # Wukong specific parameters with defaults
    sequence_length = getattr(inputparams, 'sequence_length', 200)
    max_homopolymer = getattr(inputparams, 'max_homopolymer', 4)
    min_gc = getattr(inputparams, 'min_gc', 0.4)
    max_gc = getattr(inputparams, 'max_gc', 0.6)
    rule_num = getattr(inputparams, 'rule_num', 1)
    rs_num = getattr(inputparams, 'rs_num', 0)
    add_redundancy = getattr(inputparams, 'add_redundancy', True)
    add_primer = getattr(inputparams, 'add_primer', True)
    primer_length = getattr(inputparams, 'primer_length', 20)

    # Validate parameters
    if not (0.0 <= min_gc <= 1.0 and 0.0 <= max_gc <= 1.0):
        raise ValueError("GC content values must be between 0.0 and 1.0")
    
    if min_gc > max_gc:
        raise ValueError("min_gc must be less than or equal to max_gc")
    
    if sequence_length < 1:
        raise ValueError("sequence_length must be at least 1")
    
    if max_homopolymer < 1:
        raise ValueError("max_homopolymer must be at least 1")
    
    if add_primer and primer_length < 1:
        raise ValueError("primer_length must be at least 1")

    # Build return dictionary
    attributes_dict = {
        "encoding_method": encoding_method,
        "assembly_structure": assembly_structure,
        "sequence_length": sequence_length,
        "max_homopolymer": max_homopolymer,
        "min_gc": min_gc,
        "max_gc": max_gc,
        "rule_num": rule_num,
        "rs_num": rs_num,
        "add_redundancy": add_redundancy,
        "add_primer": add_primer,
        "primer_length": primer_length,
    }

    return attributes_dict
