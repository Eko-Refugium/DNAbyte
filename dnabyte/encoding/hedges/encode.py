from dnabyte.encode import Encode
import traceback
from dnabyte.encoding.hedges.hedges_python import encode_hedges
import dnabyte.encoding.hedges.hedges_python as hedges_python

print(hedges_python.__file__)
print(hedges_python.encode_hedges)

def bitstring_to_bytes(s):
    v = int(s, 2)
    b = bytearray()
    while v:
        b.append(v & 0xff)
        v >>= 8
    return bytes(b[::-1])

class HEDGES(Encode):

    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger


    def encode(self, data):

        """
        data.data = binary string
        """

        print(f"Encoding data with HEDGES: {data.data[:50]}... (length: {len(data.data)})")

        byte_data = int(data.data, 2).to_bytes(
            (len(data.data)+7)//8,
            byteorder="big"
        )

        print(hedges_python.__file__)
        print(hedges_python.encode_hedges.__doc__)
        print(hedges_python.encode_hedges)

        dna_strands = encode_hedges(byte_data)

        info = {
            "number_of_strands": len(dna_strands),
            "strand_length": len(dna_strands[0]),
            "barcode_length": 0,  # Placeholder, update with actual barcode length if applicable
        }

        return dna_strands, info


    def decode(self, data):
        """
        Decodes DNA sequences using HEDGES decoding.
        """
        try:
            from dnabyte.encoding.hedges.decode import decode as decode_function
            result = decode_function(data, self.params, self.logger)
            return result
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    def process(self, data):
        """
        Processes DNA sequences using HEDGES processing.
        """
        try:
            from dnabyte.encoding.hedges.process import process as process_function
            return process_function(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during processing: {str(e)}")
                self.logger.error(traceback.format_exc())
            return None, {}



def attributes(inputparams):
    """
    Validates and returns HEDGES encoding attributes based on input parameters.
    """
    encoding_method = getattr(inputparams, 'encoding_method', 'hedges')
    assembly_structure = 'synthesis'

    # HEDGES specific parameters with defaults
    sequence_length = getattr(inputparams, 'sequence_length', 200)

    # Build return dictionary
    attributes_dict = {
        "encoding_method": encoding_method,
        "assembly_structure": assembly_structure,
        "sequence_length": sequence_length,
    }

    return attributes_dict
