from dnabyte.encode import Encode
import traceback
from dnabyte.encoding.hedges.hedges_python import encode_hedges
import dnabyte.encoding.hedges.hedges_python as hedges_python


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
        self.hedges_cfg = self.create_config(params)
        self.params.hedges_cfg = self.hedges_cfg


    def create_config(self, params):

        cfg = hedges_python.HedgesConfig()

        # map framework parameters -> HEDGES config

        cfg.coderate = getattr(
            params,
            "hedges_coderate",
            3
        )

        cfg.total_strand_length = getattr(
            params,
            "sequence_length",
            300
        )

        cfg.strand_id_bytes = getattr(
            params,
            "strand_id_bytes",
            2
        )

        cfg.strand_runout_bytes = getattr(
            params,
            "strand_runout_bytes",
            2
        )

        cfg.gc_window = getattr(
            params,
            "gc_window",
            12
        )

        cfg.gc_max = getattr(
            params,
            "gc_max",
            8
        )

        cfg.max_homopolymer = getattr(
            params,
            "max_homopolymer",
            4
        )

        cfg.heap_limit = getattr(
            params,
            "heap_limit",
            1000000
        )

        cfg.left_primer = getattr(
            params,
            "left_primer",
            "TCGAAGTCAGCGTGTATTGTATG"
        )

        cfg.right_primer = getattr(
            params,
            "right_primer",
            "TAGTGAGTGCGATTAAGCGTGTT"
        )


        return cfg



    def encode(self, data):

        """
        data.data = binary string
        """


        byte_data = int(data.data, 2).to_bytes(
            (len(data.data)+7)//8,
            byteorder="big"
        )


        dna_strands = encode_hedges(byte_data, self.hedges_cfg)

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
    attributes_dict = {

        "encoding_method":
            getattr(
                inputparams,
                "encoding_method",
                "hedges"
            ),

        "assembly_structure":
            "synthesis",

        "sequence_length":
            getattr(
                inputparams,
                "sequence_length",
                300
            ),

        "hedges_coderate":
            getattr(
                inputparams,
                "hedges_coderate",
                3
            ),

        "gc_window":
            getattr(
                inputparams,
                "gc_window",
                12
            ),

        "gc_max":
            getattr(
                inputparams,
                "gc_max",
                8
            ),

        "max_homopolymer":
            getattr(
                inputparams,
                "max_homopolymer",
                4
            ),
    }

    return attributes_dict
