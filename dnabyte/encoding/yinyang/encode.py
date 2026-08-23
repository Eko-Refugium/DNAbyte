import tempfile
import os
import traceback
import math

from dnabyte.encode import Encode
from dnabyte.encoding.yinyang.yyc import pipeline
from dnabyte.encoding.yinyang.yyc import scheme

from dnabyte.encoding.yinyang.yyc.utils import index_operator
from dnabyte.encoding.yinyang.yyc.utils import model_saver


class YinYang(Encode):

    def __init__(self, params, logger=None):

        self.params = params
        self.logger = logger

        _original_pow = math.pow
        
        def _int_pow(x, y):
            return int(_original_pow(x, y))

        math.pow = _int_pow

        self.method = scheme.YYC(
            support_bases=["A"],
            base_reference=[0, 1, 0, 1],
            current_code_matrix=[
                [1, 1, 0, 0],
                [1, 0, 0, 1],
                [1, 1, 0, 0],
                [1, 1, 0, 0],
            ],
            search_count=getattr(
                params,
                "yinyang_search_count",
                100
            ),
            max_homopolymer=getattr(
                params,
                "max_homopolymer",
                4
            ),
            max_content=getattr(
                params,
                "max_content",
                0.6
            )
        )


    def encode(self, data):

        temp_input = None
        temp_output = None
        temp_model = None

        try:

            # temporary files
            temp_input = tempfile.NamedTemporaryFile(
                delete=False,
                suffix=".bin"
            )

            temp_output = tempfile.NamedTemporaryFile(
                delete=False,
                suffix=".dna"
            )

            temp_model_dir = tempfile.mkdtemp()

            temp_model = os.path.join(
                temp_model_dir,
                "yyc_model.pkl"
            )
            


            # convert bit string to bytes
            bits = data.data
            bitstream = bits

            if len(bits) % 8 != 0:
                bits += "0" * (8 - len(bits) % 8)

            binary = bytes(
                int(bits[i:i+8], 2)
                for i in range(0, len(bits), 8)
            )


            temp_input.write(binary)
            temp_input.close()

            temp_output.close()


            pipeline.encode(
                method=self.method,
                input_path=temp_input.name,
                output_path=temp_output.name,
                model_path=temp_model,
                need_index=True,
                need_log=False
            )


            dna_sequences = []

            with open(
                temp_output.name,
                "r"
            ) as f:

                for line in f:

                    line=line.strip()

                    if line and not line.startswith(">"):
                        dna_sequences.append(line)

            # keep model for decoding
            self.params.yinyang_model = temp_model
            self.params.yinyang_total_bits = len(data.data)
            self.params.yingyang_method = self.method

            model = model_saver.load_model(temp_model)


            info = {
                "number_of_codewords": len(dna_sequences),
                "data_length": len(bitstream),
                "total_bits": len(bitstream),
                "barcode_length": 0,
                "metadata": "",
            }


            return dna_sequences, info


        except Exception as e:

            if self.logger:
                self.logger.error(
                    f"YYC encode error: {e}"
                )
                self.logger.error(
                    traceback.format_exc()
                )

            return None, {}



        
    def decode(self, data):
        try:
            from dnabyte.encoding.yinyang.decode import decode
            return decode(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(
                    f"YYC decode error: {e}"
                )
                self.logger.error(
                    traceback.format_exc()
                )

            raise e
        
    def process(self, data):
        try:
            from dnabyte.encoding.yinyang.process import process
            return process(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"YYC process error: {e}")
                self.logger.error(traceback.format_exc())
            return None, {}
        
def attributes(inputparams):
    """
    Return Yin-Yang encoding parameters.
    Called by Params.__init__ through the plugin system.
    """

    encoding_method = getattr(
        inputparams,
        "encoding_method",
        "yinyang"
    )

    assembly_structure = "synthesis"

    sequence_length = int(
        getattr(
            inputparams,
            "sequence_length",
            120
        )
    )

    yinyang_search_count = int(
        getattr(
            inputparams,
            "yinyang_search_count",
            100
        )
    )

    max_homopolymer = int(
        getattr(
            inputparams,
            "max_homopolymer",
            4
        )
    )

    max_content = float(
        getattr(
            inputparams,
            "max_content",
            0.6
        )
    )

    return {
        "encoding_method": encoding_method,
        "assembly_structure": assembly_structure,
        "sequence_length": sequence_length,
        "yinyang_search_count": yinyang_search_count,
        "max_homopolymer": max_homopolymer,
        "max_content": max_content,
    }