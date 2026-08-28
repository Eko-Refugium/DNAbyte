import tempfile
import os
import traceback
import math

from dnabyte.encode import Encode
from dnabyte.encoding.yinyang.yyc import pipeline
from dnabyte.encoding.yinyang.yyc import scheme
from dnabyte.encoding.yinyang.yyc.utils import model_saver


class YinYang(Encode):

    def __init__(self, params, logger=None):

        self.params = params
        self.logger = logger

        # YYC expects math.pow() to return an integer
        # in some places, e.g. random.randint().
        if not getattr(math, "_yinyang_pow_patched", False):

            _original_pow = math.pow

            def _int_pow(x, y):
                return int(_original_pow(x, y))

            math.pow = _int_pow
            math._yinyang_pow_patched = True

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
                1000
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

            # Use the YYC method created in __init__
            method = self.method


            # --------------------------------------------------
            # Temporary files
            # --------------------------------------------------

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


            # --------------------------------------------------
            # Convert binary string to bytes
            # --------------------------------------------------

            original_bits = data.data

            bits = original_bits

            if len(bits) % 8 != 0:
                bits += "0" * (
                    8 - len(bits) % 8
                )

            binary = bytes(
                int(bits[i:i + 8], 2)
                for i in range(
                    0,
                    len(bits),
                    8
                )
            )


            temp_input.write(binary)
            temp_input.close()

            temp_output.close()


            # --------------------------------------------------
            # YYC encode
            # --------------------------------------------------

            pipeline.encode(
                method=method,
                input_path=temp_input.name,
                output_path=temp_output.name,
                model_path=temp_model,
                need_index=True,
                need_log=False
            )


            # --------------------------------------------------
            # Read DNA sequences
            # --------------------------------------------------

            dna_sequences = []

            with open(
                temp_output.name,
                "r"
            ) as f:

                for line in f:

                    line = line.strip()

                    if (
                        line
                        and not line.startswith(">")
                    ):
                        dna_sequences.append(line)


            # --------------------------------------------------
            # Store everything needed by decoder
            # --------------------------------------------------

            self.params.yinyang_model = temp_model

            self.params.yinyang_total_bits = len(
                original_bits
            )

            self.params.yinyang_method = method


            # Make sure model was actually written
            model_saver.load_model(
                temp_model
            )


            # --------------------------------------------------
            # Information returned to framework
            # --------------------------------------------------

            info = {
                "number_of_codewords": len(
                    dna_sequences
                ),
                "data_length": len(
                    original_bits
                ),
                "total_bits": len(
                    original_bits
                ),
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


        finally:

            # The model must NOT be deleted because
            # the decoder needs it.

            if temp_input is not None:

                try:
                    os.unlink(
                        temp_input.name
                    )
                except Exception:
                    pass


            if temp_output is not None:

                try:
                    os.unlink(
                        temp_output.name
                    )
                except Exception:
                    pass


    def decode(self, data):

        try:

            from dnabyte.encoding.yinyang.decode import decode

            return decode(
                data,
                self.params,
                self.logger
            )

        except Exception as e:

            if self.logger:

                self.logger.error(
                    f"YYC decode error: {e}"
                )

                self.logger.error(
                    traceback.format_exc()
                )

            raise


    def process(self, data):

        try:

            from dnabyte.encoding.yinyang.process import process

            return process(
                data,
                self.params,
                self.logger
            )

        except Exception as e:

            if self.logger:

                self.logger.error(
                    f"YYC process error: {e}"
                )

                self.logger.error(
                    traceback.format_exc()
                )

            return None, {}


def attributes(inputparams):

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
            1000
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