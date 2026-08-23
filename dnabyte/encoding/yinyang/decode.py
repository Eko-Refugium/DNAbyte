import tempfile
import os
import traceback

from dnabyte.encoding.yinyang.yyc import pipeline


def decode(data, params, logger):

    temp_dna = None
    temp_output = None

    try:

        # Write DNA sequences to temporary file
        temp_dna = tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            suffix=".dna"
        )


        for seq in data.data:
            temp_dna.write(seq + "\n")

        temp_dna.close()


        with open(temp_dna.name) as f:
            lines = [l.rstrip("\n") for l in f]
        

        # Temporary decoded binary output
        temp_output = tempfile.NamedTemporaryFile(
            delete=False,
            suffix=".bin"
        )
        temp_output.close()

        # Decode using the saved YYC model
        pipeline.decode(
            model_path=params.yinyang_model,
            input_path=temp_dna.name,
            output_path=temp_output.name,
            has_index=True,
            need_log=False,
        )


        # Read decoded binary
        with open(temp_output.name, "rb") as f:
            binary = f.read()


        # Convert bytes back to bitstream
        bits = "".join(
            format(byte, "08b")
            for byte in binary
        )

        # Remove byte padding that was added during encode
        bits = bits[:params.yinyang_total_bits]

        return bits, True, {}

    except Exception as e:

        if logger:
            logger.error(f"YYC decode error: {e}")
            logger.error(traceback.format_exc())

        return None, False, {}

    finally:

        for path in [
            temp_dna.name if temp_dna else None,
            temp_output.name if temp_output else None,
        ]:
            if path and os.path.exists(path):
                try:
                    os.remove(path)
                except Exception:
                    pass