import os
import pickle
import subprocess
import tempfile
import traceback


def decode(self, data, logger):

    temp_input = None
    temp_output = None

    try:

        temp_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pkl")
        temp_output = tempfile.NamedTemporaryFile(delete=False, suffix=".pkl")

        temp_input.close()
        temp_output.close()

        with open(temp_input.name, "wb") as f:
            pickle.dump(data.data, f)

        python = getattr(
            self.params,
            "dna_aeon_python",
            os.path.join(
                os.path.dirname(__file__),
                "NOREC4DNA",
                "venv",
                "Scripts",
                "python.exe",
            ),
        )

        worker = os.path.join(
            os.path.dirname(__file__),
            "decode_worker.py",
        )

        subprocess.run(
            [
                python,
                worker,
                temp_input.name,
                temp_output.name,
            ],
            check=True,
        )

        with open(temp_output.name, "rb") as f:
            bits = pickle.load(f)

        return bits, True, {}

    except Exception as e:

        if logger:
            logger.error(e)
            logger.error(traceback.format_exc())

        return None, False, {}

    finally:

        for x in [temp_input, temp_output]:
            if x and os.path.exists(x.name):
                os.remove(x.name)