import os
import sys
import subprocess
import tempfile
import pickle
import traceback

from dnabyte.encode import Encode


class DNAAeon(Encode):

    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def encode(self, data):

        temp_bits = None
        temp_result = None
        temp_info = None

        try:

            temp_bits = tempfile.NamedTemporaryFile(delete=False, suffix=".pkl")
            temp_result = tempfile.NamedTemporaryFile(delete=False, suffix=".pkl")
            temp_info = tempfile.NamedTemporaryFile(delete=False, suffix=".pkl")

            temp_bits.close()
            temp_result.close()
            temp_info.close()

            with open(temp_bits.name, "wb") as f:
                pickle.dump(data.data, f)

            python = os.path.join(
                os.path.dirname(__file__),
                "DNA_Aeon",
                "venv",
                "Scripts",
                "python.exe",
            )

            worker = os.path.join(
                os.path.dirname(__file__),
                "encode_worker.py",
            )

            print("Python:", python)
            print("Worker:", worker)
            print("Python exists:", os.path.exists(python))
            print("Worker exists:", os.path.exists(worker))

            result = subprocess.run(
                [
                    python,
                    worker,
                    temp_bits.name,
                    temp_result.name,
                    temp_info.name,
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            print(result.stdout)
            print(result.stderr)
            print(result.returncode)
            print("Output pickle:", output_pickle)
            print("Exists:", os.path.exists(output_pickle))
            print("Size:", os.path.getsize(output_pickle))
            with open(temp_result.name, "rb") as f:
                dna = pickle.load(f)

            with open(temp_info.name, "rb") as f:
                info = pickle.load(f)

            return dna, info

        except Exception as e:

            if self.logger:
                self.logger.error(e)
                self.logger.error(traceback.format_exc())

            return None, {}

        finally:

            for x in [temp_bits, temp_result, temp_info]:
                if x and os.path.exists(x.name):
                    os.remove(x.name)

    def decode(self, data):
        from dnabyte.encoding.dna_aeon.decode import decode
        return decode(self, data, self.logger)

    def process(self, data):
        from dnabyte.encoding.dna_aeon.process import process
        return process(data, self.params, self.logger)
        
def attributes(inputparams):
    """Return (and validate) DNA-Aeon-specific encoding parameters."""
    encoding_method = getattr(inputparams, 'encoding_method', 'dna_aeon')
    assembly_structure = 'synthesis'

    dna_aeon_chunk_size = int(getattr(inputparams, 'dna_aeon_chunk_size', 10))
    dna_aeon_overhead = float(getattr(inputparams, 'dna_aeon_overhead', 0.40))
    dna_aeon_insert_header = bool(getattr(inputparams, 'dna_aeon_insert_header', False))
    dna_aeon_error_correction = str(getattr(inputparams, 'dna_aeon_error_correction', 'crc'))
    # Validate error correction
    allowed_ec = {'crc', 'nocode', 'reedsolomon', 'dna_reedsolomon'}
    if dna_aeon_error_correction not in allowed_ec:
        dna_aeon_error_correction = 'crc'
    dna_aeon_repair_symbols = int(getattr(inputparams, 'dna_aeon_repair_symbols', 2))
    dna_aeon_use_dna_rules = bool(getattr(inputparams, 'dna_aeon_use_dna_rules', True))
    dna_aeon_drop_upper_bound = float(getattr(inputparams, 'dna_aeon_drop_upper_bound', 0.5))
    sequence_length = int(getattr(inputparams, 'sequence_length', 200))

    return {
        'encoding_method': encoding_method,
        'assembly_structure': assembly_structure,
        'sequence_length': sequence_length,
        'dna_aeon_chunk_size': dna_aeon_chunk_size,
        'dna_aeon_overhead': dna_aeon_overhead,
        'dna_aeon_insert_header': dna_aeon_insert_header,
        'dna_aeon_error_correction': dna_aeon_error_correction,
        'dna_aeon_repair_symbols': dna_aeon_repair_symbols,
        'dna_aeon_use_dna_rules': dna_aeon_use_dna_rules,
        'dna_aeon_drop_upper_bound': dna_aeon_drop_upper_bound,
    }
