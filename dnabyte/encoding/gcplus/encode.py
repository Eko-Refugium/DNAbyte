import os
import sys
import math
import traceback

import numpy as np

from dnabyte.encode import Encode

# Ensure the GC+ DNA source directory is importable
_GCP_DNA_DIR = os.path.join(os.path.dirname(__file__), 'src', 'GCPdna')
if _GCP_DNA_DIR not in sys.path:
    sys.path.insert(0, _GCP_DNA_DIR)

from dnabyte.encoding.gcplus.src.GCPdna.GCP_Encode_DNA import GCP_Encode_DNA_brute, binary_to_dna
from dnabyte.encoding.gcplus.src.GCPdna.GCP_Decode_DNA import GCP_Decode_DNA_brute
from dnabyte.encoding.gcplus.src.GCPdna.preCompute_Patterns import preCompute_Patterns


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _load_codebook():
    """Load the DNA codebook from the GCP dna directory."""
    codebook_path = os.path.join(_GCP_DNA_DIR, 'codebook_DNA.txt')
    with open(codebook_path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f if line.strip()]


def _precompute(k, l, c1):
    """Pre-compute decoder pattern tables."""
    K = int(math.ceil(k / l))
    len_last = (k - 1) % l + 1
    lim = 5
    lambda_depths = [0] * lim
    lambda_depths[0] = 1
    lambda_depths[1] = 1
    # lambda_depths[2] = 0  (already 0)
    patterns = preCompute_Patterns(lambda_depths, K, len_last, lim, c1)
    return patterns, lim, len_last, K


# ---------------------------------------------------------------------------
# Plugin class
# ---------------------------------------------------------------------------

class GCPlus(Encode):
    """
    GC+ (Guess & Check Plus) encoding for DNA storage.

    Each oligo encodes *k* information bits using an inner systematic RS code
    plus a codebook-based parity tail, enabling correction of random edit
    errors (insertions, deletions, substitutions) at the per-oligo level.

    Reference
    ---------
    Serge Kas Hanna, "GC+ Code: A Systematic Short Blocklength Code for
    Correcting Random Edit Errors in DNA Storage", arXiv:2402.01244, 2025.
    """

    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    # ---- encode ----------------------------------------------------------

    def encode(self, data):
        """
        Encode a bitstream into GC+ DNA codewords.

        1. Split the bitstream into *k*-bit chunks.
        2. Encode each chunk with ``GCP_Encode_DNA_brute``.
        3. Return the list of DNA codeword strings.
        """
        try:
            k = int(getattr(self.params, 'gcplus_k', 168))
            l = int(getattr(self.params, 'gcplus_l', 8))
            c1 = int(getattr(self.params, 'gcplus_c1', 2))

            codebook = _load_codebook()

            bitstream = data.data  # string of '0'/'1'

            # Pad bitstream to a multiple of k
            total_bits = len(bitstream)
            if len(bitstream) % k != 0:
                bitstream = bitstream + '0' * (k - len(bitstream) % k)

            # Split into k-bit chunks (as lists of int)
            chunks = []
            for i in range(0, len(bitstream), k):
                chunk = [int(b) for b in bitstream[i:i + k]]
                chunks.append(chunk)

            dna_codewords = []
            n_val = None
            N_val = None
            K_val = None
            q_val = None

            for chunk in chunks:
                x, n, N, K, q, U, X, check_par = GCP_Encode_DNA_brute(
                    chunk, l, c1, codebook
                )
                dna_codewords.append(x)
                n_val = n
                N_val = N
                K_val = K
                q_val = q

            # Store metadata on params for decode
            self.params.gcplus_total_bits = total_bits
            self.params.gcplus_n = n_val
            self.params.gcplus_N = N_val
            self.params.gcplus_K = K_val
            self.params.gcplus_q = q_val

            if self.logger:
                self.logger.info(
                    f"GC+ encoded {len(chunks)} oligos, k={k}, l={l}, c1={c1}, "
                    f"n={n_val} bases per oligo"
                )

            info = {
                'number_of_codewords': len(dna_codewords),
                'data_length': total_bits,
                'barcode_length': 0,
                'metadata': '',
                'left_primer': '',
                'right_primer': '',
                'total_bits': total_bits,
                'original_fasta_metadata': '',
                'gcplus_k': k,
                'gcplus_l': l,
                'gcplus_c1': c1,
                'gcplus_n': n_val,
            }
            print(dna_codewords)

            return dna_codewords, info

        except Exception as e:
            if self.logger:
                self.logger.error(f"GC+ encode error: {e}")
                self.logger.error(traceback.format_exc())
            return None, {}

    # ---- decode ----------------------------------------------------------

    def decode(self, data):
        """Decode GC+ DNA codewords back to a bitstream."""
        try:
            from dnabyte.encoding.gcplus.decode import decode as decode_fn
            return decode_fn(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"GC+ decode error: {e}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    # ---- process ---------------------------------------------------------

    def process(self, data):
        """Process (consensus) sequenced GC+ DNA codewords."""
        try:
            from dnabyte.encoding.gcplus.process import process as process_fn
            return process_fn(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"GC+ process error: {e}")
                self.logger.error(traceback.format_exc())
            return None, {}


# ---------------------------------------------------------------------------
# attributes() — called by Params.__init__ via the plugin system
# ---------------------------------------------------------------------------

def attributes(inputparams):
    """Return (and validate) GC+-specific encoding parameters."""
    encoding_method = getattr(inputparams, 'encoding_method', 'gcplus')
    assembly_structure = 'synthesis'

    gcplus_k = int(getattr(inputparams, 'gcplus_k', 168))
    gcplus_l = int(getattr(inputparams, 'gcplus_l', 8))
    gcplus_c1 = int(getattr(inputparams, 'gcplus_c1', 2))
    sequence_length = int(getattr(inputparams, 'sequence_length', 200))

    return {
        'encoding_method': encoding_method,
        'assembly_structure': assembly_structure,
        'sequence_length': sequence_length,
        'gcplus_k': gcplus_k,
        'gcplus_l': gcplus_l,
        'gcplus_c1': gcplus_c1,
    }
