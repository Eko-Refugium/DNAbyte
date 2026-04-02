import os
import sys
import math
import traceback

import numpy as np

# Ensure the GC+ DNA source directory is importable
_GCP_DNA_DIR = os.path.join(os.path.dirname(__file__), 'src', 'GCPdna')
if _GCP_DNA_DIR not in sys.path:
    sys.path.insert(0, _GCP_DNA_DIR)

from dnabyte.encoding.gcplus.src.GCPdna.GCP_Decode_DNA import GCP_Decode_DNA_brute
from dnabyte.encoding.gcplus.src.GCPdna.preCompute_Patterns import preCompute_Patterns


def decode(data, params, logger=None):
    """
    Decode GC+ DNA codewords back into a bitstream.

    Each DNA codeword is independently decoded with ``GCP_Decode_DNA_brute``.
    Decoded *k*-bit chunks are concatenated and truncated to the original
    total-bits length.

    Args:
        data:   Data object whose ``.data`` is a list of DNA codeword strings.
        params: Parameters object carrying GC+ metadata from encoding.
        logger: Optional logger.

    Returns:
        (binary_data_str, valid, info)
    """
    binary_data = ""
    valid = False
    info = {}

    try:
        k = int(getattr(params, 'gcplus_k', 168))
        l = int(getattr(params, 'gcplus_l', 8))
        c1 = int(getattr(params, 'gcplus_c1', 2))
        total_bits = int(getattr(params, 'gcplus_total_bits', 0))

        # Derived constants
        K = int(math.ceil(k / l))
        c2 = 1
        N = K + c1 + c2
        q = 2 ** l
        len_last = (k - 1) % l + 1
        n_expected = int(getattr(params, 'gcplus_n', 0))

        # Pre-compute decoder patterns
        lim = 5
        lambda_depths = [0] * lim
        lambda_depths[0] = 1
        lambda_depths[1] = 1
        patterns = preCompute_Patterns(lambda_depths, K, len_last, lim, c1)

        # Load codebook
        codebook_path = os.path.join(_GCP_DNA_DIR, 'codebook_DNA.txt')
        with open(codebook_path, 'r', encoding='utf-8') as f:
            codebook = [line.strip() for line in f if line.strip()]
        d_min = 5

        dna_sequences = data.data if isinstance(data.data, list) else [data.data]
        clean_sequences = [
            str(seq).replace(' ', '').strip()
            for seq in dna_sequences
            if str(seq).replace(' ', '').strip()
        ]

        if not clean_sequences:
            if logger:
                logger.warning("GC+ decode: no sequences to decode")
            return "", False, {}

        if logger:
            logger.info(f"GC+ decoding {len(clean_sequences)} codewords")

        # Determine n from first codeword if not stored
        if n_expected == 0:
            n_expected = len(clean_sequences[0])

        decoded_bits = []
        success_count = 0
        fail_count = 0

        for idx, cw in enumerate(clean_sequences):
            try:
                uhat, _ = GCP_Decode_DNA_brute(
                    cw, n_expected, k, l, N, K, c1, q,
                    len_last, lim, patterns, codebook, d_min
                )
                if uhat and len(uhat) > 0:
                    bits_str = ''.join(str(b) for b in uhat)
                    # Ensure exactly k bits
                    if len(bits_str) < k:
                        bits_str += '0' * (k - len(bits_str))
                    elif len(bits_str) > k:
                        bits_str = bits_str[:k]
                    decoded_bits.append(bits_str)
                    success_count += 1
                else:
                    # Decode failure — fill with zeros
                    decoded_bits.append('0' * k)
                    fail_count += 1
            except Exception:
                decoded_bits.append('0' * k)
                fail_count += 1

        binary_data = ''.join(decoded_bits)

        # Truncate to original length
        if total_bits > 0:
            binary_data = binary_data[:total_bits]

        valid = fail_count == 0 and len(binary_data) > 0

        if logger:
            logger.info(
                f"GC+ decoded: {success_count} ok, {fail_count} failed, "
                f"{len(binary_data)} bits"
            )

        info = {
            'decoded_codewords': success_count,
            'failed_codewords': fail_count,
            'data_length': len(binary_data),
            'valid': valid,
            'total_bits': total_bits,
        }

    except Exception as e:
        if logger:
            logger.error(f"GC+ decode error: {e}")
            logger.error(traceback.format_exc())
        binary_data = ""
        valid = False
        info = {'error': str(e)}

    return binary_data, valid, info
