import os
import sys
import tempfile
import traceback

from dnabyte.encode import Encode

# Ensure the NOREC4DNA package is importable
_NOREC4DNA_DIR = os.path.join(os.path.dirname(__file__), 'DNA_Aeon', 'NOREC4DNA')
if _NOREC4DNA_DIR not in sys.path:
    sys.path.insert(0, _NOREC4DNA_DIR)

from norec4dna.Encoder import Encoder as NEncoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode, get_error_correction_decode
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from norec4dna.helper.RU10Helper import intermediate_symbols

# ---------------------------------------------------------------------------
# Format-string constants — must match between encoder and decoder
# ---------------------------------------------------------------------------
ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"


# ---------------------------------------------------------------------------
# Plugin class
# ---------------------------------------------------------------------------

class DNAAeon(Encode):
    """
    DNA-Aeon encoding for DNA storage.

    Uses the NOREC4DNA RU10 Raptor fountain code to encode binary data
    into DNA oligos.  The fountain code provides resilience against
    oligo losses (erasures) through redundancy.

    Parameters (set via ``Params``):
    - ``dna_aeon_chunk_size`` – bytes per data chunk (default 10)
    - ``dna_aeon_overhead`` – fractional overhead for fountain redundancy (default 0.40)
    - ``dna_aeon_insert_header`` – whether to store a header chunk (default False)
    - ``dna_aeon_error_correction`` – error correction method name (default 'crc')
    - ``dna_aeon_repair_symbols`` – number of RS repair symbols (default 2)
    - ``dna_aeon_use_dna_rules`` – enforce ACGT homopolymer/GC rules (default True)
    - ``dna_aeon_drop_upper_bound`` – upper bound for dropping rule-violating packets (default 0.5)

    Reference
    ---------
    DNA-Aeon: NOREC4DNA RU10 implementation of Raptor fountain codes for
    DNA data storage.
    """

    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    # ---- encode ----------------------------------------------------------

    def encode(self, data):
        """
        Encode a bitstream into DNA oligos using the NOREC4DNA RU10 fountain code.

        1. Convert the bitstream to bytes and write to a temporary file.
        2. Create an ``RU10Encoder`` with the configured parameters.
        3. Encode to packets.
        4. Extract DNA strings from each encoded packet.
        5. Return the list of DNA strings and metadata.
        """
        try:
            chunk_size = int(getattr(self.params, 'dna_aeon_chunk_size', 10))
            overhead = float(getattr(self.params, 'dna_aeon_overhead', 0.40))
            insert_header = bool(getattr(self.params, 'dna_aeon_insert_header', False))
            error_correction_name = str(getattr(self.params, 'dna_aeon_error_correction', 'crc'))
            repair_symbols = int(getattr(self.params, 'dna_aeon_repair_symbols', 2))
            use_dna_rules = bool(getattr(self.params, 'dna_aeon_use_dna_rules', True))
            drop_upper_bound = float(getattr(self.params, 'dna_aeon_drop_upper_bound', 0.5))

            # Get the error correction encode function
            error_correction = get_error_correction_encode(error_correction_name, repair_symbols)

            bitstream = data.data  # string of '0'/'1'
            total_bits = len(bitstream)

            # Convert bitstring to bytes
            # Pad to a multiple of 8 bits
            padded = bitstream
            if len(padded) % 8 != 0:
                padded = padded + '0' * (8 - len(padded) % 8)
            raw_bytes = int(padded, 2).to_bytes(len(padded) // 8, byteorder='big')

            # Write bytes to a temporary file (NOREC4DNA requires a file path)
            tmp_fd, tmp_path = tempfile.mkstemp(suffix='.bin')
            try:
                with os.fdopen(tmp_fd, 'wb') as tmp_f:
                    tmp_f.write(raw_bytes)

                # Calculate number of chunks from file size and chunk_size
                file_size = len(raw_bytes)

                # Ensure we have enough chunks for the Raptor code
                # (minimum ~4 chunks needed for intermediate block generation)
                effective_chunk_size = chunk_size
                num_chunks_estimate = NEncoder.get_number_of_chunks_for_file_with_chunk_size(
                    tmp_path, effective_chunk_size, insert_header=insert_header
                )

                # If too few chunks, reduce chunk_size to get at least 4
                min_chunks = 4
                while num_chunks_estimate < min_chunks and effective_chunk_size > 1:
                    effective_chunk_size = max(1, effective_chunk_size // 2)
                    num_chunks_estimate = NEncoder.get_number_of_chunks_for_file_with_chunk_size(
                        tmp_path, effective_chunk_size, insert_header=insert_header
                    )

                dist = RaptorDistribution(num_chunks_estimate)
                rules = FastDNARules() if use_dna_rules else None

                # ----------------------------------------------------------
                # Compute minimum overhead from the Raptor intermediate block
                # The GEPP solver needs at least L = K + S + H rows (packets)
                # to be solvable, where S and H are LDPC/Half check symbols.
                # We use 2× L to provide ample headroom: the random packet
                # selection does not guarantee linearly independent equations,
                # and in a full error pipeline some packets may be lost.
                # ----------------------------------------------------------
                L, S, H = intermediate_symbols(num_chunks_estimate, dist)
                min_packets = L * 2
                min_overhead = (min_packets / num_chunks_estimate) - 1.0
                effective_overhead = max(overhead, min_overhead)

                if effective_overhead > overhead and self.logger:
                    self.logger.info(
                        f"DNA-Aeon: overhead {overhead:.2f} too low for "
                        f"{num_chunks_estimate} chunks (L={L}, S={S}, H={H}). "
                        f"Auto-increased to {effective_overhead:.2f}"
                    )

                # Create the RU10 encoder
                # Pass chunk_size so the encoder recalculates number_of_chunks internally
                encoder = RU10Encoder(
                    tmp_path,
                    num_chunks_estimate,
                    dist,
                    insert_header=insert_header,
                    pseudo_decoder=None,
                    chunk_size=effective_chunk_size,
                    rules=rules,
                    error_correction=error_correction,
                    packet_len_format=PACKET_LEN_FORMAT,
                    crc_len_format=CRC_LEN_FORMAT,
                    number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                    id_len_format=ID_LEN_FORMAT,
                    save_number_of_chunks_in_packet=True,
                    drop_upper_bound=drop_upper_bound,
                )
                encoder.set_overhead_limit(effective_overhead)

                # Encode
                encoder.encode_to_packets()

                # Extract DNA strings from each packet
                dna_oligos = []
                for packet in encoder.encodedPackets:
                    dna_str = packet.get_dna_struct(True)
                    dna_oligos.append(dna_str)

                # Store metadata on params for decode
                actual_number_of_chunks = encoder.number_of_chunks
                self.params.dna_aeon_total_bits = total_bits
                self.params.dna_aeon_number_of_chunks = actual_number_of_chunks
                self.params.dna_aeon_number_of_packets = len(dna_oligos)

                if self.logger:
                    self.logger.info(
                        f"DNA-Aeon encoded {len(dna_oligos)} oligos, "
                        f"chunk_size={effective_chunk_size}, overhead={effective_overhead:.2f}, "
                        f"number_of_chunks={actual_number_of_chunks}, "
                        f"insert_header={insert_header}"
                    )

                info = {
                    'number_of_codewords': len(dna_oligos),
                    'data_length': total_bits,
                    'barcode_length': 0,
                    'metadata': '',
                    'left_primer': '',
                    'right_primer': '',
                    'total_bits': total_bits,
                    'original_fasta_metadata': '',
                    'dna_aeon_chunk_size': effective_chunk_size,
                    'dna_aeon_overhead': effective_overhead,
                    'dna_aeon_number_of_chunks': actual_number_of_chunks,
                    'dna_aeon_number_of_packets': len(dna_oligos),
                }

                return dna_oligos, info

            finally:
                # Clean up temporary file
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

        except Exception as e:
            if self.logger:
                self.logger.error(f"DNA-Aeon encode error: {e}")
                self.logger.error(traceback.format_exc())
            return None, {}

    # ---- decode ----------------------------------------------------------

    def decode(self, data):
        """Decode DNA-Aeon oligos back to a bitstream."""
        try:
            from dnabyte.encoding.dna_aeon.decode import decode as decode_fn
            return decode_fn(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"DNA-Aeon decode error: {e}")
                self.logger.error(traceback.format_exc())
            return None, False, {}

    # ---- process ---------------------------------------------------------

    def process(self, data):
        """Process (consensus) sequenced DNA-Aeon oligos."""
        try:
            from dnabyte.encoding.dna_aeon.process import process as process_fn
            return process_fn(data, self.params, self.logger)
        except Exception as e:
            if self.logger:
                self.logger.error(f"DNA-Aeon process error: {e}")
                self.logger.error(traceback.format_exc())
            return None, {}


# ---------------------------------------------------------------------------
# attributes() — called by Params.__init__ via the plugin system
# ---------------------------------------------------------------------------

def attributes(inputparams):
    """Return (and validate) DNA-Aeon-specific encoding parameters."""
    encoding_method = getattr(inputparams, 'encoding_method', 'dna_aeon')
    assembly_structure = 'synthesis'

    dna_aeon_chunk_size = int(getattr(inputparams, 'dna_aeon_chunk_size', 10))
    dna_aeon_overhead = float(getattr(inputparams, 'dna_aeon_overhead', 0.40))
    dna_aeon_insert_header = bool(getattr(inputparams, 'dna_aeon_insert_header', False))
    dna_aeon_error_correction = str(getattr(inputparams, 'dna_aeon_error_correction', 'crc'))
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
