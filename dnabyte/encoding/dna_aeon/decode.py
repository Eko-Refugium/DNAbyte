import os
import sys
import tempfile
import traceback

from io import BytesIO

# Ensure the NOREC4DNA package is importable
_NOREC4DNA_DIR = os.path.join(os.path.dirname(__file__), 'DNA_Aeon', 'NOREC4DNA')
if _NOREC4DNA_DIR not in sys.path:
    sys.path.insert(0, _NOREC4DNA_DIR)

from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode, get_error_correction_decode
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte

# Format-string constants — must match the encoder
ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"


def decode(data, params, logger=None):
    """
    Decode DNA-Aeon oligos back into a bitstream.

    Uses the NOREC4DNA RU10Decoder to reconstruct the original data from
    DNA oligo strings produced by the encoder.

    The decode flow:
    1. Convert each DNA string back to bytes via ``tranlate_quat_to_byte``.
    2. Feed each byte packet into the RU10Decoder via ``parse_raw_packet`` /
       ``input_new_packet``.
    3. Once enough packets are received, the Gaussian elimination solver
       reconstructs the original data.
    4. Save the decoded file and read it back as bytes.
    5. Convert the bytes back to a bitstring truncated to the original length.

    Args:
        data:   Data object whose ``.data`` is a list of DNA oligo strings.
        params: Parameters object carrying DNA-Aeon metadata from encoding.
        logger: Optional logger.

    Returns:
        (binary_data_str, valid, info)
    """
    binary_data = ""
    valid = False
    info = {}

    try:
        total_bits = int(getattr(params, 'dna_aeon_total_bits', 0))
        error_correction_name = str(getattr(params, 'dna_aeon_error_correction', 'crc'))
        repair_symbols = int(getattr(params, 'dna_aeon_repair_symbols', 2))
        insert_header = bool(getattr(params, 'dna_aeon_insert_header', False))
        number_of_chunks = getattr(params, 'dna_aeon_number_of_chunks', None)

        # Get the error correction decode function
        error_correction = get_error_correction_decode(error_correction_name, repair_symbols)

        dna_sequences = data.data if isinstance(data.data, list) else [data.data]
        clean_sequences = [
            str(seq).replace(' ', '').strip()
            for seq in dna_sequences
            if str(seq).replace(' ', '').strip()
        ]

        if not clean_sequences:
            if logger:
                logger.warning("DNA-Aeon decode: no sequences to decode")
            return "", False, {}

        if logger:
            logger.info(f"DNA-Aeon decoding {len(clean_sequences)} oligos")

        # Create the RU10Decoder
        # We pass static_number_of_chunks=None so the decoder reads it from packets
        decoder = RU10Decoder(
            file=None,
            error_correction=error_correction,
            use_headerchunk=insert_header,
            static_number_of_chunks=None,
        )
        decoder.read_all_before_decode = True

        # Feed DNA oligos into the decoder one by one
        decoded = False
        corrupt_count = 0
        success_count = 0

        for dna_str in clean_sequences:
            try:
                # Convert DNA string to bytes
                raw_bytes = tranlate_quat_to_byte(dna_str)

                # Parse the raw packet
                packet = decoder.parse_raw_packet(
                    raw_bytes,
                    crc_len_format=CRC_LEN_FORMAT,
                    number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                    packet_len_format=PACKET_LEN_FORMAT,
                    id_len_format=ID_LEN_FORMAT,
                )

                if packet == "CORRUPT":
                    corrupt_count += 1
                    continue

                success_count += 1
                decoder.input_new_packet(packet)
            except Exception:
                corrupt_count += 1
                continue

        # Now attempt to solve with ALL collected packets
        if decoder.GEPP is not None and decoder.GEPP.isPotentionallySolvable():
            decoded = decoder.GEPP.solve(partial=False)

        if logger:
            logger.info(
                f"DNA-Aeon decoder: {success_count} valid, {corrupt_count} corrupt, "
                f"decoded={decoded}"
            )

        if not decoded:
            if logger:
                logger.warning("DNA-Aeon decode: insufficient packets to reconstruct data")
            return "", False, {'error': 'insufficient packets', 'valid_packets': success_count, 'corrupt_packets': corrupt_count}

        # Save decoded file to a temp location and read back the bytes
        tmp_fd, tmp_path = tempfile.mkstemp(suffix='.dec')
        os.close(tmp_fd)
        try:
            # Set the file attribute so saveDecodedFile can generate a filename
            decoder.file = tmp_path
            output_bytes = decoder.saveDecodedFile(
                last_chunk_len_format="I",
                null_is_terminator=False,
                print_to_output=False,
                partial_decoding=True,
            )

            # output_bytes should be the concatenated decoded data
            if isinstance(output_bytes, str):
                # If saveDecodedFile returned a filename, read the file
                if os.path.exists(output_bytes):
                    with open(output_bytes, 'rb') as f:
                        output_bytes = f.read()
                    os.remove(output_bytes)
                else:
                    output_bytes = output_bytes.encode('utf-8')

            if not isinstance(output_bytes, bytes):
                output_bytes = bytes(output_bytes)

            # Convert bytes back to bitstring
            bits = bin(int.from_bytes(output_bytes, byteorder='big'))[2:]
            # Pad to full byte alignment
            bits = bits.zfill(len(output_bytes) * 8)

            # Truncate to original length
            if total_bits > 0:
                binary_data = bits[:total_bits]
            else:
                binary_data = bits

            valid = len(binary_data) > 0

            if logger:
                logger.info(f"DNA-Aeon decoded {len(binary_data)} bits")

        finally:
            # Clean up any temp files generated by the decoder
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            # The decoder saves with a "DEC_" prefix 
            dec_path = os.path.join(os.path.dirname(tmp_path), "DEC_" + os.path.basename(tmp_path))
            if os.path.exists(dec_path):
                os.remove(dec_path)
            # Also check current directory for DEC_ files
            cwd_dec_path = os.path.join(os.getcwd(), "DEC_" + os.path.basename(tmp_path))
            if os.path.exists(cwd_dec_path):
                os.remove(cwd_dec_path)

        info = {
            'valid_packets': success_count,
            'corrupt_packets': corrupt_count,
            'data_length': len(binary_data),
            'valid': valid,
            'total_bits': total_bits,
        }

    except Exception as e:
        if logger:
            logger.error(f"DNA-Aeon decode error: {e}")
            logger.error(traceback.format_exc())
        binary_data = ""
        valid = False
        info = {'error': str(e)}

    return binary_data, valid, info
