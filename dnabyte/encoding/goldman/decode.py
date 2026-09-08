import tempfile
import os
import traceback

from dnabyte.encoding.goldman.StorageD.goldmanDecode import (
   goldmanDecode
)


def decode(data, params, logger=None):
    """
    Decode Goldman DNA sequences.

    Goldman decoding:
        DNA
          ↓
        decodeNt()
          ↓
        index + payload
          ↓
        overlap-aware reconstruction
          ↓
        ternary → bytes
          ↓
        original bitstream
    """

    binary_data = ""
    valid = False
    info = {}

    try:
        if logger:
            logger.info("Starting Goldman decode")

        # --------------------------------------------------
        # 1. Get DNA sequences
        # --------------------------------------------------

        if isinstance(data.data, list):
            dna_sequences = data.data
        else:
            dna_sequences = [data.data]

        clean_sequences = [
            str(seq).replace(" ", "").strip().upper()
            for seq in dna_sequences
            if str(seq).replace(" ", "").strip()
        ]

        if not clean_sequences:
            if logger:
                logger.warning("No DNA sequences available")

            return "", False, {
                "error": "No DNA sequences available"
            }

        # --------------------------------------------------
        # 2. Read Goldman metadata
        # --------------------------------------------------

        original_metadata = getattr(
            params,
            "goldman_fasta_metadata",
            ""
        )

        if not original_metadata:
            if logger:
                logger.error(
                    "Missing Goldman FASTA metadata"
                )

            return "", False, {
                "error": "Missing Goldman FASTA metadata"
            }

        # --------------------------------------------------
        # 3. Parse metadata
        # --------------------------------------------------

        meta_clean = original_metadata.lstrip(">").strip()

        meta_parts = {}

        for part in meta_clean.split(","):

            if ":" not in part:
                continue

            key, val = part.split(":", 1)

            meta_parts[key.strip()] = val.strip()

        try:
            index_length = int(meta_parts["indexLen"])
            add_len = int(meta_parts["addLen"])

        except (KeyError, ValueError) as e:

            return "", False, {
                "error": f"Invalid Goldman metadata: {e}"
            }

        # --------------------------------------------------
        # 4. Original length
        # --------------------------------------------------

        total_bits = int(
            getattr(
                params,
                "goldman_total_bits",
                0
            )
        )

        if logger:
            logger.info(
                f"Goldman parameters: "
                f"indexLen={index_length}, "
                f"addLen={add_len}"
            )

            logger.info(
                f"Received {len(clean_sequences)} DNA sequences"
            )

        # --------------------------------------------------
        # 5. Call the NEW Goldman decoder
        #
        # IMPORTANT:
        # The substitution-error correction should happen
        # inside goldmanDecode().
        # --------------------------------------------------

        temp_output_dir = tempfile.mkdtemp()

        decode_output_file = os.path.join(
            temp_output_dir,
            "goldman_decoded.bin"
        )

        try:

            byte_list, decode_info = goldmanDecode(
                nt_list=clean_sequences,
                save_path=decode_output_file,
                idnex_len=index_length,
                add_len=add_len
            )

            # --------------------------------------------------
            # 6. Read reconstructed bytes
            # --------------------------------------------------

            if not os.path.exists(decode_output_file):

                return "", False, {
                    "error": "Goldman decoder produced no output"
                }

            with open(
                decode_output_file,
                "rb"
            ) as f:

                file_content = f.read()

            # --------------------------------------------------
            # 7. Convert bytes → bitstream
            # --------------------------------------------------

            binary_data = "".join(
                format(b, "08b")
                for b in file_content
            )

            if total_bits > 0:

                binary_data = binary_data[:total_bits]

            # --------------------------------------------------
            # 8. Validate output length
            # --------------------------------------------------

            length_ok = True

            if total_bits > 0:

                length_ok = (
                    len(binary_data) == total_bits
                )

            # --------------------------------------------------
            # 9. Determine validity
            # --------------------------------------------------

            decoder_valid = decode_info.get(
                "valid",
                True
            )

            valid = (
                decoder_valid
                and length_ok
            )

            # --------------------------------------------------
            # 10. Return information
            # --------------------------------------------------

            info = {
                **decode_info,

                "decoded_sequences":
                    decode_info.get(
                        "decoded_sequences",
                        0
                    ),

                "skipped_sequences":
                    decode_info.get(
                        "skipped_sequences",
                        0
                    ),

                "data_length":
                    len(binary_data),

                "total_bits":
                    total_bits,

                "index_length":
                    index_length,

                "add_len":
                    add_len,

                "length_valid":
                    length_ok,

                "valid":
                    valid,
            }

            if logger:

                logger.info(
                    f"Goldman decode finished: "
                    f"{len(binary_data)} bits"
                )

                logger.info(
                    f"Goldman validity: {valid}"
                )

        finally:

            # --------------------------------------------------
            # 11. Cleanup
            # --------------------------------------------------

            try:

                if os.path.exists(
                    decode_output_file
                ):
                    os.remove(
                        decode_output_file
                    )

                if os.path.exists(
                    temp_output_dir
                ):
                    os.rmdir(
                        temp_output_dir
                    )

            except Exception as cleanup_error:

                if logger:

                    logger.warning(
                        f"Cleanup failed: "
                        f"{cleanup_error}"
                    )

    except Exception as e:

        if logger:

            logger.error(
                f"Error during Goldman decode: {e}"
            )

            logger.error(
                traceback.format_exc()
            )

        binary_data = ""
        valid = False

        info = {
            "error": str(e),
            "valid": False
        }

    return binary_data, valid, info