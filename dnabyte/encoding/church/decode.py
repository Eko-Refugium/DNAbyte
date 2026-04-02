import tempfile
import os
import traceback

from dnabyte.encoding.church.StorageD.codec import ChurchDecode


def decode(data, params, logger=None):
    """
    Decodes DNA sequences from Church encoding.
    Converts DNA sequences to FASTA file, uses ChurchDecode, reads output binary,
    and cleans up temporary files.

    Church uses the same AbstractDecode pipeline as Wukong, with the same
    FASTA metadata format: totalBit, binSegLen, leftPrimer, rightPrimer,
    fileExtension, bRedundancy, RSNum.

    Args:
        data: Data object containing DNA sequences in data.data (list of strings)
        params: Parameters object with decoding settings
        logger: Optional logger

    Returns:
        tuple: (binary data string, validity flag, info dictionary)
    """
    temp_fasta_file = None
    temp_output_dir = None
    decode_output_file = None
    binary_data = ""
    valid = False
    info = {}

    def _get_total_bits():
        """Attempt to determine total bits from file or params."""
        try:
            file_paths = getattr(data, 'file_paths', None)
            if file_paths and len(file_paths) > 0 and os.path.exists(file_paths[0]):
                return os.path.getsize(file_paths[0]) * 8
        except Exception:
            pass
        if hasattr(params, 'total_bits'):
            return int(getattr(params, 'total_bits'))
        return 0

    try:
        # Create temporary output directory
        temp_output_dir = tempfile.mkdtemp()

        if logger:
            logger.info(f"Temp output directory: {temp_output_dir}")

        # Create FASTA file from DNA sequences
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        temp_fasta_file = temp_fasta.name

        if logger:
            logger.info(f"Creating FASTA file for decoding: {temp_fasta_file}")

        dna_sequences = data.data if isinstance(data.data, list) else [data.data]

        # Clean sequences
        clean_sequences = [str(seq).replace(' ', '').strip() for seq in dna_sequences if str(seq).replace(' ', '').strip()]
        if not clean_sequences:
            if logger:
                logger.warning("No DNA sequences available for decoding")
            return "", False, {}

        # Use the original FASTA metadata header preserved from encoding
        original_metadata = getattr(params, 'church_fasta_metadata', '')

        if original_metadata:
            meta_line = original_metadata if original_metadata.startswith('>') else '>' + original_metadata
            meta_line = meta_line.rstrip('\n') + '\n'
            if logger:
                logger.info(f"Using original FASTA metadata from encoding: {meta_line.strip()}")
        else:
            # Fallback: reconstruct metadata (less reliable)
            if logger:
                logger.warning("No original FASTA metadata found, reconstructing")
            total_bits = _get_total_bits()
            primer_length = int(getattr(params, 'primer_length', 20))
            left_primer = clean_sequences[0][:primer_length] if len(clean_sequences[0]) >= primer_length else ""
            right_primer = clean_sequences[0][-primer_length:] if len(clean_sequences[0]) >= primer_length else ""
            # Church uses 1:1 bit-to-base ratio, index_redundancy=0
            bin_seg_len = max(1, int(getattr(params, 'sequence_length', 200)) - (2 * primer_length))
            add_redundancy = 1 if bool(getattr(params, 'add_redundancy', False)) else 0
            rs_num = int(getattr(params, 'rs_num', 0))
            meta_line = (
                f">totalBit:{total_bits},binSegLen:{bin_seg_len},leftPrimer:{left_primer},"
                f"rightPrimer:{right_primer},fileExtension:.bin,"
                f"bRedundancy:{add_redundancy},RSNum:{rs_num}\n"
            )

        temp_fasta.write(meta_line)

        # Extract total_bits from metadata for later truncation
        total_bits = 0
        try:
            total_bits = int(meta_line.split('totalBit:')[1].split(',')[0])
        except Exception:
            total_bits = _get_total_bits()

        # Write sequences in FASTA format
        for i, sequence in enumerate(clean_sequences):
            clean_sequence = str(sequence).replace(' ', '').strip()
            if clean_sequence:
                temp_fasta.write(f">seq_{i + 1}\n{clean_sequence}\n")

        temp_fasta.close()

        if logger:
            logger.info(f"FASTA file created with {len(clean_sequences)} sequences")

        # Perform Church decoding
        if logger:
            logger.info(f"Starting ChurchDecode with FASTA: {temp_fasta_file}")

        decode_worker = ChurchDecode(
            input_file_path=temp_fasta_file,
            output_dir=temp_output_dir,
        )

        decode_output_file = decode_worker.common_decode()

        if logger:
            logger.info(f"Decoding completed. Output file: {decode_output_file}")

        # Read the decoded binary data
        if decode_output_file and os.path.exists(decode_output_file):
            file_size = os.path.getsize(decode_output_file)
            if logger:
                logger.info(f"Output file size: {file_size} bytes")

            with open(decode_output_file, 'rb') as f:
                file_content = f.read()

            if file_content:
                # Convert output bytes to bitstream expected by BinaryCode
                binary_data = ''.join(format(b, '08b') for b in file_content)
                if total_bits > 0:
                    binary_data = binary_data[:total_bits]

                if logger:
                    logger.info(f"Binary data read successfully. Length: {len(binary_data)} bits")
                valid = True
            else:
                if logger:
                    logger.warning("Output file is empty")
                valid = False
        else:
            if logger:
                logger.error(f"Decode output file missing: {decode_output_file}")
            valid = False

        # Create info dictionary
        info = {
            'decoded_file': decode_output_file,
            'data_length': len(binary_data) if binary_data else 0,
            'valid': valid,
            'total_bits': total_bits,
            'metadata_used': meta_line.strip()
        }

    except Exception as e:
        if logger:
            logger.error(f"Error during decoding: {str(e)}")
            logger.error(traceback.format_exc())

        binary_data = ""
        valid = False
        info = {'error': str(e)}

    finally:
        # Clean up temporary files
        for filepath in [temp_fasta_file, decode_output_file]:
            if filepath and os.path.exists(filepath):
                try:
                    os.remove(filepath)
                    if logger:
                        logger.info(f"Temporary file deleted: {filepath}")
                except Exception as e:
                    if logger:
                        logger.warning(f"Could not delete file {filepath}: {str(e)}")

        if temp_output_dir and os.path.exists(temp_output_dir):
            try:
                os.rmdir(temp_output_dir)
                if logger:
                    logger.info(f"Temporary output directory deleted: {temp_output_dir}")
            except Exception as e:
                if logger:
                    logger.warning(f"Could not delete directory: {str(e)}")

    return binary_data, valid, info
