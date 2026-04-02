import tempfile
import os
import traceback

from dnabyte.encoding.wukong.StorageD.goldmanDecode import (
    decodeNt, combineHuffman, huffmanToByte, saveResult
)


def decode(data, params, logger=None):
    """
    Decodes DNA sequences from Goldman encoding.
    Converts DNA sequences to FASTA file, uses GoldmanDecode, reads output binary,
    and cleans up temporary files.

    Goldman uses a different metadata format than Church/Wukong:
    >indexLen:{},addLen:{},fileExtension:{}

    Goldman overrides common_decode() entirely and uses goldmanDecode() which
    writes the output file directly.

    Args:
        data: Data object containing DNA sequences in data.data (list of strings)
        params: Parameters object with decoding settings
        logger: Optional logger

    Returns:
        tuple: (binary data string, validity flag, info dictionary)
    """
    binary_data = ""
    valid = False
    info = {}

    try:
        if logger:
            logger.info("Starting Goldman decode")

        dna_sequences = data.data if isinstance(data.data, list) else [data.data]

        # Clean sequences
        clean_sequences = [str(seq).replace(' ', '').strip() for seq in dna_sequences if str(seq).replace(' ', '').strip()]
        if not clean_sequences:
            if logger:
                logger.warning("No DNA sequences available for decoding")
            return "", False, {}

        # Use the original Goldman FASTA metadata preserved from encoding
        original_metadata = getattr(params, 'goldman_fasta_metadata', '')
        if not original_metadata:
            if logger:
                logger.error("No Goldman FASTA metadata found - cannot decode without indexLen and addLen")
            return "", False, {'error': 'Missing Goldman FASTA metadata'}

        # Parse indexLen and addLen from the metadata string
        meta_clean = original_metadata.lstrip('>').strip()
        meta_parts = {}
        for part in meta_clean.split(','):
            key, val = part.split(':')
            meta_parts[key] = val

        index_length = int(meta_parts['indexLen'])
        add_len = int(meta_parts['addLen'])

        if logger:
            logger.info(f"Goldman params: indexLen={index_length}, addLen={add_len}")
            logger.info(f"Decoding {len(clean_sequences)} DNA sequences")

        # Get total_bits for later truncation (Goldman metadata doesn't include totalBit)
        total_bits = int(getattr(params, 'goldman_total_bits', 0))

        # Decode each sequence with fault tolerance
        # Goldman's rotating code crashes on homopolymers (from sequencing errors).
        # Skip sequences that can't be decoded — the 4x overlap provides redundancy.
        huffman_str_list = []
        skipped = 0
        for nt_seq in clean_sequences:
            try:
                huffman_str = decodeNt(nt_seq)
                huffman_str_list.append(huffman_str)
            except (ValueError, KeyError, IndexError):
                skipped += 1
                continue

        if logger:
            logger.info(f"Successfully decoded {len(huffman_str_list)} sequences, skipped {skipped} corrupted")

        if not huffman_str_list:
            if logger:
                logger.error("No sequences could be decoded")
            return "", False, {'error': 'All sequences corrupted'}

        # Combine overlapping Huffman segments
        total_huffman_str = combineHuffman(huffman_str_list, index_length, add_len)

        # Convert Huffman to bytes
        byte_list = huffmanToByte(total_huffman_str)

        if not byte_list:
            if logger:
                logger.error("No bytes decoded from Huffman string")
            return "", False, {'error': 'Empty byte output'}

        # Save to temp file, read back as binary string
        temp_output_dir = tempfile.mkdtemp()
        decode_output_file = os.path.join(temp_output_dir, 'goldman_decoded.bin')
        saveResult(byte_list, decode_output_file)

        with open(decode_output_file, 'rb') as f:
            file_content = f.read()

        # Convert output bytes to bitstream
        binary_data = ''.join(format(b, '08b') for b in file_content)
        if total_bits > 0:
            binary_data = binary_data[:total_bits]

        if logger:
            logger.info(f"Binary data read successfully. Length: {len(binary_data)} bits")

        valid = True

        # Create info dictionary
        info = {
            'decoded_sequences': len(huffman_str_list),
            'skipped_sequences': skipped,
            'data_length': len(binary_data),
            'valid': valid,
            'total_bits': total_bits,
            'index_length': index_length,
            'add_len': add_len,
        }

        # Clean up
        try:
            os.remove(decode_output_file)
            os.rmdir(temp_output_dir)
        except Exception:
            pass

    except Exception as e:
        if logger:
            logger.error(f"Error during decoding: {str(e)}")
            logger.error(traceback.format_exc())

        binary_data = ""
        valid = False
        info = {'error': str(e)}

    return binary_data, valid, info
