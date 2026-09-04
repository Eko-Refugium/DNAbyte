import os
from collections import Counter, defaultdict
from dnabyte.data_classes.insilicodna import InSilicoDNA


def process(data, params, logger=None):
    """
    Process DNA strands after synthesis/sequencing for GC+ encoding.

    GC+ produces one oligo per *k*-bit chunk.  After synthesis each oligo
    is copied many times; sequencing may introduce per-copy errors.

    Processing:
    1. Extract position tags from sequences if present (variable-length safe).
    2. Group identical (or near-identical) sequences together.
    3. Majority-vote within each group to reconstruct the original oligo.
    4. Return one consensus sequence per group, ordered by position.

    GC+ does not use primers by default — processing operates on full
    sequences.

    Args:
        data:   Data object with ``data.data`` = list of DNA strings.
        params: Parameters object with optional position tag info.
        logger: Optional logger.

    Returns:
        (consensus_sequences_list, info_dict)
    """
    try:
        dna_strands = data.data
        total_count = len(dna_strands) if dna_strands else 0

        if not dna_strands:
            if logger:
                logger.warning("No DNA sequences to process")
            return InSilicoDNA([]), {}

        # Clean
        dna_strands = [seq.replace(' ', '').strip() for seq in dna_strands]

        if logger:
            logger.info(f"Processing {total_count} GC+ DNA sequences")

        # Extract position tags if embedded (variable-length safe)
        pos_bits = int(getattr(params, 'gcplus_position_bits', 0))
        tag_redundancy = int(getattr(params, 'gcplus_tag_redundancy', 1))
        tag_length = int(getattr(params, 'gcplus_tag_length', 0))
        
        if pos_bits == 0:
            tag_length = 0
        elif tag_length == 0:
            tag_length = pos_bits * tag_redundancy
        
        sequence_positions = []  # Store (position, stripped_sequence, original_sequence)
        
        for idx, seq in enumerate(dna_strands):
            position = idx
            stripped_seq = seq
            
            if tag_length > 0 and len(seq) >= tag_length:
                # Extract tag from end (works with variable-length sequences)
                tag = seq[-tag_length:]
                stripped_seq = seq[:-tag_length]
                
                # Decode redundant tag with majority voting
                pos_binary_bits = []
                for i in range(0, tag_length, tag_redundancy):
                    tag_segment = tag[i:i+tag_redundancy]
                    a_count = tag_segment.count('A')
                    t_count = tag_segment.count('T')
                    bit = '0' if a_count >= t_count else '1'
                    pos_binary_bits.append(bit)
                
                pos_binary = ''.join(pos_binary_bits[:pos_bits])
                try:
                    position = int(pos_binary, 2)
                except ValueError:
                    position = idx
            
            sequence_positions.append((position, stripped_seq, seq))
        
        # ----- Group identical sequences (by stripped sequence) -----------
        groups = defaultdict(list)
        for position, stripped_seq, original_seq in sequence_positions:
            groups[stripped_seq].append((position, original_seq))

        if logger:
            logger.info(f"Found {len(groups)} unique sequence groups")

        # ----- Separate multi-copy groups from singletons -------------------
        large_groups = {s: c for s, c in groups.items() if len(c) >= 2}
        singletons   = {s: c for s, c in groups.items() if len(c) == 1}

        if logger:
            logger.info(
                f"Multi-copy groups: {len(large_groups)}, "
                f"singletons: {len(singletons)}"
            )

        # Consensus: always keep large_groups + singletons
        consensus_sequences = []
        position_sequence_list = []
        
        # Add consensus from multi-copy groups
        for representative, position_copies in large_groups.items():
            copies = [seq for pos, seq in position_copies]
            consensus = _majority_vote(copies)
            # Use the position from the first copy
            position = position_copies[0][0]
            consensus_sequences.append(consensus)
            position_sequence_list.append((position, consensus))
        
        # Always keep singletons
        for singleton_seq, position_copies in singletons.items():
            consensus_sequences.append(singleton_seq)
            position = position_copies[0][0]
            position_sequence_list.append((position, singleton_seq))

        if logger:
            if singletons:
                logger.info(
                    f"Kept {len(singletons)} singleton sequences "
                    "(may be valid or filtered by downstream)"
                )

        # Reorder by position if tags were present
        if tag_length > 0 and position_sequence_list:
            sorted_sequences = [seq for pos, seq in sorted(position_sequence_list, key=lambda x: x[0])]
            if logger:
                logger.info(f"Reordered {len(sorted_sequences)} sequences by position tag")
            consensus_sequences = sorted_sequences

        if logger:
            logger.info(
                f"Consensus: {len(consensus_sequences)} unique sequences "
                f"from {total_count} inputs"
            )

        info = {
            'number_of_sequences_input': total_count,
            'number_of_sequences_output': len(consensus_sequences),
            'unique_groups': len(groups),
            'duplicates_removed': total_count - len(consensus_sequences),
            'status': 'consensus',
            'position_tags_extracted': tag_length > 0,
        }

        return InSilicoDNA(consensus_sequences), info

    except Exception as e:
        if logger:
            logger.error(f"Error processing DNA strands: {e}")
            import traceback
            logger.error(traceback.format_exc())
        return InSilicoDNA([]), {}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _majority_vote(sequences):
    """Position-wise majority vote across a list of DNA sequences."""
    if not sequences:
        return ""
    if len(sequences) == 1:
        return sequences[0]

    lengths = Counter(len(s) for s in sequences)
    target_len = lengths.most_common(1)[0][0]

    consensus = []
    for i in range(target_len):
        bases_at_pos = [s[i] for s in sequences if i < len(s) and s[i] in 'ACGT']
        if bases_at_pos:
            most_common = Counter(bases_at_pos).most_common(1)[0][0]
            consensus.append(most_common)
        elif any(i < len(s) for s in sequences):
            consensus.append('A')

    return ''.join(consensus)
