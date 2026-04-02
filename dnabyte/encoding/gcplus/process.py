import os
from collections import Counter, defaultdict


def process(data, params, logger=None):
    """
    Process DNA strands after synthesis/sequencing for GC+ encoding.

    GC+ produces one oligo per *k*-bit chunk.  After synthesis each oligo
    is copied many times; sequencing may introduce per-copy errors.

    Processing:
    1. Group identical (or near-identical) sequences together.
    2. Majority-vote within each group to reconstruct the original oligo.
    3. Return one consensus sequence per group.

    GC+ does not use primers by default — processing operates on full
    sequences.

    Args:
        data:   Data object with ``data.data`` = list of DNA strings.
        params: Parameters object.
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
            return [], {}

        # Clean
        dna_strands = [seq.replace(' ', '').strip() for seq in dna_strands]

        if logger:
            logger.info(f"Processing {total_count} GC+ DNA sequences")

        # ----- Group identical sequences ------------------------------------
        groups = defaultdict(list)
        for seq in dna_strands:
            groups[seq].append(seq)

        if logger:
            logger.info(f"Found {len(groups)} unique sequence groups")

        # ----- Separate multi-copy groups from singletons -------------------
        # After synthesis (many copies) + sequencing (errors on some copies):
        #   - Correct copies are identical → form large groups
        #   - Corrupted copies are unique  → singletons
        large_groups = {s: c for s, c in groups.items() if len(c) >= 2}
        singletons   = {s: c for s, c in groups.items() if len(c) == 1}

        if logger:
            logger.info(
                f"Multi-copy groups: {len(large_groups)}, "
                f"singletons: {len(singletons)}"
            )

        if large_groups:
            consensus_sequences = []
            for representative, copies in large_groups.items():
                consensus = _majority_vote(copies)
                consensus_sequences.append(consensus)
            if logger:
                logger.info(
                    f"Filtered {len(singletons)} singleton sequences "
                    "(likely corrupted)"
                )
        else:
            # No multi-copy groups — keep everything (e.g. no synthesis copies)
            consensus_sequences = list(groups.keys())

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
        }

        return consensus_sequences, info

    except Exception as e:
        if logger:
            logger.error(f"Error processing DNA strands: {e}")
            import traceback
            logger.error(traceback.format_exc())
        return None, {}


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
