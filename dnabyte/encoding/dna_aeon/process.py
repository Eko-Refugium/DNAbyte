import os
from collections import Counter, defaultdict


def process(data, params, logger=None):
    """
    Process DNA strands after synthesis/sequencing for DNA-Aeon encoding.

    DNA-Aeon uses NOREC4DNA fountain codes where each oligo is a
    self-identifying packet (carries its own ID and chunk-mapping in a
    header).  The fountain-code decoder handles corruption detection
    internally — ``parse_raw_packet`` returns ``"CORRUPT"`` for packets
    whose headers are damaged.

    Processing strategy for fountain codes:
    1. Group identical sequences to de-duplicate synthesis copies.
    2. For each group with ≥ 2 copies, take a majority-vote consensus.
    3. **Keep singletons as well** — unlike fixed-position codes, every
       unique packet potentially carries new information.  Discarding
       singletons risks losing valid packets that the decoder needs.
    4. Return ALL unique/consensus sequences and let the fountain-code
       decoder decide which are usable.

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
            logger.info(f"Processing {total_count} DNA-Aeon DNA sequences")

        # ----- Group identical sequences ------------------------------------
        groups = defaultdict(list)
        for seq in dna_strands:
            groups[seq].append(seq)

        if logger:
            logger.info(f"Found {len(groups)} unique sequence groups")

        # ----- Separate multi-copy groups from singletons -------------------
        large_groups = {s: c for s, c in groups.items() if len(c) >= 2}
        singletons  = {s: c for s, c in groups.items() if len(c) == 1}

        if logger:
            logger.info(
                f"Multi-copy groups: {len(large_groups)}, "
                f"singletons: {len(singletons)}"
            )

        # ----- Build consensus list -----------------------------------------
        # For multi-copy groups: majority-vote to get best consensus
        consensus_sequences = []
        for representative, copies in large_groups.items():
            consensus = _majority_vote(copies)
            consensus_sequences.append(consensus)

        # For fountain codes: KEEP singletons — they may carry unique
        # packet IDs that the decoder needs.  The NOREC4DNA decoder will
        # reject any packet whose header is corrupted, so there is no
        # harm in passing them through.
        for seq in singletons:
            consensus_sequences.append(seq)

        if logger:
            logger.info(
                f"Consensus: {len(consensus_sequences)} unique sequences "
                f"from {total_count} inputs "
                f"({len(large_groups)} voted, {len(singletons)} singletons kept)"
            )

        info = {
            'number_of_sequences_input': total_count,
            'number_of_sequences_output': len(consensus_sequences),
            'unique_groups': len(groups),
            'multi_copy_groups': len(large_groups),
            'singletons_kept': len(singletons),
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
