import os
from collections import Counter, defaultdict


def process(data, params, logger=None):
    """
    Processes DNA strands after synthesis/sequencing simulation for Goldman encoding.

    Goldman encoding produces many UNIQUE sequences, each with a different index
    and overlapping content. Unlike Church/Wukong, Goldman sequences are NOT
    copies of a small number of originals — each is distinct and needed for decoding.
    Goldman's error tolerance comes from its 4x overlapping segment design, not
    from having identical copies.

    Therefore, processing for Goldman:
    1. Groups identical (or near-identical) sequences together.
    2. Performs majority voting WITHIN each group (to fix copy errors).
    3. Returns ALL unique groups as separate sequences (does NOT merge groups).

    Goldman does NOT use primers, so processing operates on full sequences.

    Args:
        data: Data object containing DNA sequences in data.data
        params: Parameters object with sequence_length, etc.
        logger: Optional logger for info/error messages

    Returns:
        tuple: (consensus DNA sequences list, info dictionary)
    """
    try:
        dna_strands = data.data
        total_count = len(dna_strands) if dna_strands else 0

        if not dna_strands:
            if logger:
                logger.warning("No DNA sequences to process")
            return [], {}

        # Clean sequences
        dna_strands = [seq.replace(' ', '').strip() for seq in dna_strands]

        sequence_length = int(getattr(params, 'sequence_length', 200))

        if logger:
            logger.info(f"Processing {total_count} DNA sequences (Goldman)")

        # Goldman sequences have a specific length determined by the algorithm.
        # Do NOT normalize lengths — Goldman's rotating code is sensitive to
        # any modification of the sequence content.
        # Just group identical sequences and vote within each group.
        groups = defaultdict(list)
        for seq in dna_strands:
            groups[seq].append(seq)

        if logger:
            logger.info(f"Found {len(groups)} unique sequence groups")

        # Separate multi-copy groups from singletons.
        # After synthesis (many copies) + sequencing (errors on some copies):
        # - Correct copies are identical → form large groups
        # - Corrupted copies are unique → form singletons
        # Singletons are almost certainly corrupted and will crash Goldman's
        # rotating code decoder (homopolymers cause ValueError).
        large_groups = {seq: copies for seq, copies in groups.items() if len(copies) >= 2}
        singletons = {seq: copies for seq, copies in groups.items() if len(copies) == 1}

        if logger:
            logger.info(f"Multi-copy groups: {len(large_groups)}, singletons: {len(singletons)}")

        if large_groups:
            # Use only multi-copy groups — singletons are likely corrupted
            consensus_sequences = []
            for representative, copies in large_groups.items():
                consensus = _majority_vote(copies)
                consensus_sequences.append(consensus)
            if logger:
                logger.info(f"Filtered {len(singletons)} singleton sequences (likely corrupted)")
        else:
            # No multi-copy groups — keep everything (no synthesis copies available)
            consensus_sequences = list(groups.keys())

        if logger:
            logger.info(f"Consensus: {len(consensus_sequences)} unique sequences from {total_count} inputs")
            if consensus_sequences:
                logger.info(f"DNA sequence length: {len(consensus_sequences[0])}")

        info = {
            'number_of_sequences_input': total_count,
            'number_of_sequences_output': len(consensus_sequences),
            'unique_groups': len(groups),
            'duplicates_removed': total_count - len(consensus_sequences),
            'status': 'consensus'
        }

        return consensus_sequences, info

    except Exception as e:
        if logger:
            logger.error(f"Error processing DNA strands: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
        return None, {}


def _normalize_lengths(sequences, expected_length):
    """
    Normalize all sequences to the expected length.
    """
    normalized = []
    for seq in sequences:
        if len(seq) == expected_length:
            normalized.append(seq)
        elif len(seq) > expected_length:
            normalized.append(seq[:expected_length])
        else:
            normalized.append(seq + 'N' * (expected_length - len(seq)))
    return normalized


def _majority_vote(sequences):
    """
    Position-wise majority vote across a list of DNA sequences.
    Ignores 'N' padding characters in the vote.
    """
    if not sequences:
        return ""
    if len(sequences) == 1:
        return sequences[0]

    lengths = Counter(len(s) for s in sequences)
    target_len = lengths.most_common(1)[0][0]

    consensus = []
    for i in range(target_len):
        bases_at_pos = [s[i] for s in sequences if i < len(s) and s[i] != 'N']
        if bases_at_pos:
            most_common = Counter(bases_at_pos).most_common(1)[0][0]
            consensus.append(most_common)
        elif any(i < len(s) for s in sequences):
            consensus.append('A')

    return ''.join(consensus)
