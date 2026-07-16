import os
from collections import Counter, defaultdict


def process(data, params, logger=None):
    """
    Processes DNA strands after synthesis/sequencing simulation for Church encoding.

    Groups duplicate/similar sequences and performs majority voting per
    nucleotide position to reconstruct the most likely original sequences.

    Approach:
    1. Strip primers from each sequence to get the payload.
    2. Normalize payload lengths (trim or pad to expected length).
    3. Group payloads by similarity (cluster around unique originals).
    4. For each group, perform position-wise majority voting.
    5. Re-attach primers and return the consensus sequences.

    Args:
        data: Data object containing DNA sequences in data.data
        params: Parameters object with primer_length, sequence_length, etc.
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

        primer_length = int(getattr(params, 'primer_length', 20))
        add_primer = bool(getattr(params, 'add_primer', True))
        sequence_length = int(getattr(params, 'sequence_length', 200))

        # Expected payload length (sequence without primers)
        expected_payload_len = sequence_length - (2 * primer_length) if add_primer else sequence_length

        if add_primer and primer_length > 0:
            left_primer, right_primer = _extract_primers(dna_strands, primer_length)

            # Strip primers to get payloads
            payloads = []
            for seq in dna_strands:
                payload = seq[primer_length:-primer_length] if primer_length > 0 else seq
                payloads.append(payload)
        else:
            left_primer = ""
            right_primer = ""
            payloads = list(dna_strands)

        if logger:
            logger.info(f"Processing {total_count} DNA sequences")
            logger.info(f"Expected payload length: {expected_payload_len}")

        # Normalize payload lengths to the expected length
        normalized_payloads = _normalize_lengths(payloads, expected_payload_len)

        # Group identical payloads together
        groups = defaultdict(list)
        for payload in normalized_payloads:
            groups[payload].append(payload)

        if logger:
            logger.info(f"Found {len(groups)} unique payload groups after normalization")

        # Cluster and vote to find consensus for each original sequence
        num_original = _estimate_num_originals(groups, total_count, logger)
        consensus_payloads = _cluster_and_vote(groups, num_original, logger)

        # Re-attach primers
        consensus_sequences = []
        for payload in consensus_payloads:
            consensus_sequences.append(left_primer + payload + right_primer)

        if logger:
            logger.info(f"Consensus: {len(consensus_sequences)} sequences from {total_count} inputs")
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


def _extract_primers(dna_strands, primer_length):
    """
    Extract the most common left and right primer from the reads.
    """
    left_counts = Counter(seq[:primer_length] for seq in dna_strands)
    right_counts = Counter(seq[-primer_length:] for seq in dna_strands)
    left_primer = left_counts.most_common(1)[0][0]
    right_primer = right_counts.most_common(1)[0][0]
    return left_primer, right_primer


def _normalize_lengths(payloads, expected_length):
    """
    Normalize all payloads to the expected length.
    """
    normalized = []
    for payload in payloads:
        if len(payload) == expected_length:
            normalized.append(payload)
        elif len(payload) > expected_length:
            normalized.append(payload[:expected_length])
        else:
            normalized.append(payload + 'N' * (expected_length - len(payload)))
    return normalized


def _estimate_num_originals(groups, total_count, logger=None):
    """
    Estimate the number of original unique sequences before errors.

    If every group is a singleton (max group size == 1) it means either
    - mean=1 was used (one copy per oligo, no synthesis redundancy), or
    - sequencing errors are so severe that all copies look different.
    In both cases we cannot collapse groups, so we treat every unique
    sequence as its own original.
    """
    sizes = sorted([len(copies) for copies in groups.values()], reverse=True)

    if not sizes:
        return 0

    # If all groups are singletons, keep them all — no basis for merging.
    if sizes[0] == 1:
        estimated = len(sizes)
        if logger:
            logger.info(f"Group sizes (top 10): {sizes[:10]}")
            logger.info(f"All groups are singletons — treating each as a distinct original")
            logger.info(f"Estimated {estimated} original sequences")
        return estimated

    avg_if_3 = total_count / 3
    significant = [s for s in sizes if s >= max(2, avg_if_3 * 0.1)]
    estimated = max(1, len(significant))

    if logger:
        logger.info(f"Group sizes (top 10): {sizes[:10]}")
        logger.info(f"Estimated {estimated} original sequences")

    return estimated


def _cluster_and_vote(groups, num_expected, logger=None):
    """
    Cluster similar payloads and perform majority voting within each cluster.
    """
    sorted_groups = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)

    if not sorted_groups:
        return []

    num_centers = min(num_expected, len(sorted_groups))
    centers = {}
    for i in range(num_centers):
        key, copies = sorted_groups[i]
        centers[key] = list(copies)

    remaining = sorted_groups[num_centers:]
    orphan_count = 0
    for payload, copies in remaining:
        center_keys = list(centers.keys())
        best_center = _find_closest(payload, center_keys)
        centers[best_center].extend(copies)
        orphan_count += len(copies)

    consensus_payloads = []
    for center_key, all_copies in centers.items():
        consensus = _majority_vote(all_copies)
        consensus_payloads.append(consensus)

    if logger:
        logger.info(f"Clustering: {num_centers} clusters, "
                     f"{orphan_count} orphans merged")

    return consensus_payloads


def _find_closest(query, candidates):
    """
    Find the candidate string with the smallest distance to query.
    """
    best = candidates[0]
    best_dist = float('inf')

    for candidate in candidates:
        min_len = min(len(query), len(candidate))
        dist = sum(a != b for a, b in zip(query[:min_len], candidate[:min_len]))
        dist += abs(len(query) - len(candidate))
        if dist < best_dist:
            best_dist = dist
            best = candidate

    return best


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
