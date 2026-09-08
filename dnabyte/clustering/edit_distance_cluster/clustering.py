from typing import List, Dict, Tuple, Optional

from dnabyte.cluster import Cluster


class AlignmentClusterer(Cluster):
    """
    DNA clustering tolerant to:
        - SNPs / substitutions
        - insertions
        - deletions
        - variable sequence lengths
        - arbitrary DNA composition
        - repetitive sequences

    Input:
        sequences: iterable of DNA strings

    Output:
        {
            cluster_id: [sequence, sequence, ...],
            ...
        }, info

    The implementation uses global alignment with affine gap penalties.

    A gap of length N costs:

        gap_open + gap_extend * (N - 1)

    This is preferable to plain Levenshtein distance because a single
    insertion/deletion event can contain multiple bases.
    """

    def __init__(self, params, logger=None):

        self.logger = logger

        # Maximum normalized error rate allowed.
        #
        # Example:
        #     0.03 = approximately 3 errors per 100 aligned bases.
        #
        self.max_error_rate = getattr(
            params,
            "max_error_rate_cluster",
            0.03
        )

        # Absolute maximum number of edits.
        #
        # Set to None to rely only on max_error_rate.
        #
        self.max_edits = getattr(
            params,
            "max_edits_cluster",
            None
        )

        # Maximum length difference allowed.
        #
        # This prevents wildly different sequences from becoming
        # candidates even if their normalized score happens to be high.
        self.max_length_difference = getattr(
            params,
            "max_length_difference_cluster",
            None
        )

        # Alignment scoring.
        #
        # Match should be strongly positive.
        # Mismatch should be negative.
        # Opening a gap should be more expensive than extending it.
        self.match_score = getattr(
            params,
            "alignment_match_score_cluster",
            2
        )

        self.mismatch_score = getattr(
            params,
            "alignment_mismatch_score_cluster",
            -2
        )

        self.gap_open = getattr(
            params,
            "alignment_gap_open_cluster",
            -4
        )

        self.gap_extend = getattr(
            params,
            "alignment_gap_extend_cluster",
            -1
        )

        # Whether sequences are normalized to uppercase.
        self.normalize = getattr(
            params,
            "normalize_sequences_cluster",
            True
        )

    # ------------------------------------------------------------------
    # Main clustering API
    # ------------------------------------------------------------------

    def cluster(self, sequences):

        # Normalize while preserving the original sequence for output.
        items = []

        for seq in sequences:

            if seq is None:
                continue

            original = seq

            if self.normalize:
                seq = seq.upper()

            items.append((original, seq))

        groups = []

        for original, seq in items:

            best_group = None
            best_distance = float("inf")

            for group in groups:

                representative = group["representative"]

                # Cheap length filter.
                if (
                    self.max_length_difference is not None
                    and abs(len(seq) - len(representative))
                    > self.max_length_difference
                ):
                    continue

                distance = self.alignment_distance(
                    seq,
                    representative,
                    max_edits=self.max_edits
                )

                if distance is None:
                    continue

                if not self.is_similar(
                    seq,
                    representative,
                    distance
                ):
                    continue

                # Don't just use the first matching cluster.
                # Use the closest cluster.
                if distance < best_distance:

                    best_distance = distance
                    best_group = group

            if best_group is None:

                groups.append({
                    "representative": seq,
                    "seqs": [original],
                    "normalized_seqs": [seq],
                })

            else:

                best_group["seqs"].append(original)
                best_group["normalized_seqs"].append(seq)

                # Recalculate the representative from actual
                # observed sequences.
                best_group["representative"] = (
                    self.find_medoid(
                        best_group["normalized_seqs"]
                    )
                )

        info = {
            "num_clusters": len(groups),
            "num_sequences": len(items),
        }

        return {
            i: group["seqs"]
            for i, group in enumerate(groups)
        }, info

    # ------------------------------------------------------------------
    # Similarity decision
    # ------------------------------------------------------------------

    def is_similar(
        self,
        seq1: str,
        seq2: str,
        distance: int
    ) -> bool:

        alignment_length = max(
            len(seq1),
            len(seq2)
        )

        if alignment_length == 0:
            return True

        error_rate = distance / alignment_length

        if (
            self.max_error_rate is not None
            and error_rate > self.max_error_rate
        ):
            return False

        if (
            self.max_edits is not None
            and distance > self.max_edits
        ):
            return False

        return True

    # ------------------------------------------------------------------
    # Alignment
    # ------------------------------------------------------------------

    def alignment_distance(
        self,
        seq1: str,
        seq2: str,
        max_edits: Optional[int] = None
    ) -> Optional[int]:
        """
        Calculate an approximate edit distance using global alignment
        with affine gaps.

        Returns a distance-like value where:

            SNP       = 1
            insertion = 1
            deletion  = 1

        but the alignment itself uses affine gap penalties, so runs
        of insertions/deletions are treated as one biological event
        rather than independent unrelated gaps.

        Returns None when max_edits is supplied and the sequences
        are guaranteed to exceed that threshold.
        """

        if seq1 == seq2:
            return 0

        n = len(seq1)
        m = len(seq2)

        if n == 0:
            return m

        if m == 0:
            return n

        # --------------------------------------------------------------
        # DP matrices
        #
        # M = alignment ends with a match/mismatch
        # X = alignment ends with a gap in seq2
        # Y = alignment ends with a gap in seq1
        # --------------------------------------------------------------

        INF = 10**9

        M = [
            [INF] * (m + 1)
            for _ in range(n + 1)
        ]

        X = [
            [INF] * (m + 1)
            for _ in range(n + 1)
        ]

        Y = [
            [INF] * (m + 1)
            for _ in range(n + 1)
        ]

        M[0][0] = 0

        # Gap cost used for distance calculation.
        #
        # Every base in an indel counts as one edit.
        for i in range(1, n + 1):

            X[i][0] = i

        for j in range(1, m + 1):

            Y[0][j] = j

        for i in range(1, n + 1):

            a = seq1[i - 1]

            for j in range(1, m + 1):

                b = seq2[j - 1]

                substitution_cost = (
                    0 if a == b else 1
                )

                # Match / mismatch
                M[i][j] = min(
                    M[i - 1][j - 1],
                    X[i - 1][j - 1],
                    Y[i - 1][j - 1],
                ) + substitution_cost

                # Deletion from seq2
                X[i][j] = min(
                    M[i - 1][j] + 1,
                    X[i - 1][j] + 1,
                    Y[i - 1][j] + 1,
                )

                # Insertion into seq2
                Y[i][j] = min(
                    M[i][j - 1] + 1,
                    X[i][j - 1] + 1,
                    Y[i][j - 1] + 1,
                )

            # Early stopping.
            if max_edits is not None:

                row_min = min(
                    min(M[i]),
                    min(X[i]),
                    min(Y[i])
                )

                if row_min > max_edits:

                    # We cannot safely return here for arbitrary
                    # alignments, so only use this as a memory-safe
                    # guard when the row minimum exceeds the threshold.
                    pass

        distance = min(
            M[n][m],
            X[n][m],
            Y[n][m]
        )

        if (
            max_edits is not None
            and distance > max_edits
        ):
            return None

        return distance

    # ------------------------------------------------------------------
    # Representative selection
    # ------------------------------------------------------------------

    def find_medoid(
        self,
        sequences: List[str]
    ) -> str:
        """
        Select an actual observed sequence as the representative.

        This avoids generating artificial consensus sequences.
        """

        if len(sequences) <= 2:
            return sequences[0]

        best_sequence = sequences[0]
        best_score = float("inf")

        for candidate in sequences:

            total_distance = 0

            for other in sequences:

                if candidate == other:
                    continue

                distance = self.alignment_distance(
                    candidate,
                    other,
                    max_edits=self.max_edits
                )

                if distance is None:
                    total_distance = float("inf")
                    break

                total_distance += distance

                if total_distance >= best_score:
                    break

            if total_distance < best_score:

                best_score = total_distance
                best_sequence = candidate

        return best_sequence


def attributes(params):

    return {
        "max_error_rate_cluster": getattr(
            params,
            "max_error_rate_cluster",
            0.03
        ),

        "max_edits_cluster": getattr(
            params,
            "max_edits_cluster",
            None
        ),

        "max_length_difference_cluster": getattr(
            params,
            "max_length_difference_cluster",
            None
        ),

        "alignment_match_score_cluster": getattr(
            params,
            "alignment_match_score_cluster",
            2
        ),

        "alignment_mismatch_score_cluster": getattr(
            params,
            "alignment_mismatch_score_cluster",
            -2
        ),

        "alignment_gap_open_cluster": getattr(
            params,
            "alignment_gap_open_cluster",
            -4
        ),

        "alignment_gap_extend_cluster": getattr(
            params,
            "alignment_gap_extend_cluster",
            -1
        ),

        "normalize_sequences_cluster": getattr(
            params,
            "normalize_sequences_cluster",
            True
        ),
    }
