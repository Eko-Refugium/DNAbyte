


import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass
from collections import defaultdict
from dnabyte.sequence import SimulateSequencing

class KMER(SimulateSequencing):
    """
    KMER_k channel implementation based on the paper.
    
    Models strand-dependent insertion, deletion, and substitution errors
    with k-mer memory effects.
    """
    def simulate(self, data):
        """
        Simulate sequencing errors using the KMER channel.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        # Read parameters with safe defaults.
        k = getattr(self.params, 'kmer_k', 1)
        p_ins = getattr(self.params, 'kmer_p_ins', 0.01)
        p_del = getattr(self.params, 'kmer_p_del', 0.01)
        p_sub = getattr(self.params, 'kmer_p_sub', 0.01)
        seed = getattr(self.params, 'kmer_seed', None)

        custom_transition_probs = getattr(self.params, 'kmer_transition_probs', None)
        custom_substitution_probs = getattr(self.params, 'kmer_substitution_probs', None)

        # Build channel parameters.
        if custom_transition_probs is not None and custom_substitution_probs is not None:
            channel_params = ChannelParameters(
                k=k,
                transition_probs=custom_transition_probs,
                substitution_probs=custom_substitution_probs,
            )
        else:
            channel_params = ChannelParameters.create_default(
                k=k,
                p_ins=p_ins,
                p_del=p_del,
                p_sub=p_sub,
            )

        channel = KMERChannel(channel_params, seed=seed)

        sequenced_data = []
        event_lists = []
        error_counter = 0

        for sequence in data:
            seq_list = list(sequence)
            y, events = channel.transmit(seq_list)
            sequenced_data.append(''.join(y))
            event_lists.append(events)
            error_counter += sum(1 for ev in events if ev != Event.TRA)

        info = {
            "average_copy_number": 1.0,
            "number_of_sequencing_errors": error_counter,
            "error_dict": {"events": event_lists}
        }

        return sequenced_data, info


def check_parameter(parameter, default, min_val, max_val, inputparams, allow_none=False):
    if not hasattr(inputparams, parameter):
        return default

    parameter_value = inputparams.__dict__[parameter]

    if parameter_value is None:
        if allow_none:
            return None
        return default

    if min_val is not None and parameter_value < min_val:
        raise ValueError(
            f"{parameter} must be greater than or equal to {min_val}, got {parameter_value}"
        )

    if max_val is not None and parameter_value > max_val:
        raise ValueError(
            f"{parameter} must be less than or equal to {max_val}, got {parameter_value}"
        )

    return parameter_value


def attributes(params):
    kmer_k = check_parameter(
        parameter="kmer_k",
        default=1,
        min_val=1,
        max_val=100,
        inputparams=params,
    )

    kmer_p_ins = check_parameter(
        parameter="kmer_p_ins",
        default=0.01,
        min_val=0.0,
        max_val=1.0,
        inputparams=params,
    )

    kmer_p_del = check_parameter(
        parameter="kmer_p_del",
        default=0.01,
        min_val=0.0,
        max_val=1.0,
        inputparams=params,
    )

    kmer_p_sub = check_parameter(
        parameter="kmer_p_sub",
        default=0.01,
        min_val=0.0,
        max_val=1.0,
        inputparams=params,
    )

    if (kmer_p_ins + kmer_p_del + kmer_p_sub) > 1.0:
        raise ValueError("kmer_p_ins + kmer_p_del + kmer_p_sub must be <= 1")

    # Optional advanced parameters.
    kmer_seed = check_parameter(
        parameter="kmer_seed",
        default=None,
        min_val=None,
        max_val=None,
        inputparams=params,
        allow_none=True,
    )
    kmer_transition_probs = check_parameter(
        parameter="kmer_transition_probs",
        default=None,
        min_val=None,
        max_val=None,
        inputparams=params,
        allow_none=True,
    )
    kmer_substitution_probs = check_parameter(
        parameter="kmer_substitution_probs",
        default=None,
        min_val=None,
        max_val=None,
        inputparams=params,
        allow_none=True,
    )

    if (kmer_transition_probs is None) != (kmer_substitution_probs is None):
        raise ValueError(
            "Provide both kmer_transition_probs and kmer_substitution_probs, or neither."
        )

    return {
        "kmer_k": kmer_k,
        "kmer_p_ins": kmer_p_ins,
        "kmer_p_del": kmer_p_del,
        "kmer_p_sub": kmer_p_sub,
        "kmer_seed": kmer_seed,
        "kmer_transition_probs": kmer_transition_probs,
        "kmer_substitution_probs": kmer_substitution_probs,
    }

# DNA alphabet
SIGMA_DNA = ['A', 'C', 'G', 'T']

# Watson-Crick complement mapping
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Channel events
class Event:
    BEG = 'Beg'  # Beginning
    INS = 'Ins'  # Insertion
    DEL = 'Del'  # Deletion
    SUB = 'Sub'  # Substitution
    TRA = 'Tra'  # Transmission (no error)

@dataclass
class ChannelParameters:
    """Parameters for the KMER_k channel"""
    k: int  # Memory parameter (kmer length)
    # Transition probabilities: p(event | kmer, previous_event)
    # Format: {(kmer_tuple, prev_event): {event: probability}}
    transition_probs: Dict[Tuple[Tuple[str, ...], str], Dict[str, float]]
    # Substitution probabilities: p(symbol | original_symbol, Sub)
    substitution_probs: Dict[str, Dict[str, float]]
    
    @classmethod
    def create_default(cls, k: int = 1, 
                      p_ins: float = 0.01, 
                      p_del: float = 0.01, 
                      p_sub: float = 0.01):
        """Create default i.i.d. channel parameters"""
        p_tra = 1.0 - p_ins - p_del - p_sub
        
        # Default: uniform transition probabilities
        transition_probs = defaultdict(lambda: {
            Event.INS: p_ins,
            Event.DEL: p_del,
            Event.SUB: p_sub,
            Event.TRA: p_tra
        })
        
        # Default: uniform substitution (equally likely wrong symbols)
        substitution_probs = {
            symbol: {s: 1.0/3 if s != symbol else 0.0 
                    for s in SIGMA_DNA}
            for symbol in SIGMA_DNA
        }
        
        return cls(k=k, 
                  transition_probs=transition_probs,
                  substitution_probs=substitution_probs)


class KMERChannel:
    """
    KMER_k channel implementation based on the paper.
    
    Models strand-dependent insertion, deletion, and substitution errors
    with k-mer memory effects.
    """
    
    def __init__(self, params: ChannelParameters, seed: int = None):
        self.params = params
        self.k = params.k
        self.rng = np.random.RandomState(seed)
        
    def get_kmer(self, x: List[str], t: int) -> Tuple[str, ...]:
        """
        Extract k-mer at position t with symmetric window.
        
        Args:
            x: Input sequence
            t: Position (0-indexed)
        
        Returns:
            Tuple of symbols forming the k-mer
        """
        L = len(x)
        radius = self.k // 2
        
        # Determine actual window based on boundaries
        left = max(0, t - radius)
        right = min(L, t + radius + 1)
        
        # Ensure symmetric window when possible
        if t < radius:
            right = min(L, 2 * t + 1)
        elif t >= L - radius:
            left = max(0, 2 * t - L + 1)
            
        return tuple(x[left:right])
    
    def get_transition_prob(self, kmer: Tuple[str, ...], 
                           prev_event: str, 
                           event: str) -> float:
        """Get transition probability for given kmer and events"""
        key = (kmer, prev_event)
        try:
            # [] access triggers the defaultdict factory for missing keys;
            # for a plain dict this raises KeyError and falls through below.
            return self.params.transition_probs[key].get(event, 0.0)
        except KeyError:
            pass

        # Fallback: scan for any entry with a matching prev_event
        for (k, e), probs in self.params.transition_probs.items():
            if e == prev_event:
                return probs.get(event, 0.0)

        return 0.0
    
    def transmit_symbol(self, x: List[str], t: int, 
                       prev_event: str) -> Tuple[List[str], str, bool]:
        """
        Transmit one symbol through the channel.
        
        Args:
            x: Input sequence
            t: Current position
            prev_event: Previous channel event
        
        Returns:
            (output_symbols, new_event, symbol_consumed)
        """
        kmer = self.get_kmer(x, t)
        key = (kmer, prev_event)

        # Single lookup: [] triggers the defaultdict factory for unseen keys.
        # For a plain dict without the key, fall back to any entry whose
        # prev_event matches, guaranteeing a non-zero probability vector.
        try:
            prob_dict = self.params.transition_probs[key]
        except KeyError:
            prob_dict = next(
                (p for (_, e), p in self.params.transition_probs.items()
                 if e == prev_event),
                None
            )
            if prob_dict is None:
                raise ValueError(
                    f"No transition probabilities found for key={key}"
                )

        # Sample event
        events = [Event.INS, Event.DEL, Event.SUB, Event.TRA]
        probs = np.array([prob_dict.get(ev, 0.0) for ev in events])
        total = probs.sum()
        if total <= 0:
            raise ValueError(
                f"All transition probabilities sum to zero for key={key}: {prob_dict}"
            )
        probs = probs / total
        
        event = self.rng.choice(events, p=probs)
        
        if event == Event.INS:
            # Insert random symbol, don't consume input
            symbol = self.rng.choice(SIGMA_DNA)
            return [symbol], Event.INS, False
            
        elif event == Event.DEL:
            # Delete symbol (no output)
            return [], Event.DEL, True
            
        elif event == Event.SUB:
            # Substitute with another symbol
            orig = x[t]
            sub_probs = self.params.substitution_probs[orig]
            symbols = [s for s in SIGMA_DNA if s != orig]
            probs = [sub_probs[s] for s in symbols]
            probs = np.array(probs) / sum(probs)  # Normalize
            
            symbol = self.rng.choice(symbols, p=probs)
            return [symbol], Event.SUB, True
            
        else:  # Event.TRA
            # Transmit correctly
            return [x[t]], Event.TRA, True
    
    def transmit(self, x: List[str]) -> Tuple[List[str], List[str]]:
        """
        Transmit entire sequence through the KMER_k channel.
        
        Args:
            x: Input DNA sequence (list of symbols)
        
        Returns:
            (y, events): Output sequence and list of events
        """
        y = []
        events = []
        prev_event = Event.BEG
        
        t = 0
        while t < len(x):
            output, event, consumed = self.transmit_symbol(x, t, prev_event)
            
            y.extend(output)
            events.append(event)
            prev_event = event
            
            if consumed:
                t += 1
        
        return y, events


class SamplingKMERChannel:
    """
    Sampling KMER channel implementation.
    
    Models the full DNA storage channel including:
    - Random drawing with replacement
    - Forward/backward (reverse-complement) reads
    - KMER_k channel transmission
    - Random permutation
    """
    
    def __init__(self, 
                 M: int,  # Number of strands
                 L: int,  # Length of each strand
                 c: float,  # Coverage depth
                 p_rc: float,  # Reverse-complement probability
                 channel_params_forward: ChannelParameters,
                 channel_params_backward: ChannelParameters = None,
                 draw_probs: np.ndarray = None,
                 seed: int = None):
        
        self.M = M
        self.L = L
        self.c = c
        self.N = int(c * M)  # Total number of reads
        self.p_rc = p_rc
        
        self.rng = np.random.RandomState(seed)
        
        # Drawing probabilities (default: uniform)
        if draw_probs is None:
            self.draw_probs = np.ones(M) / M
        else:
            assert len(draw_probs) == M
            self.draw_probs = draw_probs / draw_probs.sum()
        
        # Forward and backward channels
        self.channel_forward = KMERChannel(channel_params_forward, seed)
        
        if channel_params_backward is None:
            channel_params_backward = channel_params_forward
        self.channel_backward = KMERChannel(channel_params_backward, seed)
    
    def reverse_complement(self, x: List[str]) -> List[str]:
        """Compute reverse complement of a DNA strand"""
        return [COMPLEMENT[s] for s in reversed(x)]
    
    def transmit(self, X: List[List[str]]) -> Tuple[List[List[str]], 
                                                     List[int], 
                                                     List[bool]]:
        """
        Transmit through the sampling KMER channel.
        
        Args:
            X: List of M input DNA strands
        
        Returns:
            (Y, indices, is_forward): 
                - Y: List of N output reads
                - indices: Original strand index for each read
                - is_forward: True if forward read, False if backward
        """
        assert len(X) == self.M
        
        # Step 1: Draw N strands with replacement
        drawn_indices = self.rng.choice(self.M, size=self.N, 
                                       p=self.draw_probs)
        
        # Step 2 & 3: Reverse-complement and transmit
        Z = []
        is_forward = []
        
        for idx in drawn_indices:
            # Decide forward or backward
            forward = self.rng.rand() > self.p_rc
            is_forward.append(forward)
            
            # Get strand
            x = X[idx]
            
            if forward:
                # Transmit through forward channel
                y, _ = self.channel_forward.transmit(x)
            else:
                # Reverse-complement then transmit through backward channel
                x_rc = self.reverse_complement(x)
                y, _ = self.channel_backward.transmit(x_rc)
            
            Z.append(y)
        
        # Step 4: Random permutation
        perm = self.rng.permutation(self.N)
        Y = [Z[i] for i in perm]
        indices = [drawn_indices[i] for i in perm]
        is_forward_perm = [is_forward[i] for i in perm]
        
        return Y, indices, is_forward_perm


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def compute_error_statistics(x: List[str], y: List[str]) -> Dict[str, float]:
    """
    Compute error statistics using Levenshtein distance alignment.
    
    Returns:
        Dictionary with insertion, deletion, substitution rates
    """
    # Simple alignment (for demonstration - use proper alignment for accuracy)
    from difflib import SequenceMatcher
    
    matcher = SequenceMatcher(None, x, y)
    
    insertions = 0
    deletions = 0
    substitutions = 0
    matches = 0
    
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'insert':
            insertions += j2 - j1
        elif tag == 'delete':
            deletions += i2 - i1
        elif tag == 'replace':
            substitutions += max(i2 - i1, j2 - j1)
        elif tag == 'equal':
            matches += i2 - i1
    
    total = len(x)
    
    return {
        'insertion_rate': insertions / total,
        'deletion_rate': deletions / total,
        'substitution_rate': substitutions / total,
        'error_rate': (insertions + deletions + substitutions) / total,
        'matches': matches
    }


def sequences_to_string(seqs: List[List[str]]) -> List[str]:
    """Convert list of symbol lists to strings"""
    return [''.join(seq) for seq in seqs]


def string_to_sequence(s: str) -> List[str]:
    """Convert string to list of symbols"""
    return list(s)


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("KMER CHANNEL MODEL IMPLEMENTATION")
    print("=" * 80)
    
    # Example 1: Simple KMER_k channel
    print("\n### Example 1: KMER_1 Channel (i.i.d. errors) ###\n")
    
    params = ChannelParameters.create_default(
        k=1, 
        p_ins=0.02, 
        p_del=0.02, 
        p_sub=0.01
    )
    
    channel = KMERChannel(params, seed=42)
    
    # Test sequence
    x = string_to_sequence("ACGTACGTACGT")
    print(f"Input:  {''.join(x)}")
    
    y, events = channel.transmit(x)
    print(f"Output: {''.join(y)}")
    print(f"Events: {' '.join(events)}")
    
    stats = compute_error_statistics(x, y)
    print(f"\nError Statistics:")
    for key, val in stats.items():
        print(f"  {key}: {val:.4f}")
    
    # Example 2: Sampling KMER Channel
    print("\n\n### Example 2: Sampling KMER Channel ###\n")
    
    M = 10  # Number of strands
    L = 20  # Length of each strand
    c = 3.0  # Coverage
    
    # Generate random input strands
    X = []
    for i in range(M):
        strand = [np.random.choice(SIGMA_DNA) for _ in range(L)]
        X.append(strand)
    
    print(f"System parameters: M={M}, L={L}, c={c}, N={int(c*M)}")
    print(f"\nInput strands:")
    for i, x in enumerate(X):
        print(f"  Strand {i}: {''.join(x)}")
    
    # Create sampling channel
    sampling_channel = SamplingKMERChannel(
        M=M, L=L, c=c, p_rc=0.5,
        channel_params_forward=params,
        seed=42
    )
    
    Y, indices, is_forward = sampling_channel.transmit(X)
    
    print(f"\nOutput reads (N={len(Y)}):")
    for i, (y, idx, fwd) in enumerate(zip(Y, indices, is_forward)):
        direction = "FWD" if fwd else "BWD"
        print(f"  Read {i} [{direction}, from strand {idx}]: {''.join(y)}")
    
    # Coverage statistics
    print(f"\nCoverage statistics:")
    unique, counts = np.unique(indices, return_counts=True)
    print(f"  Strands drawn: {len(unique)}/{M}")
    print(f"  Average copies per drawn strand: {counts.mean():.2f}")
    print(f"  Min/Max copies: {counts.min()}/{counts.max()}")
    
    print("\n" + "=" * 80)


# Additional utility for channel parameter estimation
class ChannelEstimator:
    """Estimate KMER_k channel parameters from experimental data"""
    
    @staticmethod
    def estimate_from_alignments(alignments: List[Tuple[List[str], List[str]]],
                                 k: int = 1) -> ChannelParameters:
        """
        Estimate channel parameters from aligned sequence pairs.
        
        Args:
            alignments: List of (input, output) sequence pairs
            k: Memory parameter
        
        Returns:
            Estimated ChannelParameters
        """
        # Count events for each (kmer, prev_event) context
        transition_counts = defaultdict(lambda: defaultdict(int))
        substitution_counts = defaultdict(lambda: defaultdict(int))
        
        # Process each alignment
        for x, y in alignments:
            # Here you would implement proper alignment and counting
            # This is a placeholder for the concept
            pass
        
        # Convert counts to probabilities
        transition_probs = {}
        for context, event_counts in transition_counts.items():
            total = sum(event_counts.values())
            transition_probs[context] = {
                event: count / total 
                for event, count in event_counts.items()
            }
        
        substitution_probs = {}
        for symbol, sub_counts in substitution_counts.items():
            total = sum(sub_counts.values())
            substitution_probs[symbol] = {
                s: count / total 
                for s, count in sub_counts.items()
            }
        
        return ChannelParameters(
            k=k,
            transition_probs=transition_probs,
            substitution_probs=substitution_probs
        )




