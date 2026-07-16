#!/usr/bin/env python
"""Debug script to trace clustering + recovery in detail"""

from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerGrouper
from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder
from dnabyte.params import Params

# Use actual church-encoded sequences
seqs_original = [
    'GGAAACACGCCTTCACCCAGAACAGGAACGAAGGAACGAAGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACCGAAGGCACCCTGGTACTT',
    'GGAAACACGCCTTCACCCAGAGAACAGGAACAACGGAACAACGGAACAACGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACGAAGGAACGAAGGCCGAAGGCACCCTGGTACTT',
]

print("="*80)
print("CLUSTERING")
print("="*80)

params = Params(primer_length=20)
grouper = PrimerGrouper()

print(f"\nInput: {len(seqs_original)} sequences of {len(seqs_original[0])}bp each")
for i, seq in enumerate(seqs_original):
    print(f"  [{i}] {seq[:30]}...{seq[-30:]}")

grouped = grouper.cluster(seqs_original)

print(f"\nGrouped into {len(grouped)} groups:")
for key, payloads in grouped.items():
    left, right, *rest = key
    hash_val = rest[0] if rest else "N/A"
    print(f"  Key: ({left[:10]}..., ...{right[-10:]}, {hash_val})")
    print(f"    Payloads: {len(payloads)} items")
    for i, payload in enumerate(payloads[:2]):  # Show first 2
        print(f"      [{i}] Length: {len(payload)}bp, Content: {payload[:40]}...")
    if len(payloads) > 2:
        print(f"      ... and {len(payloads)-2} more")

print("\n" + "="*80)
print("RECOVERY")
print("="*80)

params2 = Params(debruijn_kmer_size=21, debruijn_min_coverage=1)
recovery = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=1)

print(f"\nInput from clustering: {len(grouped)} groups")
result = recovery.recover(grouped)
consensus_seqs, stats = result

print(f"\nResult: {len(consensus_seqs)} consensus sequences")
for i, seq in enumerate(consensus_seqs[:3]):
    print(f"  [{i}] Length: {len(seq)}bp")
    print(f"      {seq[:50]}...{seq[-30:]}")

print(f"\nStats: {stats}")
